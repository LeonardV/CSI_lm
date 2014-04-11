#Code adjusted 04-11-2014 by L.Vanbrabant
#This code contains experimental parts. 

##to do##
# add weights

source('my.quadprog.R')

CSI <- function(model, data, weights, ui=NULL, meq=0, 
                pvalue=TRUE, bootstrap=FALSE, R=9999, 
                parallel=c("no", "multicore", "snow"), ncpus=1L, cl=NULL, 
                p.distr="N", df=7, R2=FALSE, 
                double.bootstrap=FALSE, double.bootstrap.R=9999,
                seed=NULL, verbose=FALSE) {
  
  #checks
  parallel <- tolower(parallel)
  stopifnot(parallel %in% c("no", "multicore", "snow"),
            p.distr %in% c("N", "T", "Chi"))#,
            #double.bootstrap %in% c("no", "standard"))
  
  if(!is.data.frame(data)) stop("\n The data should be a dataframe.") 
  if(is.null(ui)) stop("\n No constraints matrix has been specified.") 
    
  T.obs <- vector("numeric", 6)
  if(double.bootstrap == "standard") {
    F-bar.plugin.pvalue.A <- rep(as.numeric(NA), R)
    F-bar.plugin.pvalue.B <- rep(as.numeric(NA), R)
    E-bar.plugin.pvalue.A <- rep(as.numeric(NA), R)
    E-bar.plugin.pvalue.B <- rep(as.numeric(NA), R)
    LRT.A <- rep(as.numeric(NA), R)
    LRT.B <- rep(as.numeric(NA), R)
  }
      
  #prepare for parallel processing
  #only used when bootstrap=TRUE
  if(missing(parallel)) parallel <- "no"
  #if(missing(double.bootstrap)) double.bootstrap <- "no"
  
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if(parallel != "no" && ncpus > 1L) {
    if(parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if(parallel == "snow") 
      have_snow <- TRUE
    if(!have_mc && !have_snow) 
      ncpus <- 1L
  }
  
  if(missing(weights)) weights <- NULL
  fit.lm <- lm(as.formula(model), data, weights = weights)
  mfit <- fit.lm$model
  Y <- model.response(mfit)
  X <- model.matrix(fit.lm)[,,drop=FALSE]
  resid.lm <- fit.lm$residuals
  w.model <- fit.lm$weights
  n = length(Y) 
  p = length(coef(fit.lm))
  
  #weigths specified
  if(!is.null(w.model)) {
    W <- diag(w.model)
    XX <- t(X) %*% W %*% X
    Xy <- t(X) %*% W %*% Y
  }
  #no weights specified
  else {
    XX <- crossprod(X) 
    Xy <- t(X) %*% Y   
  }
    
  ##fit models
  
  #fit h0
  RSS.h0 <- sum((Y - mean(predict(fit.lm)))^2)
          
#  optim.h0 <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(ui), meq = 0L)
#  optim.h0$solution[abs(optim.h0$solution) < sqrt(.Machine$double.eps)] <- 0L  
#  par.h0 <- optim.h0$solution
#  RSS.h0 <- sum((Y - (X %*% par.h0))^2)
  
  
  #fit h1 
  optim.h1 <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(ui), meq = meq)
    optim.h1$solution[abs(optim.h1$solution) < sqrt(.Machine$double.eps)] <- 0L  
  par.h1 <- optim.h1$solution
  RSS.h1 <- sum((Y - (X %*% par.h1))^2)
  #number of active order constraints
  iact <- optim.h1$iact 
  
  #fit h2
  RSS.h2 <- sum(resid(fit.lm)^2)
  
  #Compute LRT
  LRT.A <- 2*((mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h1))) -
              (mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h0))))
  
  LRT.B <- 2*((mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h2))) -
              (mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h1))))

    
  #summary(fit.lm)$sigma^2
    #v = sum(n, -p)
    #s2 <- 1 / v * RSS.h2 
  
  s2 <- summary(fit.lm)$sigma^2
    
  #compute observed Fbar values
  T.obs[1] <-  (RSS.h0 - RSS.h1) / s2
  T.obs[2] <-  (RSS.h1 - RSS.h2) / s2 
  
  #compute observed Ebar-square values
  T.obs[3] <-  (RSS.h0 - RSS.h1) / RSS.h0
  T.obs[4] <-  (RSS.h1 - RSS.h2) / RSS.h1

  T.obs[5] <- LRT.A
  T.obs[6] <- LRT.B
  
  #Fix negative and very small values to zero? 
  #ind.zero <- which(T.obs < 1e-24) 
  #T.obs <- replace(T.obs, ind.zero, 0)  
  
  #Compute R-square for models with and without weights and with and without intercept
  Rsq <- NULL
  
  if(R2) {
  #constrained residuals
  residuals <- Y - (X %*% par.h1)
  #Code taken from package ic.infer
  Rsq <- 1 - sum(residuals^2) / sum((Y - mean(Y))^2)
  if(is.null(weights(fit.lm)) && !attr(fit.lm$terms, "intercept")) 
    Rsq <- 1 - sum(residuals^2) / sum(Y^2)
  else if(attr(fit.lm$terms, "intercept") && !is.null(weights(fit.lm))) 
    Rsq <- 1 - sum(weights(fit.lm) * residuals^2) / sum(weights(fit.lm) * 
                                                         (Y - weighted.mean(Y, w = w.model))^2)
  else if(!(attr(fit.lm$terms, "intercept") || is.null(weights(fit.lm)))) 
    Rsq <- 1 - sum(weights(fit.lm) * residuals^2) / sum(weights(fit.lm) * Y^2)
  }
  
  ##Compute p-value##
  p.value <- rep(as.numeric(NA), 6)
  Rboot.tot <- as.numeric(NA)
  
  if (pvalue) { 
    T.boot <- matrix(as.numeric(NA), R, 6)
        
    fn <- function(b) { 
      if (verbose) cat("R =", b)
      if(!is.null(seed)) set.seed(seed + b)
      if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)              
      RNGstate <- .Random.seed
      
      #bootstrapped p-value based on normal-, T-, or Chi-square distribution
      if(p.distr=="N")
        Yboot <- rnorm(n=n) 
      else if(p.distr=="T")
        Yboot <- rt(n=n, df=df)
      else if(p.distr=="Chi")
        Yboot <- rchisq(n=n, df=df)
      
      boot.data <- data.frame(Yboot, X[,-1])
        
      out <- Fbar(model=formula(boot.data), data=boot.data, weights=w.model, 
                  meq=meq, 
                  pvalue=FALSE, R=1L, ui=ui,  
                  p.distr=p.distr, df=df, R2=R2,
                  seed=NULL, verbose=verbose) 
      
      if(verbose) cat("  ... ... T.obs  = ", out$T.obs, "\n")
      
      out <- out$T.obs
      #out <- out$iact
      
      #standard double boostrap
      if(double.bootstrap) {
                
        if (verbose) cat("   ... ... ... calibrating p.value - ")
        
        plugin.pvalue <- Fbar(model, data=boot.data, weights=w.model, meq=meq, 
                              pvalue=TRUE, R=double.bootstrap.R, ui=ui, 
                              p.distr=p.distr, df=df, Rsquare=Rsquare,
                              seed=NULL, verbose=FALSE,
                              double.bootstrap=FALSE)$p.value
        
        if (verbose) cat(sprintf("%5.3f", plugin.pvalue), "\n\n")
        attr(out, "F-bar plugin.pvalue.A") <- plugin.pvalue[1]
        attr(out, "F-bar plugin.pvalue.B") <- plugin.pvalue[2]
        attr(out, "E^2-bar plugin.pvalue.A") <- plugin.pvalue[3]
        attr(out, "E^2-bar plugin.pvalue.B") <- plugin.pvalue[4]
        attr(out, "LRT.A") <- plugin.pvalue[5]
        attr(out, "LRT.B") <- plugin.pvalue[6]
        
      }
      
      out
    }
    
    RR <- sum(R)
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
      if (have_mc) {
        parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
      }
      else if (have_snow) {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus)) 
          if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
            parallel::clusterSetRNGStream(cl) 
          res <- parallel::parLapply(cl, seq_len(RR), fn)  
          parallel::stopCluster(cl) 
          res
        }
        else parallel::parLapply(cl, seq_len(RR), fn)
      }
    } 
    else lapply(seq_len(RR), fn)
        
    error.idx <- integer(0)
      for (b in seq_len(RR)) {
        if (!is.null(res[[b]])) {
          T.boot[b, 1:6] <- res[[b]]
          if (double.bootstrap) { 
            F-bar.plugin.pvalue.A[b] <- attr(res[[b]], "F-bar plugin.pvalue.A")
            F-bar.plugin.pvalue.B[b] <- attr(res[[b]], "F-bar plugin.pvalue.B")
            E-bar.plugin.pvalue.A[b] <- attr(res[[b]], "E^2-bar plugin.pvalue.A")
            E-bar.plugin.pvalue.B[b] <- attr(res[[b]], "E^2-bar plugin.pvalue.B")
            LRT.A[b] <- attr(res[[b]], "LRT.A")
            LRT.B[b] <- attr(res[[b]], "LRT.B")
          } 
        } else {
          error.idx <- c(error.idx, b)
        }
      }
    
    if(double.bootstrap) {  
      attr(p.value, "plugin.pvalue.A") <- F-bar.plugin.pvalue.A
      attr(p.value, "plugin.pvalue.B") <- F-bar.plugin.pvalue.B
      attr(p.value, "plugin.pvalue.B") <- E-bar.plugin.pvalue.A
      attr(p.value, "plugin.pvalue.B") <- E-bar.plugin.pvalue.B
      attr(p.value, "plugin.pvalue.A") <- LRT.A
      attr(p.value, "plugin.pvalue.B") <- LRT.B
    }
    #remove NA values / Inf values
    na.boot.ind   <- which(is.na(T.boot), arr.ind = TRUE)
    inf.boot.ind  <- which(T.boot == Inf, arr.ind = TRUE)
      ind <- c(na.boot.ind[,1], inf.boot.ind[,1])
    
    ind.unique <- unique(ind)
    Rboot.tot <- (R - length(ind.unique))
    
    if(length(ind.unique) > 0) 
      T.boot <- T.boot[-ind.unique,]
    
    p.value[1]  <- sum(T.boot[,1] > T.obs[1]) / Rboot.tot
    p.value[2]  <- sum(T.boot[,2] > T.obs[2]) / Rboot.tot
    p.value[3]  <- sum(T.boot[,3] > T.obs[3]) / Rboot.tot
    p.value[4]  <- sum(T.boot[,4] > T.obs[4]) / Rboot.tot
    p.value[5]  <- sum(T.boot[,5] > T.obs[5]) / Rboot.tot
    p.value[6]  <- sum(T.boot[,6] > T.obs[6]) / Rboot.tot
  }

  names(T.obs) <- c("Fbar.A","Fbar.B","Ebar^2.A","Ebar^2.B","LRT.A","LRT.B")
  names(p.value) <- c("Fbar.A","Fbar.B","Ebar^2.A","Ebar^2.B","LRT.A","LRT.B")


  out <- list(T.obs=T.obs, iact=iact, p.value=round(p.value,4), 
              Rboot.tot=Rboot.tot,
              weights=w.model, ui=ui, meq=meq, R2=Rsq,
              par.h0=mean(predict(fit.lm)),
              par.h1=par.h1, 
              par.h2=optim.h1$unconstrainted.solution)
  
  class(out) <- "Fbar"
  
  return(out)
  
}



