#CHANGED
#Code adjusted 06-20-2014 by L.Vanbrabant
#For now, when bootstrap=FALSE, the p-value is only computed for the E-bar-square.

##To do##
#- fix issue with weights=NULL
#- add equality constraints to null-distribution
#- add overall=FALSE option

##############################
## explaining the arguments ##
##############################
#model                : A description of the user-specified model. Typically, the model is described using the lm model syntax.
#data                 : Data frame containing the observed variables used in the model.
#ui                   : Matrix (or vector in case of one single restriction only) defining the left-hand side of the restriction, ui%*%beta >= ci, where beta is the parameter vector.
#meq                  : Integer number (default 0) giving the number of rows of ui that are used for equality restrictions instead of inequality restrictions.
#overall              : Currently not used
#pvalue               : If TRUE (default), a p-value is computed
#bootstrap            : If TRUE, the p-value is computed based on simulation. Otherwise, the p-value is computed based on the calculated weights.
#p.distr              : Assumed error-distribution (normal by default, "N") for computing a bootstrapped p-value. Two other options are the t-distribution ("T") and the chi^2-distributions ("Chi").
#df                   : Degrees of freedom, when p.distr="T" of p.distr="Chi".
#R                    : Integer; number of bootstrap draws. The default value is set to 99999.
#double.bootstrap     : If TRUE the genuine double bootstrap is used to compute an additional set of plug-in p-values for each bootstrap sample.
#double.bootstrap.R   : Integer; number of double bootstrap draws. The default value is set to 9999.
#R2                   : Computes the R-squared based on the constrained residuals.
#parallel             : The type of parallel operation to be used (if any). If missing, the default is set "no".
#ncpus                : Integer: number of processes to be used in parallel operation: typically one would chose this to the number of available cores.
#cl                   : An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the InformativeTesting call.
#seed                 : Seed value
#verbose              : Logical; if TRUE, information is shown at each bootstrap draw.
#...                  : Currently not used. 

source('my.quadprog.R')

CSI <- function(model, data, ui=NULL, meq=0, overall=TRUE,
                pvalue=TRUE, bootstrap=FALSE, p.distr="N", df=7, R=99999, 
                double.bootstrap=FALSE, double.bootstrap.R=9999, R2=FALSE,
                parallel=c("no", "multicore", "snow"), ncpus=1L, cl=NULL,                 
                seed=NULL, verbose=FALSE, ...) {
  
  #checks
  parallel <- tolower(parallel)
  stopifnot(parallel %in% c("no", "multicore", "snow"),
            p.distr %in% c("N", "T", "Chi"))
  
  if(!is.data.frame(data)) stop("\n The data should be a dataframe.") 
  if(is.null(ui)) stop("\n No constraints matrix has been specified.")  
  
  T.obs <- rep(as.numeric(NA), 6)
  p.value <- rep(as.numeric(NA), 6)
  Rboot.tot <- as.numeric(NA)
  
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
  
  #Weights FIXME!
  #if(is.missing(weights)) { weights <- NULL }
  fit.lm <- lm(as.formula(model), data, weights=NULL)
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
  #Overall test
#  if(overall) {
#    ui.eq <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1)))
#    optim.h0 <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(ui.eq), meq=nrow(ui.eq))
#  }else {  
#    optim.h0 <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(ui), meq=nrow(ui))
#  }
  
#  optim.h0$solution[abs(optim.h0$solution) < sqrt(.Machine$double.eps)] <- 0L  
#  par.h0 <- optim.h0$solution
#  RSS.h0 <- sum((Y - (X %*% par.h0))^2)
  
  par.h0 <- c(mean(predict(fit.lm)), rep(0, (p-1)))
  RSS.h0 <- sum((Y - mean(predict(fit.lm)))^2)
  
  #fit h1 
  optim.h1 <- my.solve.QP(Dmat=XX, dvec=Xy, Amat=t(ui), meq=meq)
    optim.h1$solution[abs(optim.h1$solution) < sqrt(.Machine$double.eps)] <- 0L  
  par.h1 <- optim.h1$solution
  RSS.h1 <- sum((Y - (X%*%par.h1))^2)
    
  #number of active order constraints
  iact <- optim.h1$iact 
  
  #fit h2
  RSS.h2 <- sum(resid(fit.lm)^2)
  par.h2 <- optim.h1$unconstrainted.solution
  
  #Transform RSS into LRT
  LRT.A <- 2*((mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h1))) -
              (mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h0))))
  
  LRT.B <- 2*((mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h2))) -
              (mllik <- -(n/2)*log(2*pi/n) - (n/2) - ((n/2)*log(RSS.h1))))

    
  #summary(fit.lm)$sigma^2
    #v = sum(n, -p)
    #s2 <- 1 / v * RSS.h2 
  
  s2 <- summary(fit.lm)$sigma^2
  df.error <- summary(fit.lm)$fstatistic[3]
  
    
  #compute observed Fbar values
  #T.obs[1] <- t(par.h1 - par.h0) %*% solve(vcov(fit.lm), par.h1 - par.h0)
  #T.obs[2] <- t(par.h2 - par.h1) %*% solve(vcov(fit.lm), par.h2 - par.h1)
  T.obs[1] <-  (RSS.h0 - RSS.h1) / s2
  T.obs[2] <-  (RSS.h1 - RSS.h2) / s2 
  
  #compute observed Ebar-square values
  #T1 <- t(par.h1 - par.h0) %*% solve(vcov(fit.lm)/summary(fit.lm)$sigma^2, par.h1 - par.h0)
  #T2 <- t(par.h2 - par.h1) %*% solve(vcov(fit.lm)/summary(fit.lm)$sigma^2, par.h2 - par.h1)
  #T.obs[3] <- T1/(df.error * s2 + T1)
  #T.obs[4] <- T2/(df.error * s2 + T2)
  
  T.obs[3] <-  (RSS.h0 - RSS.h1) / RSS.h0
  T.obs[4] <-  (RSS.h1 - RSS.h2) / RSS.h1

  T.obs[5] <- LRT.A
  T.obs[6] <- LRT.B
  
  #Fix negative and very small values to zero? 
  ind.zero <- which(T.obs < 1e-14) 
  T.obs <- replace(T.obs, ind.zero, 0)  
  
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
  if(pvalue) {
    if(bootstrap) {    
      T.boot <- matrix(as.numeric(NA), R, 6)
          
      fn <- function(b) { 
        if (verbose) cat("R =", b)
        if(!is.null(seed)) set.seed(seed + b)
        if(!exists(".Random.seed", envir = .GlobalEnv))
          runif(1)              
        RNGstate <- .Random.seed
        
        #bootstrapped p-value based on normal-, T-, or Chi-square distribution
        #Additional distributions can be added if needed.
        if(p.distr=="N")
          Yboot <- rnorm(n=n, 0, 1) 
        else if(p.distr=="T")
          Yboot <- rt(n=n, df=df)
        else if(p.distr=="Chi")
          Yboot <- rchisq(n=n, df=df)
        
        boot.data <- data.frame(Yboot, X[,-1])
        
        out <- CSI(model=formula(boot.data), data=boot.data, ui=ui, meq=meq, 
                   overall=overall, pvalue=FALSE, R2=R2, seed=NULL, 
                   verbose=verbose) 
        
        if(verbose) cat("  ... ... T.obs  = ", out$T.obs, "\n")
        
        out <- out$T.obs
        #out <- out$iact
        
        #standard double boostrap
        if(double.bootstrap) {
                  
          if (verbose) cat("   ... ... ... calibrating p.value - ")
          
          plugin.pvalue <- CSI(model, data=boot.data, meq=meq, ui=ui,
                               overall=TRUE, pvalue=TRUE, R=double.bootstrap.R,  
                               p.distr=p.distr, df=df, R2=R2,
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
    }else if(!bootstrap){
        cov <- vcov(fit.lm)
        
        #only inequality constraints 
        if (meq==0) {
          wt.bar <- ic.infer:::ic.weights(ui %*% cov %*% t(ui)) 
        }
        #equality and inequality constraints
        else if (meq > 0) {
          wt.bar <- ic.infer:::ic.weights(solve(solve(ui %*% cov %*% t(ui))[-(1:meq),-(1:meq)])) 
          
          stop("equality constraints not yet implemented, use bootstrap=TRUE")
        }
      
        ##Compute p-values for hypothesis test Type A, see Silvapulle and Sen, 2005, p99-100 or
        r=Matrix:::rankMatrix(ui)[1]
        q=(nrow(ui) - meq) 
        i=0:q
        wt.bar2=rev(wt.bar)
        
        #df2a=(n-p+q-i)/2
        df2a=rep(df.error, r+1)/2
        
        #less than the max. number of constraints
        #In case of the overall test
        if(overall) { 
          df1a <- (((length(par.h1)-1) - nrow(ui)):((length(par.h1)-1) - meq))/2  #CHECK IF CORRECT!
        }else{        
          df1a=(r-q+i)/2
                
          rank <- ic.infer:::RREF(t(ui), tol=sqrt(.Machine$double.eps))$rank
          if(rank == (p-1)) { 
  #          #max number of constraints
            df1a <- (0:(nrow(ui) - meq))/2 
          }
        }
        #These are adjusted functions of the pbetabar() function from the ic.infer package.
        pbar.A <- function(x, df1a, df2a, wt) {
            if (x <= 0) { 
              return(0)
            }
          zed <- df1a == 0
          cdf <- ifelse(any(zed), wt[zed], 0)
          
          cdf <- cdf + sum(pbeta(x, df1a[!zed], df2a[!zed])*wt[!zed])
          
          return(cdf)
        }
        
        Ebar.pA <- 1 - pbar.A(x=T.obs[3], df1a=df1a, df2a=df2a, wt=wt.bar2)
        
                
        #Hypothesis test Type B
        #df1b=i/2
        df1b <- meq:nrow(ui)/2
        df2b=(n-p)/2
                
        pbar.B <- function(x, df1b, df2b, wt) {
           if (x <= 0) {
              return(0)
          }
          zed <- df1b == 0
          cdf <- ifelse(any(zed), wt[zed], 0)
        
          cdf <- cdf + sum(pbeta(x, df1b[!zed], df2b)*wt[!zed])
          
          return(cdf)
        }
        
        Ebar.pB <- 1 - pbar.B(x=T.obs[4], df1b=df1b, df2b=df2b, wt=wt.bar)
        
        p.value[3:4] <- c(Ebar.pA, Ebar.pB)
        
      }      
            
    }
    
    #When a bootstrapped p-value is computed and parallel processing is used.
    if(bootstrap){
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
    
    #Assign names to the output
    names(T.obs) <- c("Fbar.A","Fbar.B","Ebar^2.A","Ebar^2.B","LRT.A","LRT.B")
    names(p.value) <- c("Fbar.A","Fbar.B","Ebar^2.A","Ebar^2.B","LRT.A","LRT.B")
  
  
  out <- list(T.obs=round(T.obs,4), iact=iact, p.value=round(p.value,4), 
              Rboot.tot=if(bootstrap) {Rboot.tot},
              weights=w.model, ui=ui, meq=meq, R2=if(R2) {Rsq},
              par.h0=par.h0,
              par.h1=par.h1, 
              par.h2=optim.h1$unconstrainted.solution)
  
  class(out) <- "CSI"
  
  return(out)
  
}



