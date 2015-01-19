# Code adjusted 09-19-2014 by L.Vanbrabant
# Parts of the code are taken from the R package ic.infer (Grömping, 2010)


## to do ##
# no intercept model
# add confidence interval constrained estimates
# add residual bootstrap
# compute bootstrapped std. errors

##############################
## explaining the arguments ##
##############################
# model               : A description of the user-specified model. Typically, the model is described using the lm model syntax.
# data                : Data frame containing the observed variables used in the model.
# ui                  : Matrix (or vector in case of one single restriction only) defining the left-hand side of the restriction, ui%*%beta >= ci, where beta is the parameter vector.
# meq                 : Integer number (default 0) giving the number of rows of ui that are used for equality restrictions instead of inequality restrictions.
# bvec                : Vector holding the values of c (defaults to zero).
# pvalue              : If TRUE (default), a p-value is computed
# mix.weights         : The procedure of compuinge the p-value. If "none" (first approach), the p-value is computed directly without first calculating the mixing weights. 
#                       If "boot" (second approach) the weights are computed based on a simulation procedure. 
#                       If "mvtnorm" (third approach), the weights are computed based on the multivariate normal probability distribution. 
# R                   : Integer; number of bootstrap draws. The default value is set to 99999.
# double.bootstrap    : If a double bootstrap procedure should be used. (can be savely ignored in linear models)
# double.bootstrap.R  : Integer; number of double bootstrap draws. The default value is set to 9999
# p.distr             : Assumed error-distribution (normal by default, "N") for computing a bootstrapped p-value. Two other options are the t-distribution ("T") and the chi^2-distributions ("Chi").
# df                  : Degrees of freedom, when p.distr="T" of p.distr="Chi".
# R2                  : Computes the R-squared based on the constrained residuals.
# parallel            : The type of parallel operation to be used (if any). If missing, the default is set "no".
# ncpus               : Integer: number of processes to be used in parallel operation: typically one would chose this to the number of available cores.
# cl                  : An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the InformativeTesting call.
# seed                : Seed value
# verbose             : Logical; if TRUE, information is shown at each bootstrap draw.
# ...                 : Currently not used. 

source('my.quadprog.R')

csi.lm <- function(model, data, ui = NULL, bvec = NULL, meq = 0, pvalue = TRUE, 
                   mix.weights = c("mvtnorm", "none", "boot"), R=9999, 
                   double.bootstrap = FALSE, double.bootstrap.R = 9999, 
                   p.distr = c("n", "t", "chi"), df = 7, R2 = TRUE,
                   parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL,                 
                   seed = NULL, verbose = FALSE, ...) {
    
  # checks
  if (qr(ui)$rank < nrow(ui) && mix.weights == "mvtnorm") {
    stop("Matrix ui must have full row-rank.")
  }  
  
  parallel <- tolower(parallel)
  p.distr <- tolower(p.distr)
  mix.weights <- tolower(mix.weights)
  
  stopifnot(parallel %in% c("no", "multicore", "snow"),
            p.distr %in% c("n", "t", "chi"),
            mix.weights %in% c("none", "boot", "mvtnorm"))          
    
  #if (missing(p.distr)) { p.distr <- "n" }
  #if (missing(mix.weights)) { mix.weights <- "boot" }
  #if (missing(parallel)) { parallel <- "no" }
  if (!is.data.frame(data)) { stop("the data should be a dataframe.") }
  if (is.null(ui)) { stop("no constraints matrix has been specified.") }
  if (meq == nrow(ui)) { stop("test not applicable with equality restrictions only.") }
  if (is.null(bvec)) { bvec <- rep(0, nrow(ui)) }
  if (!is.vector(bvec)) { stop("bvec must be a vector.") }
  
  p.distr <- match.arg(p.distr)
  mix.weights <- match.arg(mix.weights)
    
  T.obs <- vector("numeric", 2)
  p.value <- vector("numeric", 2)
  Rboot.tot <- as.numeric(NA)
  wt.bar <- vector("numeric", nrow(ui) + 1)
  
  # prepare for parallel processing
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
  }
  
  # fit model
  fit.lm <- lm(as.formula(model), data, ...)
  # model
  mfit <- fit.lm$model
  # response variable
  Y <- model.response(mfit)
  # model.matrix
  X <- model.matrix(fit.lm)[,,drop=FALSE]
  # weights
  w <- NULL
  n = length(Y) 
  b <- coef(fit.lm)
  p <- length(b)
  cov <- vcov(fit.lm)
  s2 <- summary(fit.lm)$sigma^2
  df.error <- as.numeric(summary(fit.lm)$fstatistic[3])
    
  #weigths specified
  if(!is.null(w)) {
    W <- diag(w.model)
    XX <- t(X) %*% W %*% X
    Xy <- t(X) %*% W %*% Y
    #no weights specified
  } else {
    XX <- crossprod(X) 
    Xy <- t(X) %*% Y   
  }
    
  ### fit models ###
  if (attr(fit.lm$terms, "intercept")) {
    ui.eq <- cbind(rep(0, (p-1)), diag(rep(1, p-1)))
  } else {
    ui.eq <- cbind(diag(rep(1, p)))
  }  
  RSS.h0 <- sum((Y - mean(predict(fit.lm)))^2)
  par.h0 <- rep(mean(predict(fit.lm)), p)
  # full model RSS (possibly constrained) 
  out.ic <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(ui), 
                        bvec = bvec, meq = meq)
  out.ic$solution[abs(out.ic$solution) < sqrt(.Machine$double.eps)] <- 0L  
  par.h1 <- out.ic$solution
  RSS.h1 <- sum((Y - (X%*%par.h1))^2)
  # number of active order constraints
  iact <- out.ic$iact 
  # completely unconstrained RSS (not always the same...)
  par.h2 <- out.ic$unconstrained.solution
  RSS.h2 <- sum(resid(fit.lm)^2)
    
  # compute the F-bar test-statistics
  # Type A
  T.obs[1] <- (RSS.h0 - RSS.h1) / s2
  # Type B
  T.obs[2] <- (RSS.h1 - RSS.h2) / s2
  # fix negative and very small values to zero? 
  ind.zero <- which(T.obs < 1e-14) 
  T.obs <- replace(T.obs, ind.zero, 0)  
  
  # compute R-square for models with and without weights and with and without 
  # intercept
  Rsq <- NULL
  if (R2) {
    # constrained residuals
    residuals <- Y - (X %*% par.h1)
    # code taken from package ic.infer
    Rsq <- 1 - sum(residuals^2) / sum((Y - mean(Y))^2)
    if (is.null(weights(fit.lm)) && !attr(fit.lm$terms, "intercept")) 
      Rsq <- 1 - sum(residuals^2) / sum(Y^2)
    #    else if (attr(fit.lm$terms, "intercept") && !is.null(weights(fit.lm))) 
    #      Rsq <- 1 - sum(weights(fit.lm) * residuals^2) / 
    #      sum(weights(fit.lm) * (Y - weighted.mean(Y, w = w.model))^2)
    #    else if (!(attr(fit.lm$terms, "intercept") || is.null(weights(fit.lm)))) 
    #      Rsq <- 1 - sum(weights(fit.lm) * residuals^2) / sum(weights(fit.lm) * Y^2)
  }
  
  # compute p-value based on bootstrapping or mixing weights  
  if (pvalue) {
    # bootstrapped p value
    if (mix.weights == "none") {    
      T.boot <- matrix(as.numeric(NA), R, 2)
      
      fn <- function(b) { 
        if (verbose) cat("R =", b)
        if (!is.null(seed)) set.seed(seed + b)
        if (!exists(".Random.seed", envir = .GlobalEnv))
          runif (1)              
        RNGstate <- .Random.seed
        
        # bootstrapped p-value based on normal-, T-, or Chi-square distribution
        # Additional distributions can be added if needed.
        # The null-distribution does not depend on mu and sigma. So any value
        # can be used as long as the same values are used over all bootstrap runs.
        if (p.distr == "n") {
          Yboot <- rnorm(n, 0, 1) 
        } else if (p.distr == "t") {
          Yboot <- rt(n, df = df)
        } else if (p.distr == "chi") {
          Yboot <- rchisq(n=n, df=df)
        }
        boot.data <- data.frame(Yboot, X[,-1])
        
        out <- csi.lm(model = formula(boot.data), data = boot.data, ui = ui, 
                      bvec = bvec, meq = meq, mix.weights = "none", 
                      pvalue = FALSE, R = 0L, 
                      p.distr = p.distr, df = df, R2 = R2,
                      parallel = "no", ncpus = 1L, cl = NULL,
                      seed = seed, verbose = verbose) 
        
        if (verbose) cat("  ... ... T.obs  = ", out$T.obs, "\n")
        
        out <- out$T.obs
        
        out 
      }
    } else if (mix.weights == "mvtnorm" | mix.weights == "boot") {
      if (mix.weights == "boot") {
        if (meq != 0L) stop("not yet implemented. set mix.weights = mvtnorm")
        posPar <- matrix(as.numeric(NA), R, nrow(ui))
        
        fn <- function(b) { 
          if (verbose) cat("R =", b, "\n")
          if (!is.null(seed)) set.seed(seed + b)
          if (!exists(".Random.seed", envir = .GlobalEnv))
            runif (1)              
          RNGstate <- .Random.seed
          
          # bootstrapped weights based on normal-, T-, or Chi-square distribution
          # Additional distributions can be added if needed.
          # The null-distribution does not depend on mu and sigma. So any value
          # can be used as long as the same values are used over all bootstrap runs.
          if (p.distr=="n") {
            Yboot <- rnorm(n, 0, 1) 
          } else if (p.distr == "t") {
            Yboot <- rt(n, df = df)
          } else if (p.distr == "chi") {
            Yboot <- rchisq(n, df = df)
          }
          X <- model.matrix(fit.lm)[,,drop=FALSE]
          #weigths specified
          if(!is.null(w)) {
            W <- diag(w.model)
            XX <- t(X) %*% W %*% X
            Xy <- t(X) %*% W %*% Y
            #no weights specified
          } else {
            XX <- crossprod(X) 
            Xy <- t(X) %*% Yboot   
          }
          start.idx <- min(sapply(1:nrow(ui), function(x) which(ui[x,] == 1)))
          idx <- start.idx:(nrow(ui)+(start.idx-1))
          uiw <- rbind(diag(p)[idx,])
          out.ic <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(uiw), meq = 0L)
          out.ic$solution[abs(out.ic$solution) < sqrt(.Machine$double.eps)] <- 0L  
          par <- out.ic$solution
          idx <- sapply(1:nrow(ui), function(x) which(uiw[x,] == 1))
          out <- par[idx]
          
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
        for (b in seq_len(R)) {
          if (!is.null(res[[b]])) {
            posPar[b, 1:nrow(ui)] <- res[[b]]
          } else {
            error.idx <- c(error.idx, b)
          }
        }
        posPar <- sapply(1:R, function(x) sum(posPar[x,] > 0L))
        wt.bar <- sapply(0:nrow(ui), function(x) sum(posPar == x) / R)
        #wt.bar <- rev(wt.bar)
          names(wt.bar) <- nrow(ui):0        
      }
      else if (mix.weights == "mvtnorm") {
        # only inequality constraints 
        if (meq == 0L) {
          wt.bar <- ic.infer:::ic.weights(ui %*% cov %*% t(ui))
          # equality and inequality constraints
        } else if (meq > 0) {
          wt.bar <- ic.infer:::ic.weights(solve(solve(ui %*% cov %*% t(ui))[-(1:meq),-(1:meq)]))
        }
      }
      
      #r = qr(ui.eq)$rank
      #q = qr(ui)$rank - meq 
      #i = 0:q
      #df1a = r-q+i
      df1a <- ((length(par.h1)-1) - nrow(ui)):((length(par.h1)-1) - meq)
      df1b <- meq:nrow(ui)
      
      ## these are adjusted functions of the pbetabar() function from the ##
      ## ic.infer package ##
      
      #Hypothesis test Type A
      pbarA <- function(x, df1a, df2a, wt) {
        if (x <= 0) { 
          return(0)
        }
        zed <- df1a == 0
        cdf <- ifelse(any(zed), wt[zed], 0)
        cdf <- cdf + sum(pf(x/df1a[!zed], df1a[!zed], df2a)*wt[!zed])
        
        return(cdf)
      }
      Fbar.pA <- 1 - pbarA(x=T.obs[1], df1a=df1a, df2a=df.error, wt=rev(wt.bar))
      
      #Hypothesis test Type B
      pbarB <- function(x, df1b, df2b, wt) {
        if (x <= 0) {
          return(0)
        }
        zed <- df1b == 0
        cdf <- ifelse(any(zed), wt[zed], 0)
        cdf <- cdf + sum(pf(x/df1b[!zed], df1b[!zed], df2b)*wt[!zed])
        
        return(cdf)
      }
      Fbar.pB <- 1 - pbarB(x=T.obs[2], df1b=df1b, df2b=df.error, wt=wt.bar)
      
      p.value[1:2] <- c(Fbar.pA, Fbar.pB)
    }      
  }
  
  # when a bootstrapped p-value is computed and parallel processing is used.
  if (mix.weights == "none" & pvalue) {
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
        T.boot[b, 1:2] <- res[[b]]
      } else {
        error.idx <- c(error.idx, b)
      }
    }
    # remove NA values / Inf values
    na.boot.ind   <- which(is.na(T.boot), arr.ind = TRUE)
    inf.boot.ind  <- which(T.boot == Inf, arr.ind = TRUE)
    ind <- c(na.boot.ind[,1], inf.boot.ind[,1])
    
    ind.unique <- unique(ind)
    Rboot.tot <- (R - length(ind.unique))
    
    if (length(ind.unique) > 0) { 
      T.boot <- T.boot[-ind.unique,]
    }
    
    #compue bootstrap p-value
    p.value[1]  <- sum(T.boot[,1] > T.obs[1]) / Rboot.tot
    p.value[2]  <- sum(T.boot[,2] > T.obs[2]) / Rboot.tot
  }
  names(p.value) <- names(T.obs) <- c("Fbar.A","Fbar.B")
  
  out <- list(T.obs = T.obs, iact = iact, p.value = p.value, 
              Rboot.tot = if (mix.weights == "none") { Rboot.tot },
              ui = ui, meq = meq,
              wt.bar = wt.bar,
              R2 = if (R2) {Rsq},
              par.h0 = par.h0,
              par.h1 = par.h1,
              par.h2 = par.h2)
  
  class(out) <- "CSI"
  
  return(out)
  
}
