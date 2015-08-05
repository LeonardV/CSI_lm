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
# object              : lm object
# Amat                : Matrix (or vector in case of one single restriction only) defining the left-hand side of the restriction, Amat%*%beta >= ci, where beta is the parameter vector.
# bvec                : Vector holding the values of ci (defaults to zero).
# meq                 : Integer number (default 0) giving the number of rows of Amat that are used for equality restrictions instead of inequality restrictions.
# pvalue              : If TRUE (default), a p-value is computed
# mix.weights         : The procedure of compuinge the p-value. If "none" (first approach), a bootstrapped p-value is computed. 
#                       If "boot" (second approach) the weights are computed based on a simulation procedure. 
#                       If "mvtnorm" (third approach), the weights are computed based on the multivariate normal probability distribution. 
# R                   : Integer; number of bootstrap draws. The default value is set to 9999.
# p.distr             : Assumed error-distribution (normal by default, "n") for computing a bootstrapped p-value. Two other options are the t-distribution ("t") and the chi^2-distributions ("chi").
# df                  : Degrees of freedom, when p.distr="t" of p.distr="chi".
# parallel            : The type of parallel operation to be used (if any). If missing, the default is set "no".
# ncpus               : Integer: number of processes to be used in parallel operation: typically one would chose this to the number of available cores.
# cl                  : An optional parallel or snow cluster for use if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the InformativeTesting call.
# seed                : Seed value
# verbose             : Logical; if TRUE, information is shown at each bootstrap draw.
# ...                 : Currently not used. 

source('my.quadprog.R')

csi.lm <- function (object, Amat = NULL, bvec = NULL, meq = 0, overall = FALSE,
                    pvalue = TRUE, mix.weights = c("mvtnorm", "boot", "none"), 
                    R = 9999, p.distr = c("n", "t", "chi"), df = 7, 
                    parallel = c("no", "multicore", "snow"), ncpus = 1L, 
                    cl = NULL, seed = NULL, verbose = FALSE, ...) 
{
  if (qr(Amat)$rank < nrow(Amat) && mix.weights != "none") {
    stop("Matrix Amat must have full row-rank. Set mixing weights to 'none'.")
  }
  parallel <- tolower(parallel)
  p.distr <- tolower(p.distr)
  mix.weights <- tolower(mix.weights)
  stopifnot(parallel %in% c("no", "multicore", "snow"), p.distr %in% 
              c("n", "t", "chi"), mix.weights %in% c("none", "boot", 
                                                     "mvtnorm"))
  if (is.null(Amat)) {
    stop("no constraints matrix has been specified.")
  }
  if (meq == nrow(Amat)) {
    stop("test not applicable with equality restrictions only.")
  }
  if (is.null(bvec)) {
    bvec <- rep(0, nrow(Amat))
  }
  if (!is.vector(bvec)) {
    stop("bvec must be a vector.")
  }
  p.distr <- match.arg(p.distr)
  mix.weights <- match.arg(mix.weights)
  stat <- vector("numeric", 2)
  p_value <- vector("numeric", 2)
  Rboot.tot <- as.numeric(NA)
  wt.bar <- vector("numeric", nrow(Amat) + 1)
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
  fit.lm <- object
  mfit <- fit.lm$model
  Y <- model.response(mfit)
  X <- model.matrix(fit.lm)[,,drop = FALSE]
  w <- NULL
  n = length(Y)
  b <- coef(fit.lm)
  p <- length(b)
  cov <- vcov(fit.lm)
  s2 <- summary(fit.lm)$sigma^2
  df.error <- summary(fit.lm)$fstatistic[[3]]
  if (!is.null(w)) {
    W <- diag(w.model)
    XX <- t(X) %*% W %*% X
    Xy <- t(X) %*% W %*% Y
  }
  else {
    XX <- crossprod(X)
    Xy <- t(X) %*% Y
  }
  
  if(overall) {
    RSS.h0 <- sum((Y - mean(predict(fit.lm)))^2)
    par.h0 <- rep(mean(predict(fit.lm)), p)  
  } else {
      out.eq <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat), bvec = bvec, 
                            meq = nrow(Amat))
      out.eq$solution[abs(out.eq$solution) < sqrt(.Machine$double.eps)] <- 0L
      par.h0 <- out.eq$solution
      RSS.h0 <- sum((Y - (X %*% par.h0))^2)
    }  
  out.ic <- my.solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat), 
                        bvec = bvec, meq = meq)
  out.ic$solution[abs(out.ic$solution) < sqrt(.Machine$double.eps)] <- 0L
  par.h1 <- out.ic$solution
  RSS.h1 <- sum((Y - (X %*% par.h1))^2)
  iact <- out.ic$iact
  par.h2 <- out.ic$unconstrained.solution
  RSS.h2 <- sum(resid(fit.lm)^2)
  stat[1] <- (RSS.h0 - RSS.h1)/s2
  stat[2] <- (RSS.h1 - RSS.h2)/s2
  ind.zero <- which(stat < 1e-14)
  stat <- replace(stat, ind.zero, 0)
  residuals <- Y - (X %*% par.h1)
  Rsq <- 1 - sum(residuals^2)/sum((Y - mean(Y))^2)
  if (is.null(weights(fit.lm)) && !attr(fit.lm$terms, "intercept")) {
    Rsq <- 1 - sum(residuals^2)/sum(Y^2)
  }
  
  if (pvalue) {
    if (mix.weights == "none") {
      Tboot <- matrix(as.numeric(NA), R, 2)
      fn <- function(b) {
        if (verbose) 
          cat("R =", b)
        if (!is.null(seed)) 
          set.seed(seed + b)
        if (!exists(".Random.seed", envir = .GlobalEnv)) 
          runif(1)
        RNGstate <- .Random.seed
        if (p.distr == "n") {
          Yboot <- rnorm(n, 0, 1)
        }
        else if (p.distr == "t") {
          Yboot <- rt(n, df = df)
        }
        else if (p.distr == "chi") {
          Yboot <- rchisq(n = n, df = df)
        }        
        X <- model.matrix(fit.lm)[,,drop=TRUE]
        object <- lm(Yboot ~ X -1)
        
        out <- csi.lm(object, Amat = Amat, bvec = bvec, meq = meq, 
                      overall = overall, mix.weights = "none", pvalue = FALSE, 
                      R = 0L, p.distr = p.distr, df = df,  
                      parallel = "no", ncpus = 1L, cl = NULL, 
                      seed = seed, verbose = verbose)
        if (verbose) 
          cat(" ...FbarA = ", format(out$stat[1], digits=4, nsmall=3), 
              "...FbarB = ", format(out$stat[2], digits=4, nsmall=3), "\n")
        out <- out$stat
        out
      }
    }
    else if (mix.weights == "mvtnorm" | mix.weights == "boot") {
      if (mix.weights == "boot") {
        #if (meq != 0L) 
        #  stop("not yet implemented, set mix.weights to \"mvtnorm\" or \"none\"")
        start.idx <- meq+2
        idx <- start.idx:ncol(Amat)# + (start.idx - 1))
        Amatw <- rbind(diag(p)[idx,])
        posPar <- matrix(as.numeric(NA), R, nrow(Amat)-meq)
        fn <- function(b) {
          if (verbose) 
            cat("R =", b, "\n")
          if (!is.null(seed)) 
            set.seed(seed + b)
          if (!exists(".Random.seed", envir = .GlobalEnv)) 
            runif(1)
          RNGstate <- .Random.seed
          if (p.distr == "n") {
            Yboot <- rnorm(n, 0, 1)
          }
          else if (p.distr == "t") {
            Yboot <- rt(n, df = df)
          }
          else if (p.distr == "chi") {
            Yboot <- rchisq(n, df = df)
          }
          
          X <- model.matrix(fit.lm)[,,drop = FALSE]
                    
          if (!is.null(w)) {
            W <- diag(w.model)
            XX <- t(X) %*% W %*% X
            Xy <- t(X) %*% W %*% Y
          }
          else {
            XX <- crossprod(X)
            Xy <- t(X) %*% Yboot
          }
          #start.idx <- min(sapply(1:nrow(Amat), function(x) which(Amat[x,]==1)))
          out.ic <- my.solve.QP(Dmat = XX, dvec = Xy, 
                                Amat = t(Amatw), meq=0L)
          out.ic$solution[abs(out.ic$solution) < sqrt(.Machine$double.eps)] <- 0L
          par <- out.ic$solution
          idx <- sapply(1:nrow(Amatw), function(x) which(Amatw[x,] == 1))
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
              cl <- parallel::makePSOCKcluster(rep("localhost", 
                                                   ncpus))
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
            posPar[b, 1:nrow(Amatw)] <- res[[b]]
          }
          else {
            error.idx <- c(error.idx, b)
          }
        }
        posPar <- sapply(1:R, function(x) sum(posPar[x,] > 0L))
        wt.bar <- sapply(0:nrow(Amatw), function(x) sum(posPar == x)/R)
          names(wt.bar) <- nrow(Amatw):0
      }
      else if (mix.weights == "mvtnorm") {
        if (meq == 0L) {
          wt.bar <- ic.infer:::ic.weights(Amat %*% cov %*% t(Amat))
        }
        else if (meq > 0) {
          wt.bar <- ic.infer:::ic.weights(solve(solve(Amat %*% 
                                                        cov %*% t(Amat))[-(1:meq), 
                                                                         -(1:meq)]))
        }
      }
      if(overall) {
        df1 <- ((length(par.h1) - 1) - nrow(Amat)):((length(par.h1) - 1) - meq)  
      } 
      else {
        df1 <- 0:(nrow(Amat) - meq)
      }
            
      pbar <- function(x, df1, df2, wt) {
        if (x <= 0) {
          return(0)
        }
        zed <- df1 == 0
        cdf <- ifelse(any(zed), wt[zed], 0)
        cdf <- cdf + sum(pf(x/df1[!zed], df1[!zed], 
                            df2) * wt[!zed])
        return(cdf)
      }
      
      p_value[1] <- 1 - pbar(x=stat[1], df1=df1, df2=df.error, wt = rev(wt.bar))
      
      df1 <- meq:nrow(Amat)
      p_value[2] <- 1 - pbar(x=stat[2], df1=df1, df2=df.error, wt = wt.bar)
    }
  }
  if (mix.weights == "none" & pvalue) {
    RR <- sum(R)
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
      if (have_mc) {
        parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
      }
      else if (have_snow) {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost", 
                                               ncpus))
          if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
            parallel::clusterSetRNGStream(cl)
          res <- parallel::parLapply(cl, seq_len(RR), 
                                     fn)
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
        Tboot[b, 1:2] <- res[[b]]
      }
      else {
        error.idx <- c(error.idx, b)
      }
    }
    na.boot.ind <- which(is.na(Tboot), arr.ind = TRUE)
    inf.boot.ind <- which(Tboot == Inf, arr.ind = TRUE)
    ind <- c(na.boot.ind[,1], inf.boot.ind[,1])
    ind.unique <- unique(ind)
    Rboot.tot <- (R - length(ind.unique))
    if (length(ind.unique) > 0) {
      Tboot <- Tboot[-ind.unique, ]
    }
    p_value[1] <- sum(Tboot[,1] > stat[1])/Rboot.tot
    p_value[2] <- sum(Tboot[,2] > stat[2])/Rboot.tot
  }
  names(p_value) <- names(stat) <- c("FbarA", "FbarB")
  out <- list(stat = stat, iact = iact, p_value = p_value, 
              Rboot.tot = if (mix.weights == "none") {Rboot.tot}, Amat = Amat, 
              meq = meq, wt.bar = if (mix.weights != "none") wt.bar, R2 = Rsq, 
              par.h0 = par.h0, par.h1 = par.h1, par.h2 = par.h2)
  class(out) <- "CSI"
    return(out)
}
