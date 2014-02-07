simulate.sd <- function(R2, n, sigma=diag(3), sims=100, start=1, seed=3013073, tol=0.01) {
  
  set.seed(seed)
  sd.value <- start
  results <- rep(as.numeric(NA), sims)
  rsquare <- 1:sims
  
  for (i in 1:sims) {
    
    iter <- 0 
    
    while (rsquare[i] > (R2 + tol) | rsquare[i] < (R2 - tol)) {
      sd.value <- sd.value + 0.001
      rsquare[i] <- simulate.sd.iter(sd.value, n, sigma)
      iter <- iter + 1
      if (iter > 1000) { break }
    }
    results[i] <- sd.value  
    sd.value <- 1
    
    cat("sim =", i, "...", "sd =", results[i], "...", "rsq =", rsquare[i], "\n")
  }
  cbind(results, rsquare)
}


simulate.sd.iter <- function(sd.value, n, sigma) {  
    
  X <- matrix(mvtnorm:::rmvnorm(n * 1, sigma = sigma), nrow = n)
  beta <- seq(0.1, 0.1*ncol(X), 0.1)
  Y <- rnorm(n, (0 + X%*%beta), sd.value)
  
  simdata <- data.frame(Y,X)
    
    return(summary(lm(Y ~ ., data=simdata))$r.squared)
}

