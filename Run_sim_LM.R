run_sim <- function(p=p, n=n, R2=R2, cov.x=0, model, ui, meq=0, sims=10000, 
                    pwr=0.80, R=10000, 
                    parallel="multicore", ncpus, cl=NULL,
                    method=c("weights", "bootstrap", "robust"),
                    double.bootstrap.R=double.bootstrap.R,
                    verbose=FALSE, seed=NULL) {
  
  #create data with pre-defined R^2
  data.lm <- function(n=n, X=X, beta=beta, var.epsilon=var.epsilon, seed=NULL) {
    
    if(!is.null(seed)) {
      set.seed(seed)
    } 
    epsilon <- rnorm(n, sd=sqrt(var.epsilon))
    y <- X%*%beta + epsilon
    
    data <- data.frame(y=y, x=X)
    
      return(data)
  }  
       
  #power function
  getpower <- function() {
    
    if(!is.null(seed)) {
      set.seed(seed)
    } 
    
    F.sigA <- F.sigB <- rep(NA, sims) ; E.sigA <- E.sigB <- rep(NA, sims)
    F.sigAB <- E.sigAB <- rep(NA, sims)
    T.obs <- pvalue <- matrix(as.numeric(NA), sims, 4)
    
    #covariance matrix between independent variables
    sigma <- matrix(rep(cov.x, p*p), p, p)
      diag(sigma) <- 1
    #data and parameters are fixed
    X <- mvtnorm:::rmvnorm(n, sigma = sigma)
    beta <- seq(0.1, 0.1*ncol(X), 0.1)
    #compute residual variance according to a pre-defined R^2
    var.epsilon <- (t(beta) %*% sigma %*% beta)  * (1 - R2) / R2
        
    for(i in 1:sims) {
      #simulate data with pre-defined R^2   
      sim.data <- data.lm(n=n, X=X, beta=beta, var.epsilon=var.epsilon, 
                          seed=seed+i)
      
      #fit model
      fit.lm <- lm(y ~ ., data=sim.data)
      data <- fit.lm$model
      
      #methods for calculating p-values
      
      #robust estimations (currently not used)
      if(method == "robust") {  
        fit <- rFbar(model, data, ui=ui, meq=meq, pvalue=TRUE, R=R,
                     cl=cl, parallel=parallel, ncpus=ncpus, 
                     scale.est="proposal 2",
                     verbose=verbose, psi="psi.bisquare", c=4.685)
      #bootstrap approach  
      } else if(method == "bootstrap"){
        fit <- Fbar(model, data, R=R, double.bootstrap.R=double.bootstrap.R,
                    parallel=parallel, ncpus=ncpus, cl=cl,
                    pvalue=TRUE, verbose=verbose, meq=meq, 
                    ui=ui, p.distr="N", df=4)
      #compute weights  
      } else if(method == "weights"){
        #Type A
        r=p 
        q=nrow(ui)
        x=(p+1) #include intercept
        idx=q:0
        N2=nrow(fit.lm$model)
        nu=fit.lm$df.residual
        df1a=r-q+idx
        df2a=N2-x+q-idx
        
        #Type B
        df1b=idx
        df2b=N2-p 
        
        #model implied covariance matrix
        cov <- vcov(fit.lm)
        
        #only inequality constraints 
        if (meq==0) {
          wt.bar <- ic.infer:::ic.weights(ui %*% cov %*% t(ui)) 
        }
        #equality and inequality constraints
        else if (meq > 0) {
          wt.bar <- ic.weights(solve(solve(uiw %*% cov %*% t(uiw))[-(1:meq), 
                                                                   -(1:meq)])) 
        }
        
        #compute observed F-bar and E-bar-square
        fit <- Fbar(model='y ~ .', data=data, R=0, 
                    parallel="multicore", ncpus=1,
                    pvalue=FALSE, verbose=TRUE, meq=0, 
                    ui=ui, p.distr="N", df=4)
        
        Fbar.A <- fit$T.obs[1] ; Fbar.B <- fit$T.obs[2]
        Ebar.A <- fit$T.obs[3] ; Ebar.B <- fit$T.obs[4]
        
        #################################################################################
        ##Compute p-values for hypothesis test Type A, see Silvapulle and Sen, 2005, p99-100
        #These are adjusted functions of the pbetabar() function from the ic.infer package.
        pbar.A <- function(x1, x2, df1a, df2a, nu, wt) {
          if (x1 <= 0) { 
            cdf <- list(cdf.ebarA=0L, cdf.fbarA=0L)
            return(cdf)
          }
          zed <- df1a == 0
          cdf <- ifelse(any(zed), wt[zed], 0)
          
          cdf.fbarA <- cdf + sum(pf(x1/df1a[!zed], df1a[!zed], nu)*wt[!zed])
          cdf.ebarA <- cdf + sum(pbeta(x2, df1a[!zed]/2, df2a[!zed]/2)*wt[!zed])
          
          cdf <- list(cdf.ebarA=cdf.ebarA, cdf.fbarA=cdf.fbarA)
          
          return(cdf)
        }
        
        out.A <- pbar.A(Fbar.A, Ebar.A, df1a, df2a, nu, wt.bar)
        
        Fbar.pA <- 1-out.A$cdf.fbarA
        Ebar.pA <- 1-out.A$cdf.ebarA
        
        #Hypothesis test Type B
        pbar.B <- function(x1, x2, df1b, df2b, nu, wt) {
          if (x1 <= 0) {
            cdf <- list(cdf.ebarB=0L, cdf.fbarB=0L)
            return(cdf)
          }
          zed <- df1b == 0
          cdf <- ifelse(any(zed), rev(wt)[zed], 0)
          
          cdf.fbarB <- cdf + sum(pf(x1/df1b[!zed], df1b[!zed], nu)*rev(wt)[!zed])
          cdf.ebarB <- cdf + sum(pbeta(x2, df1b[!zed]/2, df2b/2)*rev(wt)[!zed])
          
          cdf <- list(cdf.ebarB=cdf.ebarB, cdf.fbarB=cdf.fbarB)
          
          return(cdf)
        }
        
        out.B <- pbar.B(Fbar.B, Ebar.B, df1b, df2b, nu, wt.bar)
        Fbar.pB <- 1-out.B$cdf.fbarB
        Ebar.pB <- 1-out.B$cdf.ebarB
        
        fit$p.value <- attr(fit, "p.value") <- c(Fbar.pA, Fbar.pB, Ebar.pA, Ebar.pB)
        fit$T.obs <- attr(fit, "T.obs") <- c(Fbar.A, Fbar.B, Ebar.A, Ebar.B)
      }
      
      
      ##########################################################################
        pvalue[i,1:4] <- fit$p.value[1:4] ; T.obs[i,1:4] <- fit$T.obs[1:4]
        
        F.sigA[i] <- ifelse(pvalue[i,1] <= .05, 1, 0)
        F.sigB[i] <- ifelse(pvalue[i,2] <= .05, 1, 0)
        E.sigA[i] <- ifelse(pvalue[i,3] <= .05, 1, 0)
        E.sigB[i] <- ifelse(pvalue[i,4] <= .05, 1, 0)
        
        F.sigAB[i] <- ifelse((pvalue[i,1] <= 0.05 && pvalue[i,2] > 0.05), 1, 0) 
        E.sigAB[i] <- ifelse((pvalue[i,3] <= 0.05 && pvalue[i,4] > 0.05), 1, 0)
    }
    
    out <- list("pvalue"=pvalue, "T.obs"=T.obs,
                "F_bar_power_A"=mean(F.sigA), "E_bar_power_A"=mean(E.sigA),
                "F_bar_power_AB"=mean(F.sigAB), "E_bar_power_AB"=mean(E.sigAB))
  }
  
  Ebar.power.A <- Fbar.power.A <- Fbar.power.AB <- Ebar.power.AB <- 0
  Fbar.power.AA <- Ebar.power.AA <- data.frame()
  Fbar.power.BB <- Ebar.power.BB <- data.frame()
  Fbar.power.AB <- Ebar.power.AB <- data.frame()
  
  TT.obs <- ppvalue <- data.frame()
  
  power <- 0
  
  #This is where the magic happens
  while(power < pwr) { 
    result <- getpower()
    pvalue <- result$pvalue
    Fbar.powerA <- result$F_bar_power_A
    Ebar.powerA <- result$E_bar_power_A
    
    Fbar.powerAB <- result$F_bar_power_AB
    Ebar.powerAB <- result$E_bar_power_AB
    
    power <- result$F_bar_power_AB
    
    T.obs <- result$T.obs
    
    cat("Fbar.powerA =", Fbar.powerA, "Ebar.powerA =", Ebar.powerA, 
        "Fbar.powerAB =", Fbar.powerAB, "Ebar.powerAB =", Ebar.powerAB,
        "group size =", n,"\n")
    
    Fbar.power.AA <- rbind(Fbar.power.AA, data.frame("group_n"=n, "Fbar.powerA"=Fbar.powerA)) 
    Ebar.power.AA <- rbind(Ebar.power.AA, data.frame("group_n"=n, "Ebar.powerA"=Ebar.powerA)) 
    Fbar.power.AB <- rbind(Fbar.power.AB, data.frame("group_n"=n, "Fbar.powerB"=Fbar.powerAB)) 
    Ebar.power.AB <- rbind(Ebar.power.AB, data.frame("group_n"=n, "Ebar.powerB"=Ebar.powerAB)) 
    
    ppvalue  <- rbind(ppvalue, data.frame("group_n"=n, "pvalue"=pvalue))
    TT.obs  <- rbind(TT.obs, data.frame("group_n"=n, "T.obs"=T.obs))
    
    n <- n + 2
  }
  
  out <- list("Fbar.powerA" = Fbar.power.AA, "Ebar.powerA" = Ebar.power.AA, 
              "Fbar.powerAB"= Fbar.power.AB, "Ebar.powerAB"= Ebar.power.AB,
              "pvalue"= ppvalue, "T.obs"=T.obs)
  return(out) 
}
