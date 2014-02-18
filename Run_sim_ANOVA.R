run_sim <- function(k=k, n=n, means=means, sigma=diag(k), model=model, ui=ui, 
                    meq=0, sims=sims, pwr=0.80, R=10000, 
                    parallel="multicore", ncpus, cl=NULL, 
                    method=c("weights", "bootstrap", "robust"),
                    double.bootstrap.R=double.bootstrap.R,
                    verbose=FALSE, seed=12345) {
     
  getpower <- function() {
    
    E.sigAB <- F.sigAB <- F.sigB <- F.sigA <- E.sigB <- E.sigA <- rep(NA, sims)
    T.obs <- pvalue <- matrix(as.numeric(NA), sims, 4)
        
    for(i in 1:sims) {
      #create ANOVA data, sigma is an identity-matrix
      sim.data <- cbind(c(matrix(mvtnorm:::rmvnorm(n[1], mean=means, sigma=sigma), nrow=n[1])))
      
      #fit model
      fit.lm <- lm(sim.data ~ factor(group))
      data <- fit.lm$model
        names(data) <- paste(c("y", "grp"))
      
      #methods for calculation p-values
      
      #robust estimation (currently not used)
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
      #determine weights  
      } else if(method == "weights") {
        #Type A, see Silvapulle and Sen, 2005, p99-100
        r=k-1
        q=nrow(ui)
        p=k
        idx=q:0
        N2=nrow(fit.lm$model)
        nu=fit.lm$df.residual
        df1a=r-q+idx
        df2a=N2-p+q-idx
        
        #Type B
        df1b=idx
        df2b=N2-k
        
        #model implied covariance matrix
        cov <- vcov(fit.lm)
    
        #in case of only inequallity constraints
        #the weights are computed using the R function ic.weights() from the 
        #ic.infer packages
        if(meq==0) {
          wt.bar <- ic.infer:::ic.weights(ui %*% cov %*% t(ui)) 
        } 
        #in case of inequality and equality constraints
        else if(meq > 0) {
          wt.bar <- ic.weights(solve(solve(uiw %*% cov %*% t(uiw))[-(1:meq), 
                                                                   -(1:meq)])) 
        }
        
        #compute observed F-bar and E-bar-square
        fit <- Fbar(model='y ~ .', data, R=0, 
                    parallel="multicore", ncpus=1,
                    pvalue=FALSE, verbose=TRUE, meq=0, 
                    ui=ui, p.distr="N", df=4)
          
        Fbar.A <- fit$T.obs[1] ; Fbar.B <- fit$T.obs[2]
        Ebar.A <- fit$T.obs[3] ; Ebar.B <- fit$T.obs[4]
        
        ###################################################################################################################
        ##Compute p-values for hypothesis test Type A, see Silvapulle and Sen, 2005, p99-100
        #These are adjusted functions of the pbetabar() function from the ic.infer package.
        pbar.A <- function(x1, x2, df1a, df2a, nu, wt) {
          if (x1 <= 0) { ##also x2??
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
          Fbar.pA <- 1-out.A$cdf.fbarA ; Ebar.pA <- 1-out.A$cdf.ebarA
        
        
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
          Fbar.pB <- 1-out.B$cdf.fbarB ; Ebar.pB <- 1-out.B$cdf.ebarB
        
        fit$p.value <- attr(fit, "p.value") <- c(Fbar.pA, Fbar.pB, Ebar.pA, Ebar.pB)
        fit$T.obs <- attr(fit, "T.obs")     <- c(Fbar.A, Fbar.B, Ebar.A, Ebar.B)
      }
      
      #store p-values and observed test statistics
      pvalue[i,1:4] <- fit$p.value[1:4] ; T.obs[i,1:4] <- fit$T.obs[1:4]
      
      #assign 1 to significant result and 0 otherwise
      F.sigA[i] <- ifelse(pvalue[i,1] <= .05, 1, 0)  
      F.sigB[i] <- ifelse(pvalue[i,2] <= .05, 1, 0)
      E.sigA[i] <- ifelse(pvalue[i,3] <= .05, 1, 0)
      E.sigB[i] <- ifelse(pvalue[i,4] <= .05, 1, 0)
      
      #combine hypothesis Test A and Type B into Type A*
      F.sigAB[i] <- ifelse((pvalue[i,1] <= 0.05 && pvalue[i,2] > 0.05), 1, 0) 
      E.sigAB[i] <- ifelse((pvalue[i,3] <= 0.05 && pvalue[i,4] > 0.05), 1, 0)
      
      
    }
    
    out <- list("pvalue"=pvalue, "T.obs"=T.obs,
                "F_bar_power_A"=mean(F.sigA), "F_bar_power_B"=mean(F.sigB),
                "E_bar_power_A"=mean(E.sigA), "E_bar_power_B"=mean(E.sigB),
                "F_bar_power_AB"=mean(F.sigAB), "E_bar_power_AB"=mean(E.sigAB))
  }
  

  Ebar.power.A <- Fbar.power.A <- Ebar.power.B <- Fbar.power.B <- 0
  Fbar.power.AB <- Ebar.power.AB <- 0
  Fbar.power.AA <- Fbar.power.BB <- Ebar.power.AA <- Ebar.power.BB <- data.frame()
  Fbar.power.AB <- Ebar.power.AB <- data.frame()
  ppvalue <- TT.obs <- data.frame()
  
  power <- 0
  
  #This is where the magic happens
  while(power < pwr) { 
    group <- rep(1:k, each=n[1])
    result <- getpower()
    pvalue <- result$pvalue
    
    Fbar.powerA <- result$F_bar_power_A
    Fbar.powerB <- result$F_bar_power_B
    Ebar.powerA <- result$E_bar_power_A
    Ebar.powerB <- result$E_bar_power_B
    
    #since the power for the F-bar and E-bar-square are approximately the same
    #we focus only on the F-bar.
    power <- Fbar.powerAB <- result$F_bar_power_AB
    Ebar.powerAB <- result$E_bar_power_AB
    T.obs <- result$T.obs
    
    cat("Fbar.powerA =", Fbar.powerA, "Fbar.powerB =", Fbar.powerB, 
        "Ebar.powerA =", Ebar.powerA, "Ebar.powerB =", Ebar.powerB,
        "Fbar.powerAB =", Fbar.powerAB, "Ebar.powerAB =", Ebar.powerAB,
        "group size =", n[1],"\n")
    
    Fbar.power.AA <- rbind(Fbar.power.AA, data.frame("group_n"=n[1], "Fbar.powerA"=Fbar.powerA)) 
    Fbar.power.BB <- rbind(Fbar.power.BB, data.frame("group_n"=n[1], "Fbar.powerB"=Fbar.powerB)) 
    Ebar.power.AA <- rbind(Ebar.power.AA, data.frame("group_n"=n[1], "Ebar.powerA"=Ebar.powerA)) 
    Ebar.power.BB <- rbind(Ebar.power.BB, data.frame("group_n"=n[1], "Ebar.powerB"=Ebar.powerB)) 
    
    Fbar.power.AB <- rbind(Fbar.power.AB, data.frame("group_n"=n[1], "Fbar.powerAB"=Fbar.powerAB)) 
    Ebar.power.AB <- rbind(Ebar.power.AB, data.frame("group_n"=n[1], "Ebar.powerAB"=Ebar.powerAB)) 
    
    ppvalue  <- rbind(ppvalue, data.frame("group_n"=n[1], "pvalue"=pvalue))
    TT.obs  <- rbind(TT.obs, data.frame("group_n"=n[1], "T.obs"=T.obs))
    
    n <- n + 2
  }
  
  out <- list("Fbar.powerA" = Fbar.power.AA, "Fbar.powerB"= Fbar.power.BB, 
              "Ebar.powerA" = Ebar.power.AA, "Ebar.powerB"= Ebar.power.BB,
              "Fbar.powerAB"= Fbar.power.AB, "Ebar.powerAB"= Ebar.power.AB,
              "pvalue"= ppvalue, "T.obs"=T.obs)
  
  return(out) 
}
