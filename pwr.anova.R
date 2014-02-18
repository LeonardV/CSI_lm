pwr.anova <- function(f=f, n=group_size, ui, seed=54321, means=NULL, order=c(1:k), sigma=diag(k)) {
  
  power <- list()
  FFbar.powerAB <- EEbar.powerAB <- data.frame()
  
  for(e in 1:length(f)) {
    
    es <- f[e]
    
    if(is.null(means)) { 
      i=1:k
      #compute equal difference scores, d, between the means
      d=(2*sqrt(k)*es) / sqrt(sum((2*i-1-k)^2))#*
      #compute first (read lowest) mean mu
      mu.i=( -(k-1)*d ) / 2
      #compute k means
      means=c(mu.i, mu.i + 1:(k-1)*d)
      
      #ordering means
      means <- means[order]
    }  
    
    out <- run_sim(k=k, n=rep(group_size,k), means=means,  
                   model='y ~ grp', ui=ui, 
                   meq=0, sims=sims, pwr=group_size, method="weights",
                   verbose=TRUE, seed=seed) 
    
    
    power[[e]] <- out$Fbar.powerAB      
    
  } 
  
  return(power)
}
