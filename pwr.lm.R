pwr.lm <- function(p=p, n=n, ui=ui, cov.x=cov.x, f2=f2, seed=3013073) {
  
  power <- list()
  FFbar.powerAB <- EEbar.powerAB <- data.frame()
  
  for(e in 1:length(f2)) {
    
    es <- f2[e]
    
    out <- run_sim(p=p, cov.x=cov.x, n=n, R2=(es/(1+es)), 
                   model='y ~ grp', ui=ui, 
                   meq=0, sims=sims, pwr=group_size, method="weights",
                   verbose=TRUE, seed=seed) 
    
    
    power[[e]] <- out$Fbar.powerAB      
    
  } 
  
  return(power)
}


