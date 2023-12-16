#' estimation for area proportions based on clustered intercept model
#' @param obj model object
#' @param indexy area index
#' @param y response
#' @param z covariates common regression coef
#' @param sweights sampling weights
#' @param Ni population size in each unit (area); a vector of m (number of units)
#' @param zpopd In each area, the number of units for each combination,  
#' first column area, second column number of units in each domain Nc
#' same as original z in a dummies form with same column names
#' 
#' 
estSAE_intercept <- function(obj,indexy, y, z, sweights, 
                              zpopd, Ni)
{
  betaest <- obj$beta
  etaest <- obj$eta
  clusterest <- obj$cluster
  
  uindex <- unique(indexy)
  m <- length(uindex)
  
  betaest_area <- data.frame(area = unique(indexy), 
                             intercept = betaest) 
  
  # phatc in each 
  phatc <- zpopd %>% left_join(betaest_area, by = "area")
  phatc <- phatc %>%
    mutate(phat = 1/(1+exp(-intercept - as.numeric(as.matrix(zpopd[,-(1:2)]) %*% etaest))))
  
  
  Ybarest_SR <- Ybarest_syn <- Ybarest_comp <-  rep(0,m)

  for(i in 1:m)
  {
    indexi <- indexy == uindex[i]
    indexi_pop <- zpop[,1] ==uindex[i]
    zi_sample <- z[indexi,]
    phati_sample <- 1/(1+exp(-betaest[i] - as.matrix(zi_sample) %*% etaest))
    wtsi <- sweights[indexi]
    
    Ybarest_SR[i] <- sum(wtsi * (y[indexi] - phati_sample))/Ni[i]
  }
  
  phatcs <- phatc %>% group_by(area) %>%
    summarise(phatc = sum(Nc*phat))
  
  Ybarest_SR <- phatcs$phatc/Ni + Ybarest_SR
  Ybarest_syn <- phatcs$phatc/Ni
  
  datay <- data.frame(area = indexy, y = y, z)
  sample_ag <- datay %>% group_by_at(vars(area,colnames(z))) %>%
    summarise(nc = n(), sumy = sum(y)) %>% 
    right_join(phatc, by = c("area",colnames(z))) %>%
    replace(is.na(.), 0)%>%
    ungroup() %>% group_by(area) %>%
    summarise(est = sum((Nc - nc)*phat),sumy = sum(sumy))
  
  Ybarest_comp <- sample_ag$sumy/Ni + sample_ag$est/Ni
  
  
  return(list(SR = Ybarest_SR, SYN = Ybarest_syn, Comp = Ybarest_comp))
  
}


