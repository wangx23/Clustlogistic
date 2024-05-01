#' estimation for area proportions based on model with common regression coefficients
#' @param obj model object
#' @param indexy area index
#' @param y response
#' @param x covariates common regression coef
#' @param sweights sampling weights
#' @param N population size in each unit (area); a vector of m (number of units)
#' @param xunit In each area, population level
#' @import survey
#' @return
#' @export

lgreg <- function(indexy, y, x, sweights, xunit, N)
{

  uindex <- unique(indexy)
  m <- length(uindex)
  res_trad <- refit_Lm2(indexy, y, cbind(1,x), rep(1,m), wts=sweights , N)
  betaest <- res_trad$beta[1,]

  xunit <- as.matrix(xunit)

  muhat_unit <- rep(0, nrow(xunit))
  for(ii in 1:m)
  {
    xuniti <- xunit[xunit[,1]==uindex[ii],-1]
    muhat_unit[xunit[,1]==uindex[ii]] <- 1/(1+exp(-as.matrix(xuniti) %*% betaest[-1]-betaest[1]))
  }

  muhat_unit_mean <- aggregate(muhat_unit,by = list(area = xunit[,1]),FUN = mean)

  muhat <- res_trad$muhat
  subfun <- function(ii)
  {
    indexii <- indexy == uindex[ii]
    yii <- y[indexii]
    muhatii <- muhat[indexii]
    sum((yii - muhatii)*sweights[indexii])/N[ii]
  }
  estp1 <- unlist(lapply(1:length(uindex),subfun))


  estY_SR <- estp1 + muhat_unit_mean$x
  estY_syn <- muhat_unit_mean$x

  estp2 <- aggregate(y-muhat, by = list(area = indexy), sum)[,2]
  estY_comp <- muhat_unit_mean$x + estp2/N


  return(list(estY_SR = estY_SR, estY_syn = estY_syn, estY_comp = estY_comp,
              fit = res_trad))

}
