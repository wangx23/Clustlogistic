#' estimation for area proportions based on models three estimator
#' @param obj model object
#' @param indexy area index
#' @param y response
#' @param x covariates
#' @param xunit population level x
#' @param sweights sampling weights
#' @param N population size in each unit (area); a vector of m (number of units)
#'
#'
#'
#' @return
#' @export
estSAE_logistic <- function(obj,indexy, y, x, xunit, wts,N,
                               model, group = NULL)
{

  cluster_est <- obj$cluster
  nivec <- as.numeric(table(indexy))
  nt <- sum(nivec)
  uindexy <- unique(indexy)
  m <- length(nivec)

  x <- as.matrix(x)
  xunit <- as.matrix(xunit)

  if(model == "intercept")
  {
    xintercept <- matrix(1, nt)
    res_refit <- refit_Lm1(indexy = indexy,y = y, z = x, x = xintercept,
                          cluster = cluster_est,wts = wts,N = N)
    eta_est <- res_refit$eta
    beta_est <- res_refit$beta
    ## muhat_unit

    muhat_unit <- rep(0, nrow(xunit))
    for(ii in 1:m)
    {
      xuniti <- xunit[xunit[,1]==uindexy[ii],-1]
      muhat_unit[xunit[,1]==uindexy[ii]] <- 1/(1+exp(-as.matrix(xuniti) %*% eta_est - beta_est[ii]))
    }

    muhat_unit_mean <- aggregate(muhat_unit,by = list(area = xunit[,1]),FUN = mean)
  }

  if(model == "creg")
  {
    res_refit <- refit_Lm2(indexy = indexy,y = y, x = cbind(1,x),
                           cluster = cluster_est,wts = wts,N = N)

    eta_est <- res_refit$eta
    beta_est <- res_refit$beta
    ## muhat_unit

    muhat_unit <- rep(0, nrow(xunit))
    for(ii in 1:m)
    {
      xuniti <- xunit[xunit[,1]==uindexy[ii],-1]
      muhat_unit[xunit[,1]==uindexy[ii]] <- 1/(1+exp(-cbind(1,as.matrix(xuniti)) %*% beta_est[ii,]))
    }

    muhat_unit_mean <- aggregate(muhat_unit,by = list(area = xunit[,1]),FUN = mean)
  }


  if(model == "ccreg")
  {
    res_refit <- refit_Lm3(indexy = indexy,y = y, x = cbind(1,x),
                           clustermat = cluster_est, wts = wts,N = N,group = group)

    eta_est <- res_refit$eta
    beta_est <- res_refit$beta
    ## muhat_unit

    muhat_unit <- rep(0, nrow(xunit))
    for(ii in 1:m)
    {
      xuniti <- xunit[xunit[,1]==uindexy[ii],-1]
      muhat_unit[xunit[,1]==uindexy[ii]] <- 1/(1+exp(-cbind(1,as.matrix(xuniti)) %*% beta_est[ii,]))
    }

    muhat_unit_mean <- aggregate(muhat_unit,by = list(area = xunit[,1]),FUN = mean)
  }

  muhat <- res_refit$muhat
  subfun <- function(ii)
  {
    indexii <- indexy == uindexy[ii]
    yii <- y[indexii]
    muhatii <- muhat[indexii]

    sum((yii - muhatii)*wts[indexii])/N[ii]
  }
  estp1 <- unlist(lapply(1:length(uindexy),subfun))

  # three estimators
  estY_SR <- estp1 + muhat_unit_mean$x
  estY_syn <- muhat_unit_mean$x

  estp2 <- aggregate(y-muhat, by = list(area = indexy), sum)[,2]
  estY_comp <- muhat_unit_mean$x + estp2/N


  out <- list(estY_SR = estY_SR, estY_syn = estY_syn, estY_comp = estY_comp,
              refit = res_refit, model = model)

  return(out)

}


