 #' Get estimates  based on given data, model and lambda search
#' @param area area vector
#' @param y response, sampled data
#' @param x covariates without intercept, sampled data
#' @param xunit population level data, first column the area index
#' @param wts sampling weights
#' @param N population size in each unit (area); a vector of m (number of units)
#' @param model three models, "intercept","creg" and "ccreg"
#' @param group group information for ccreg model
#' @param lambda a vector of lambda values
#' @return
#' @export
#'
#' @examples



Clogistic_SAE <- function(area, y, x, xunit,
                          wts, N, model, group=NULL, lambda,
                          standardize = TRUE, center = TRUE)
{

  nt <- length(y)
  m <- length(unique(area))
  nivec <- as.numeric(table(area))

  x <- as.matrix(x)
  xunit <- as.matrix(xunit)

  xintercept <- matrix(1, nrow = nt)

  if(standardize)
  {
    #y_sc <- scale(y, center = center, scale = sd(y))
    x_sc <- scale(x, center = center, scale = apply(x, 2, sd, na.rm = TRUE))
  }else{
    #y_sc <- y
    x_sc <- x
  }

  ncx <- ncol(x)

  nlam <- length(lambda)

  bic <- rep(0, nlam)
  ngest <- rep(0, nlam)

  if(model == "intercept")
  {
    Cm <- log(m)
    init0 <- initial1(indexy = area, y = y,z = x_sc, x = xintercept,
                      sweights = wts,
                      Ni = N)
    betam0 <- init0$betam
    eta00 <- init0$eta

    for(j in 1:nlam)
    {
      resj <- Clogistic1(indexy = area, y = y, z = x_sc, x = xintercept,
                          sweights = wts, Ni = N,
                          betam0 = betam0, eta00 = eta00,
                          lam = lambda[j])
      ng <- length(unique(resj$cluster))# number of clusters

      ngest[j] <- ng
      bic[j] <- 2*resj$nll/m + Cm *ng*log(nt)/nt
    }

    BICm <- data.frame(ng = ngest, bic = bic)

    resm <- Clogistic1(indexy = area, y = y, z = x_sc, x = xintercept,
                       sweights = wts, Ni = N,
                       betam0 = betam0, eta00 = eta00,
                       lam = lambda[which.min(bic)])
    # based on cluster get the estimate


    refit_est <- estSAE_logistic(obj = resm, indexy = area, y = y,
                          x = x, xunit = xunit, wts = wts,
                          N = N, model = "intercept", group = group)

  }

  if(model == "creg")
  {
    Cm <- log(m*(1+ncol(x)))
    init0 <- initial2(indexy = area, y = y, x = cbind(1,x_sc),
                      sweights = wts, Ni = N)
    betam0 <- init0$betam

    for(j in 1:nlam)
    {
      resj <- Clogistic2(indexy = area, y = y, x = cbind(1,x_sc),
                         sweights = wts, Ni = N,
                         betam0 = betam0,
                         lam = lambda[j])
      ng <- length(unique(resj$cluster))# number of clusters

      ngest[j] <- ng
      bic[j] <- 2*resj$nll/m + Cm *ng*(1+ncx)*log(nt)/nt
    }

    BICm <- data.frame(ng = ngest, bic = bic)

    resm <- Clogistic2(indexy = area, y = y, x = cbind(1,x_sc),
                       sweights = wts, Ni = N,
                       betam0 = betam0,
                       lam = lambda[which.min(BICm$bic)])

    refit_est <- estSAE_logistic(obj = resm, indexy = area, y = y,
                                 x = x, xunit = xunit, wts = wts,
                                 N = N, model = "creg", group = group)

  }


  if(model == "ccreg")
  {
    Cm <- log(m*(1+ncol(x)))
    init0 <- initial2(indexy = area, y = y, x = cbind(1,x_sc),
                      sweights = wts, Ni = N)
    betam0 <- init0$betam

    for(j in 1:nlam)
    {
      resj <- Clogistic_group(indexy = area, y = y, x = cbind(1,x_sc),
                         sweights = wts, Ni = N,group = group,
                         betam0 = betam0,
                         lam = lambda[j])
      ng <- sum(apply(resj$cluster,2, function(x) length(unique(x))))

      ngest[j] <- ng
      bic[j] <- 2*resj$nll/m + Cm *ng*log(nt)/nt
    }

    BICm <- data.frame(ng = ngest, bic = bic)

    resm <- Clogistic_group(indexy = area, y = y, x = cbind(1,x_sc),
                       sweights = wts, Ni = N, group = group,
                       betam0 = betam0,
                       lam = lambda[which.min(BICm$bic)])

    refit_est <- estSAE_logistic(obj = resm, indexy = area, y = y,
                                 x = x, xunit = xunit, wts = wts,
                                 N = N, model = "ccreg", group = group)

  }


  out <- list(resm = resm, BICm = BICm,
              refit = refit_est,
              cluster = resm$cluster,
              estY_SR = refit_est$estY_SR,
              estY_syn = refit_est$estY_syn,
              estY_comp = refit_est$estY_comp,
              init = betam0)


  return(out)
}


