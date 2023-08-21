#' initial values for Clogistic2 and Clogistic_group
#'
#' @param indexy index for y (area), should be same for same areas or units
#' @param y response
#' @param x covariates with clustered effects
#' @param sweights sampling weights
#' @param Ni population size in each unit (area); a vector of m (number of units)
#' @param lam0
#' @param tol
#' @param maxiter
#'
#' @return
#' @export
#'
#' @examples
#'
initial2 <- function(indexy,y, x, sweights, Ni,lam0  = 0.0001, tol = 1e-5, maxiter = 500)
{
  ns <- as.numeric(table(indexy))
  wt <- rep(ns/Ni, ns) * sweights

  ### first fit a logistic regression with all are
  #in one group to get the starting values
  x <- as.matrix(x)

  ncx <- ncol(x)
  n0 <- length(indexy)
  uindex <- unique(indexy)
  nr <- length(uindex)

  Ip <- diag(1,ncx,ncx)
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }
  Am <- D %x% Ip
  AtA <- t(D)%*%D %x% Ip

  Xm <- matrix(0, n0, nr*ncx)
  for(i in 1:nr)
  {
    Xm[indexy == uindex[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uindex[i],]
  }

  res0 <- glm(y ~ 0+x, family ="binomial")
  beta00 <- rep(coef(res0),nr)

  lin <- Xm %*% beta00
  mu <- 1/(1+exp(-lin))
  V <- mu * (1-mu)

  beta_cur <- beta00
  for(r in 1:maxiter)
  {
    beta_r <- beta_cur +
      solve(t(as.numeric(wt*V)*Xm) %*% Xm + lam0 *AtA) %*% (t(Xm) %*% (wt*(y-mu)) - lam0*AtA %*%beta_cur)

    diffnorm <-  norm(beta_r - beta_cur)
    beta_cur <- beta_r
    lin <-  Xm %*% beta_cur
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)

    if(diffnorm <= tol)
    {
      break
    }

  }

  est <- list( beta = beta_cur, betam = matrix(beta_cur, nr, byrow = TRUE))
  return(est)
}
