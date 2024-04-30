#' initial values for Clogistic1
#'
#' @param indexy index for y (area), should be same for same areas or units
#' @param y response
#' @param z covariates with common effects
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

initial1 <- function(indexy,y, z, x, sweights, Ni, lam0  = 0.0001, tol = 1e-5, maxiter = 500)
{

  ns <- as.numeric(table(indexy))
  #wt <- rep(ns/Ni, ns) * sweights

  Nbar <- rep(mean(Ni), length(Ni))
  wt <- rep(1/Nbar, ns) * sweights

  ### first fit a logistic regression with all are
  #in one group to get the starting values
  x <- as.matrix(x)
  z <- as.matrix(z)

  ncz <- ncol(z)
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

  res0 <- glm(y ~ 0+cbind(z,x), family ="binomial")


  eta00 <- coef(res0)[1:ncz]
  beta00 <- rep(coef(res0)[(ncz+1):(ncz+ncx)],nr)

  lin <- z %*% eta00 + Xm %*% beta00
  mu <- 1/(1+exp(-lin))
  V <- mu * (1-mu)

  eta_cur <- eta00
  beta_cur <- beta00
  for(r in 1:maxiter)
  {
    eta_r <- eta_cur + solve(t(as.numeric(wt*V)*z) %*% z) %*% t(z) %*% (wt*(y-mu))
    beta_r <- beta_cur +
      solve(t(as.numeric(wt*V)*Xm) %*% Xm + lam0 *AtA) %*% (t(Xm) %*% (wt*(y-mu)) - lam0*AtA %*%beta_cur)


    diffnorm <- norm(eta_cur - eta_r) + norm(beta_r - beta_cur)
    eta_cur <- eta_r
    beta_cur <- beta_r
    lin <- z %*% eta_cur + Xm %*% beta_cur
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)

    if(diffnorm <= tol)
    {
      break
    }

  }

  est <- list(eta = eta_cur, beta = beta_cur, betam = matrix(beta_cur, nr, byrow = TRUE))
  return(est)
}



