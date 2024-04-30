#' Clustered coefficient logistic regression model with group covariates
#'
#' @param indexy index for y (area), should be same for same areas or units
#' @param y response
#' @param z covariates with common effects
#' @param x covariates with clustered effects
#' @param sweights sampling weights
#' @param Ni population size in each unit (area); a vector of m (number of units)
#' @param group group index for covariates
#' @param betam0 initial values for coefficient of x
#' @param eta00 initial values for z
#' @param nu a fixed value in ADMM
#' @param gam a fixed value in SCAD
#' @param lam tuning parameter
#' @param inner_loop
#' @param maxiter
#' @param tol
#' @param tol1
#' @param tol2
#'
#' @return
#' @export
#'
#' @examples

Clogistic_group <- function(indexy, y, x, sweights, Ni, group, betam0,
                       nu = 1, gam = 3, lam = 0.5, inner_loop= 500, maxiter=1000,
                       tol = 0.001, tol1 = 0.001, tol2 = 0.001)
{

  ### data info#####
  ns <- as.numeric(table(indexy))
  #wt <- rep(ns/Ni, ns) * sweights
  Nbar <- rep(mean(Ni), length(Ni))
  wt <- rep(1/Nbar, ns) * sweights

  x <- as.matrix(x)
  ncx <- ncol(x)
  n0 <- length(indexy)
  uindex <- unique(indexy)
  nr <- length(uindex)
 #nbar <- mean(ns)
  nbar <- 1

  ugroup <- unique(group)
  ngroup <- length(ugroup)

  ### prepare
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

  ##### initial values ####

  deltam <- D %*% betam0
  vm <-  matrix(0, nr*(nr-1)/2 , ncx)

  beta_cur <- c(t(betam0))

  lin <- Xm %*% beta_cur
  mu <- 1/(1+exp(-lin))
  V <- mu * (1-mu)

  for(m in 1:maxiter)
  {

    ### inner eta and beta
    for(r in 1:inner_loop)
    {
      temp1 <-  1/nbar * t(Xm) %*% (wt*(y-mu))
      temp2 <- nu * AtA %*% beta_cur - nu* c(t(deltam -  vm/nu) %*% D)
      beta_r <- beta_cur +
        solve(1/nbar*t(as.numeric(wt*V)*Xm) %*% Xm + nu *AtA) %*% (temp1 - temp2)

      diffinner <- norm(beta_r - beta_cur)
      beta_cur <- beta_r
      lin <-  Xm %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)

      if(diffinner < tol)
      {break}
    }

    betam <- matrix(beta_r, nr, byrow = TRUE)

    ### update deltam
    psim <- D %*% betam + vm/nu

    deltam_new <- deltam

    ### update group ##

    for(j in 1:ngroup)
    {
      indexj <- group == ugroup[j]
      deltam_new[,indexj] <- t(sapply(1:nrow(psim),function(xx) scad(psim[xx,indexj],lam,nu,gam)))
    }

    vm <- vm + nu * (D %*% betam - deltam_new)
    sm <- norm(deltam_new - deltam)
    rm <- norm(D %*% betam  - deltam_new)

    deltam <- deltam_new

    if(sm <= tol1 & rm <= tol2)
    {break}

  }

  groupest <- matrix(0, nr, ngroup)
  for(j in 1:ngroup)
  {
    indexj <- group == ugroup[j]
    groupest[,j] <- getgroup(t(deltam_new[,indexj]),nr)
  }


  lin <- Xm %*% beta_r
  muest  <- 1/(1+exp(-lin))

  nllvalue <-  -sum(wt* (y * lin - log(1+exp(lin))))

  output <- list(beta = beta_cur, betam = betam,
                 deltam = deltam_new, cluster = groupest,
                 lin = lin, mu  = muest, nll = nllvalue,
                 sm = sm, rm = rm, niter = m )
  return(output)
}
