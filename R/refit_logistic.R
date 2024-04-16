#' refit based on given cluster information and sampling weights three models
#' @param obj model object
#' @param indexy area index
#' @param y response
#' @param cluster cluster structure
#' @param wts sampling weights
#' @param N population size in each unit (area); a vector of m (number of units)
#' @return
#' @export
#'
#' @import survey

refit_Lm1<- function(indexy, y, z, x, cluster, wts, N)
{
  x <- as.matrix(x)
  ncx <- ncol(x)
  nr <- length(cluster)
  Xm <- matrix(0, nrow(x), nr*ncx)
  uniqxy <- unique(indexy)
  ng <- length(unique(cluster))

  ns <- as.numeric(table(indexy))
  wtilde <- rep(1/N, ns) * wts

  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }

  W <- matrix(0, nr, ng)
  W[cbind(1:nr,cluster)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- cbind(Xm%*%W, z)

  datn <- data.frame(y = y, x = Ux, wts = wtilde)
  colnames(datn) <- c("y",paste("x", 1:ncol(Ux),sep=""),"wts")

  design <- svydesign(ids = ~1, weights = ~wts, data = datn )
  for1 <- as.formula(paste("y","~0+",paste(paste("x", 1:ncol(Ux),sep=""),collapse = "+")))
  resw <- svyglm(for1,design, family = "binomial")

  est <- coef(resw)

  eta_est <- est[-(1:(ng*ncx))]
  beta_est <- matrix(est[1:(ng*ncx)],ng, byrow = TRUE)
  beta_est <- beta_est[cluster,]
  muhat <- 1/ (1+exp(-as.matrix(Ux) %*% as.numeric(est)))

  out <- list(eta = eta_est, beta = beta_est, muhat = muhat)
  return(out)
}


#### model 2 ####
## all x have groups, no common group
refit_Lm2 <- function(indexy, y, x, cluster, wts, N)
{
  x <- as.matrix(x)
  ncx <- ncol(x)
  nr <- length(cluster)
  Xm <- matrix(0, nrow(x), nr*ncx)
  uniqxy <- unique(indexy)
  ng <- length(unique(cluster))

  ns <- as.numeric(table(indexy))
  wtilde <- rep(1/N, ns) * wts

  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }

  W <- matrix(0, nr, ng)
  W[cbind(1:nr,cluster)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- Xm%*%W

  datn <- data.frame(y = y, x = Ux, wts = wtilde)
  colnames(datn) <- c("y",paste("x", 1:ncol(Ux),sep=""),"wts")

  design <- svydesign(ids = ~1, weights = ~wts, data = datn )
  for1 <- as.formula(paste("y","~0+",paste(paste("x", 1:ncol(Ux),sep=""),collapse = "+")))
  resw <- svyglm(for1,design, family = "binomial")
  est <- coef(resw)
  glm(for1, family = "binomial", data = datn)

  beta_est <- matrix(est[1:(ng*ncx)],ng, byrow = TRUE)
  beta_est <- beta_est[cluster,]
  muhat <- 1/ (1+exp(-as.matrix(Ux) %*% as.numeric(est)))

  out <- list(beta = beta_est, muhat = muhat)

  return(out)

}


### model 3 ###
## each coordinate has its own group information
# groupmat is the true group structure
# group is the group index for coordinate
refit_Lm3 <- function(indexy, y, x, group, clustermat,wts, N)
{
  x <- as.matrix(x)
  ncx <- ncol(x)
  nr <- nrow(clustermat)
  clustermat <- clustermat[,group]
  ngest <- apply(clustermat, 2, function(x){length(unique(x))})
  ngtotal <- sum(ngest)

  ns <- as.numeric(table(indexy))
  wtilde <- rep(1/N, ns) * wts

  Xm <- matrix(0, nrow(x), nr*ncx)
  uniqxy <- unique(indexy)

  for(i in 1:nr)
  {
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]
  }

  ngestcum <- c(0,cumsum(ngest)[-ncx])
  wfun <- function(i)
  {
    W <- matrix(0, ncx, ngtotal)
    W[cbind(1:ncx,clustermat[i,] + ngestcum)] <- 1
    return(W)
  }
  Wmat <- do.call("rbind",lapply(1:nr, wfun))

  Ux <- Xm %*% Wmat

  datn <- data.frame(y = y, x = Ux, wts = wtilde)
  colnames(datn) <- c("y",paste("x", 1:ncol(Ux),sep=""),"wts")

  design <- svydesign(ids = ~1, weights = ~wts, data = datn )
  for1 <- as.formula(paste("y","~0+",paste(paste("x", 1:ncol(Ux),sep=""),collapse = "+")))
  resw <- svyglm(for1,design, family = "binomial")
  est <- coef(resw)

  muhat <- 1/ (1+exp(-as.matrix(Ux) %*% as.numeric(est)))

  estm <- Wmat %*% est
  beta <- matrix(estm, nr, byrow=TRUE)

  out <- list(beta=beta, muhat = muhat)

  return(out)
}



