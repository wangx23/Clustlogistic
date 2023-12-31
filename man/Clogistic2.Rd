% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clogistic2.R
\name{Clogistic2}
\alias{Clogistic2}
\title{Clustered coefficient logistic regression model}
\usage{
Clogistic2(
  indexy,
  y,
  x,
  sweights,
  Ni,
  betam0,
  nu = 1,
  gam = 3,
  lam = 0.5,
  inner_loop = 500,
  maxiter = 1000,
  tol = 0.001,
  tol1 = 0.001,
  tol2 = 0.001
)
}
\arguments{
\item{indexy}{index for y (area), should be same for same areas or units}

\item{y}{response}

\item{x}{covariates with clustered effects}

\item{sweights}{sampling weights}

\item{Ni}{population size in each unit (area); a vector of m (number of units)}

\item{betam0}{initial values for coefficient of x}

\item{nu}{a fixed value in ADMM}

\item{gam}{a fixed value in SCAD}

\item{lam}{tuning parameter}

\item{tol2}{}
}
\description{
Clustered coefficient logistic regression model
}
