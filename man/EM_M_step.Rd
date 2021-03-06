% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_M_step.R
\name{EM_M_step}
\alias{EM_M_step}
\title{EM M-step Function}
\usage{
EM_M_step(
  Ymatrix,
  Xmatrix,
  split,
  t,
  N,
  pden.post,
  VarMatrix,
  Beta,
  K,
  Lambda = 0
)
}
\arguments{
\item{Ymatrix}{matrix of dependent variables, multiple dimensions eligible}

\item{Xmatrix}{matrix of independent variables}

\item{split}{index of columns with no penalty}

\item{t}{length of trajectory}

\item{N}{number of trajectories}

\item{pden.post}{probablity density from E-step}

\item{VarMatrix}{covariance matrix}

\item{Beta}{coefficient matrix}

\item{K}{number of clusters}

\item{Lambda}{penalty coefficient}
}
\value{
list of outputs used for M-step
}
\description{
M-step to calculate the coefficients and covariance matrix in EM algorithm
}
