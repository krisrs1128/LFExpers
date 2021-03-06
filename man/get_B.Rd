% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_B}
\alias{get_B}
\title{Given Y, H, and W, optimize B coordinate descent}
\usage{
get_B(Y, H, W, B0, nu = 1e-05, n_iter = 100L)
}
\arguments{
\item{Y}{An N x J real matrix, with censored values set to NA. The standard
application we have in mind is the sample x OTU count matrix.}

\item{H}{The spline basis matrix, from which the latent sources arise (as
the linear mixture HB.)}

\item{W}{The samples-to-latent-source coefficients matrix, assumed known.}

\item{nu}{The learning rate for the coordinate descent. Arbitrarily
defaults to 1e-3.}

\item{n_iter}{The number of sweeps over all entries of B via coordinate
descent.}
}
\value{
A list with the following elements, \cr
  $obj The RSS after each update to an entry of B. \cr
  $B The optimized value of B.
}
\description{
This is minimizing ||Y - HBW^{T}||_{2}^{2} over B, with all
other matrices assumed known.
}
\examples{
# generate data
N <- 150
P <- 20
K <- 5
L <- 6

library("splines")
H <- bs(1:N, df = L, degree = 1)
W <- matrix(rnorm(P * K), P, K)
B <- matrix(rnorm(L * K), L, K)
E <- matrix(.5 * rnorm(N * P), N, P)

Y <- H \%*\% B \%*\% t(W) + E
Y[sample(N * P, N * P * .4)] <- NA # 40\% missing at random

# fit the model
B0 <- matrix(rnorm(L * K), L, K)
B_res <- get_B(Y, H, W, B0)
}

