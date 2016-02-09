
################################################################################
# A latent factor model using a spline basis: This minimizes
# ||Y - H B W^{T}||_{2}^{2}, where H is a known basis matrix, over B and W.
# It only optimizes over nonmissing entries in Y.
################################################################################

#' @title Given Y, H, and W, optimize B stochastic gradient descent
#' @description This is minimizing ||Y - HBW^{T}||_{2}^{2} over B, with all
#' other matrices assumed known.
#' @param Y An N x J real matrix, with censored values set to NA. The standard
#' application we have in mind is the sample x OTU count matrix.
#' @param H The spline basis matrix, from which the latent sources arise (as
#' the linear mixture HB.)
#' @param W The samples-to-latent-source coefficients matrix, assumed known.
#' @param nu The learning rate for the stochastic gradient descent. Arbitrarily
#' defaults to 1e-3.
#' @param n_iter The number of sweeps over all entries of B via stochastic
#' gradient descent.
#' @return A list with the following elements, \cr
#'   $obj The RSS after each update to an entry of B. \cr
#'   $B The optimized value of B.
#' @export
#' @examples
#' # generate data
#' N <- 150
#' P <- 20
#' K <- 5
#' L <- 6
#'
#' H <- bs(1:N, df = L, degree = 1)
#' W <- matrix(rnorm(P * K), P, K)
#' B <- matrix(rnorm(L * K), L, K)
#' E <- matrix(.5 * rnorm(N * P), N, P)
#'
#' Y <- H %*% B %*% t(W) + E
#' Y[sample(n * p, n * p * .4)] <- NA # 40% missing at random
#'
#' # fit the model
#' B0 <- matrix(rnorm(L * K), L, K)
#' B_res <- get_B(Y, H, W, B0)
get_B <- function(Y, H, W, B0, nu = .0001, n_iter = 100) {
  # get problem dimensions
  N <- nrow(Y)
  P <- ncol(Y)
  K <- ncol(W)
  L <- ncol(H)

  # intialize results
  obj <- vector(length = n_iter * L * K)
  ix <- 1
  B <- B0

  # SGD sweep over entries of B, seeing each entry n_iter times
  for(iter in seq_len(n_iter)) {
    if(iter %% 20 == 0) cat(sprintf("iter %d\n", iter))
    for(l in seq_len(L)) {
      for(k in seq_len(K)) {
        A <- 2 * (Y - H %*% B %*% t(W))
        A[is.na(Y)] <- 0
        B[l, k] <- B[l, k] + nu * t(H[, l]) %*% A %*% W[, k]
        obj[ix] <- sum((Y - H %*% B %*% t(W))^2, na.rm = T)
        ix <- ix + 1
      }
    }
  }
  list(B = B, obj = obj)
}
