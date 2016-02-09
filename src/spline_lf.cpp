#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// helpers ---------------------------------------------------------------------

// @title Replace NAs with a number
// @description Equivalent to X[is.na(X)] <- impute in R. Just needed so
// everything in get_B() will be in C++.
arma::mat replace_nas(arma::mat X, double impute) {
  for(int i = 0; i < X.n_rows; i++ ) {
    for(int j = 0; j < X.n_cols; j++) {
      if(!arma::is_finite(X(i, j))) {
	X(i, j) = impute;
      }
    }
  }
  return X;
}

// optimizing-terms ------------------------------------------------------------
//' @title Given Y, H, and W, optimize B coordinate descent
//' @description This is minimizing ||Y - HBW^{T}||_{2}^{2} over B, with all
//' other matrices assumed known.
//' @param Y An N x J real matrix, with censored values set to NA. The standard
//' application we have in mind is the sample x OTU count matrix.
//' @param H The spline basis matrix, from which the latent sources arise (as
//' the linear mixture HB.)
//' @param W The samples-to-latent-source coefficients matrix, assumed known.
//' @param nu The learning rate for the coordinate descent. Arbitrarily
//' defaults to 1e-3.
//' @param n_iter The number of sweeps over all entries of B via coordinate
//' descent.
//' @return A list with the following elements, \cr
//'   $obj The RSS after each update to an entry of B. \cr
//'   $B The optimized value of B.
//' @export
//' @examples
//' # generate data
//' N <- 150
//' P <- 20
//' K <- 5
//' L <- 6
//'
//' library("splines")
//' H <- bs(1:N, df = L, degree = 1)
//' W <- matrix(rnorm(P * K), P, K)
//' B <- matrix(rnorm(L * K), L, K)
//' E <- matrix(.5 * rnorm(N * P), N, P)
//'
//' Y <- H %*% B %*% t(W) + E
//' Y[sample(N * P, N * P * .4)] <- NA # 40% missing at random
//'
//' # fit the model
//' B0 <- matrix(rnorm(L * K), L, K)
//' B_res <- get_B(Y, H, W, B0)
//' @importFrom Rcpp evalCpp
//' @useDynLib LFExpers
// [[Rcpp::export]]
Rcpp::List get_B(arma::mat Y, arma::mat H, arma::mat W, arma::mat B0,
		     double nu = 1e-5, int n_iter = 100) {
  // get problem dimensions
  int L = H.n_cols;
  int K = W.n_cols;

  // initialize results
  arma::vec obj = arma::zeros<arma::vec>(n_iter * L * K);
  int ix = 0;
  arma::mat B = B0;

  for(int i = 0; i < n_iter; i++) {
    arma::mat A = 2 * (Y - H * B * W.t());
    for(int l = 0; l < L; l++ ) {
      for(int k = 0; k < K; k++) {
	A = replace_nas(A, 0);
	B(l, k) += nu * arma::as_scalar(H.col(l).t() * A * W.col(k));
	obj(ix) = arma::accu(arma::square(replace_nas(Y - H * B * W.t(), 0)));
	ix++;
      }
    }
  }

  return List::create(Named("B") = B, Named("obj") = obj);
}
