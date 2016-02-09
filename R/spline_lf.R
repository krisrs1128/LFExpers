
################################################################################
# A latent factor model using a spline basis: This minimizes
# ||Y - H B W^{T}||_{2}^{2}, where H is a known basis matrix, over B and W.
# It only optimizes over nonmissing entries in Y.
################################################################################

# helpers ----------------------------------------------------------------------
#' @title Merge in default options for spline_lf()
#' @param opts A potentially partially specified list, which we will merge in
#' with defaults. Please include new parameters in alphabetical order. The
#' possible options, and their defaults are, \cr
#'  $alpha The tradeoff between l1 and l2 penalties in the elastic net. \cr
#'  $lambda The overall regularization parameter in the elastic net. \cr
#'  $nu The learning rate / step-size for the coordinate descent.
#'   Arbitrarily defaults to 1e-5. If things are not converging, try making this
#'   smaller. \cr
#'  $n_iter The number of steps of the alternating minimization to perform. \cr
#'  $n_iter_B The number of sweeps over all entries of B via coordinate descent.
#' @return The modified opts list, with defaults filled in for entries that had
#' been missing.
merge_spline_lf_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$alpha <- 1
  default_opts$lambda <- 0
  default_opts$n_iter <- 100
  default_opts$n_iter_B <- 10
  default_opts$nu <- 1e-5
  modifyList(default_opts, opts)
}

#' @title Compute terms in the objective for the spline latent factor model
spline_lf_objective <- function(Y, H, param) {
  c(sum((Y - H %*% param$B %*% t(param$W))^2, na.rm = T), sum(abs(param$W)),
    sum(param$W ^ 2), sum(abs(param$B)), sum(param$B ^ 2))
}

# optimization -----------------------------------------------------------------
#' @title Alternating minimization for spline latent source model
#' @description Minimize (over B and W)
#' ||M(Y - HBW')||_{2}^{2} + lambda * ((1 - alpha) ||W||_{2}^{2} + (alpha) * ||W||_{1}),
#' where M is a masking matrix that ignores terms in Y that are missing.
#' This is just the usual RSS with a basis matrix H plus an elastic net penalty
#' on W. The interpretation is that HB are latent sources and W are coefficients
#' on those latent sources. The sources are expanded in terms of a basis H, which
#' ensures that it is smooth / piecewise linear, ...
#' The masking matrix and the fact that H is not necessarily orthogonal is the
#' reason we have to optimize B using coordinate descent (and can't just
#' reformulate it as another elastic net problem, like in the lsam() function).
#' @param Y An N x J real matrix, with censored values set to NA. The standard
#' application we have in mind is the sample x OTU count matrix.
#' @param H The spline basis matrix, from which the latent sources arise (as
#' the linear mixture HB.)
#' @param W The samples-to-latent-source coefficients matrix, assumed known.
#' @param B0 The initial matrix for B in the optimization.
#' @param W0 The initial matrix for W in the optimization. Usually initialize
#' eigenvectors from SVD.
#' @param opts A (potentially only partially specified) list of options to use
#' during the optimization. See merge_spline_lf_opts() for choices.
#' gradient descent.
#' @return A list with the following elements, \cr
#'    $obj The matrix of terms of the objective function, after each iteration. \cr
#'     of the alternating minimization. \cr
#'    $W The optimized coefficients matrix. \cr
#'    $B The optimized mixing matrix for the sources.
#' @export
spline_lf <- function(Y, H, B0, W0, opts = list()) {
  opts <- merge_spline_lf_opts(opts)

  # get problem dimensions
  N <- nrow(Y)
  P <- ncol(Y)
  K <- ncol(W)
  L <- ncol(H)

  # initialize results
  param <- list(W = W0, B = B0)
  obj_names <- c("iter", "RSS", "W1", "W2", "B1", "B2")
  obj <- matrix(NA, opts$n_iter, 6, dimnames = list(1:opts$n_iter, obj_names))

  # perform optimization
  elnet_W <- elnet_fun(opts$lambda, opts$alpha)
  for(i in seq_len(opts$n_iter)) {
    param$W <- independent_models(H %*% param$B, Y, elnet_W)
    param$B <- get_B(Y, H, param$W, param$B, opts$nu, opts$n_iter_B)$B
    obj[i, ] <- c(i, spline_lf_objective(Y, H, param))
    cat(sprintf("Iteration %g | RSS: %f | W1: %f | W2: %f | | B1: %f | B2: %f\n",
                obj[i, 1], obj[i, 2], obj[i, 3], obj[i, 4], obj[i, 5], obj[i, 6]))
  }

  c(param, list(obj = obj))
}
