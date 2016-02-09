
################################################################################
# The most elementary latent factor model: This minimizes
# ||Y - Phi W^{T}||_{2}^{2} over Phi and W through alternating least squares.
# Note that it only optimizes over nonmissing entries in Y.
################################################################################

# modeling ---------------------------------------------------------------------
#' @title Perform alternating minimization for a latent source model
#' @description This is minimizing ||Y - Phi * W^T||_{2}^{2} + lambda * ||W||_{1}
#' over W and Phi, by alternatively minimizing over each term. Note that Phi
#' is fixed across people for each time, which is why it is a shorter dimension
#' than Y. I.e., this model can also be written as
#' y_{ijt} = phi_{t} ^ T w_{j} + noise.
#' @param Y An n x J real matrix, with censored values set to NA. The standard
#' application we have in mind is the sample x OTU count matrix.
#' @param time_mask An n x T matrix, where T is the unique number of time points
#' across all samples, where the it^th entry is 1 is sample i was taken at time
#' i.
#' @param Phi0 The initial value of Phi, the basis matrix.
#' @param n_iter The number of iterations to run the alternating minimization.
#' @param model_funs A list of functions that, when given a matrix x and vector
#' y, will return the coefficients of x on y. Must have two elemets "Phi" and
#' "W" for fitting the two parts of the alternating minimization.
#' @return A list with the following elements, \cr
#'   $Phi The learned basis matrix across times. \cr
#'   $W The learned coefficients matrix. \cr
#'   $obj The RSS, and l1 / l2 of W across iterations.
#' @importFrom Matrix bdiag
#' @export
lsam <- function(Y, time_mask, Phi0, n_iter = 10, model_funs = NULL) {
  # checks / defaults
  if(ncol(time_mask) != nrow(Phi0)) {
    stop("Number of unique times in time mask must equal number of rows in Phi0")
  }
  if(is.null(model_funs)) {
    model_funs <- list(W = elnet_fun(.5, 1), Phi = elnet_fun(.5, 1))
  }

  # initialize results
  param <- list(Phi = time_mask %*% Phi0)
  obj <- matrix(NA, n_iter, 6,
                dimnames = list(1:n_iter, c("iter", "RSS", "W1", "W2", "Phi1", "Phi2")))

  # alternation
  for(i in seq_len(n_iter)) {
    param$W <- independent_models(param$Phi, Y, model_funs$W)
    param$Phi <- independent_models(param$W, t(Y), model_funs$Phi)

    # get objective
    obj[i, ] <- c(i, lsam_objective(Y, param))
    cat(sprintf("Iteration %g | RSS: %f | W1: %f | W2: %f | | Phi1: %f | Phi2: %f\n",
                obj[i, 1], obj[i, 2], obj[i, 3], obj[i, 4], obj[i, 5], obj[i, 6]))
  }

  param$Phi <- param$Phi[!duplicated(time_mask), ] # return to original unique times
  colnames(param$Phi) <- paste0("source", 1:ncol(Phi0))
  dimnames(param$W) <- list(colnames(Y), paste0("source_", 1:ncol(Phi0)))
  c(param, list(obj = obj))
}

#' @title Compute objective for alternating minimization
#' @param Y The data to which we aim to fit.
#' @param param A list with the elements Phi and W, for the source and
#' coefficients in the latent source model.
#' @return The RSS for the difference, and the l1 + l2 norms of the coefficients
#' @export
lsam_objective <- function(Y, param) {
  c(sum((Y - param$Phi %*% t(param$W))^2, na.rm = T), sum(abs(param$W)),
    sum(param$W ^ 2), sum(abs(param$Phi)), sum(param$Phi ^ 2))
}
