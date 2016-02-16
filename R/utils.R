
################################################################################
# Utilities for general modeling problems.
################################################################################

# logit ------------------------------------------------------------------------
logit <- function(x) {
  exp(x) / (1 + exp(x))
}

# elnet ------------------------------------------------------------------------
#' @title Wrapper for glmnet at a specific lambda
#' @param lambda What lambda to use in the optimization?
#' @return A function that fits the elastic net given only two arguments (x and y)
#' at the specified lambda.
#' @export
#' @importFrom glmnet glmnet
elnet_fun <- function(lambda, alpha, family = "gaussian") {
  function(x, y) {
    B <- glmnet(x, y, intercept = F, family = family,
                lambda = lambda, alpha = alpha)$beta
    B[, 1]
  }
}

logit_fun <- function(lambda, alpha){
  function(x, y) {
    if(length(unique(y)) > 1 & all(table(y) > 1)) {
      B <- glmnet(x, y, intercept = F, family = family,
                  lambda = lambda, alpha = alpha)$beta
    } else {
      B <- matrix(0, ncol(x), 1)
    }
    B[, 1]
  }
}

#' @title Compute the elastic net penalty
#' @export
elnet_pen <- function(lambda, alpha, w1, w2) {
  lambda * ((1 - alpha) / 2 * w2 + alpha * w1)
}

# indep ------------------------------------------------------------------------
#' @title Fit an independent model for each response coordinate
#' @param X The N x P covariates matrix. NA values will be dropped.
#' @param Y an N x J response matrix
#' @param model_fun A function that returns coefficients Y[, j] ~ X for each j.
#' @return The J X P coefficients matrix.
#' @export
independent_models <- function(X, Y, model_fun) {
  B <- matrix(nrow = ncol(Y), ncol = ncol(X))
  for(j in seq_len(ncol(Y))) {
    non_na <- !is.na(Y[, j])
    B[j, ] <- model_fun(X[non_na, ], Y[non_na, j])
  }
  B
}
