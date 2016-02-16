
################################################################################
# Latent factor modeling for binary data, fit using alternating logistic
# regressions.
################################################################################

merge_logit_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$n_iter <- 10
  default_opts$lambda_Phi <- 1e-3
  default_opts$alpha_Phi <- 0
  default_opts$lambda_W <- 1e-3
  default_opts$alpha_W <- 0
  modifyList(default_opts, opts)
}

logit_lf <- function(Y, Phi0, W0, opts = list()) {
  opts <- merge_logit_opts(opts)

  # get problem dimensions
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Phi0)

  # initialize results
  param <- list(Phi = Phi0, W = W0)
  obj_names <- c("iter", "ll", "W1", "W2", "Phi1", "Phi2", "obj")
  obj <- matrix(NA, opts$n_iter, length(obj_names),
                dimnames = list(1:opts$n_iter, obj_names))

  # perform optimization
  logit_W <- logit_fun(opts$lambda_W, opts$alpha_W)
  logit_Phi <- logit_fun(opts$lambda_Phi, opts$alpha_Phi)
  for(i in seq_len(opts$n_iter)) {
    param$W <- independent_models(param$Phi, Y, logit_W)
    param$Phi <- independent_models(param$W, t(Y), logit_Phi)
    cur_obj <- lf_logit_obj(Y, param$Phi, param$W, opts)

    obj[i, ] <- c(i, cur_obj)
    cat(sprintf("iter %g | ll: %f | W1: %f | W2: %f | Phi1: %f | Phi2: %f | obj %f \n",
                obj[i, 1], obj[i, 2], obj[i, 3], obj[i, 4],
                obj[i, 5], obj[i, 6], obj[i, 7]))
  }
  c(param, list(obj = obj))
}

lf_logit_obj <- function(Y, Phi, W, opts) {
  ll <- mean(Y * Phi %*% t(W) - log(1 + exp(Phi %*% t(W))))
  W1 <- sum(abs(W))
  W2 <- sum(W ^ 2)
  Phi1 <- sum(abs(Phi))
  Phi2 <- sum(Phi ^ 2)
  obj <- ll - elnet_pen(opts$lambda_W, opts$alpha_W, W1, W2) -
    elnet_pen(opts$lambda_Phi, opts$alpha_Phi, Phi1, Phi2)
  c(ll, W1, W2, Phi1, Phi2, obj)
}
