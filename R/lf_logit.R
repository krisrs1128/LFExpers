
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
  obj_names <- c("iter", "ll", "01", "10", "W1", "W2", "Phi1", "Phi2")
  obj <- matrix(NA, opts$n_iter, length(obj_names),
                dimnames = list(1:opts$n_iter, obj_names))

  # perform optimization
  logit_W <- elnet_fun(opts$lambda_W, opts$alpha_W, family = "binomial")
  logit_Phi <- elnet_fun(opts$lambda_Phi, opts$alpha_Phi, family = "binomial")
  for(i in seq_len(opts$n_iter)) {
    param$W <- independent_models(param$Phi, Y, logit_W)
    param$Phi <- independent_models(param$W, t(Y), logit_Phi)
    cat(sprintf("iter %f", i))
    #obj[i, ] <- c(i, lf_logit_obj)
  }
  c(param, list(obj = obj))

}

lf_logit_obj <- function(Y, Phi, W) {

}
