
################################################################################
# binary model with underlying splines
################################################################################

merge_spline_logit_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$n_iter <- 10
  default_opts$lambda_B <- 1e-3
  default_opts$alpha_B <- 0
  default_opts$lambda_W <- 1e-3
  default_opts$alpha_W <- 0
  modifyList(default_opts, opts)
}

lf_spline_logit <- function(Y, H, B0, W0, opts = list()) {
  opts <- merge_spline_logit_opts(opts)

  # get problem dimensions
  N <- nrow(Y)
  J <- ncol(Y)
  L <- ncol(H)
  K <- ncol(B0)

  # initialize results
  param <- list(B = B0, W = W0)
  obj_names <- c("iter", "ll", "W1", "W2", "B1", "B2", "obj")
  obj <- matrix(NA, opts$n_iter, length(obj_names),
                dimnames = list(1:opts$n_iter, obj_names))

  # perform optimization
  logit_W <- elnet_fun(opts$lambda_W, opts$alpha_W, family = "binomial")
  for(i in seq_len(opts$n_iter)) {
    param$W <- independent_models(H %*% param$B, Y, logit_W)
    param$B <- get_logit_B(Y, H, param$W, param$B)
    cur_obj <- lf_spline_logit_obj(Y, H, param$B, param$W, opts)

    obj[i, ] <- c(i, cur_obj)
    cat(sprintf("iter %g | ll: %f | W1: %f | W2: %f | B1: %f | B2: %f | obj %f \n",
                obj[i, 1], obj[i, 2], obj[i, 3], obj[i, 4],
                obj[i, 5], obj[i, 6], obj[i, 7]))
  }
  c(param, list(obj = obj))
}

lf_spline_logit_obj <- function(Y, H, B, W, opts) {
  ll <- mean(Y * H %*% B %*% t(W) - log(1 + exp(H %*% B %*% t(W))))
  W1 <- sum(abs(W))
  W2 <- sum(W ^ 2)
  B1 <- sum(abs(B))
  B2 <- sum(B ^ 2)
  obj <- ll - elnet_pen(opts$lambda_W, opts$alpha_W, W1, W2) -
    elnet_pen(opts$lambda_B, opts$alpha_B, B1, B2)
  c(ll, W1, W2, B1, B2, obj)
}

get_logit_B <- function(Y, H, W, B) {
  stop("get_logit_B() not implemented yet")
}
