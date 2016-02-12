
################################################################################
# R script accompanying spline_sim.Rmd
################################################################################

## ---- libraries ----
library("LFExpers")
library("reshape2")
library("plyr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("splines")
library("stringr")
theme_set(theme_bw())

## ---- example-reconstruction ----
example(spline_lf)

Y_hat <- H %*% spline_fit$B %*% t(spline_fit$W)
mY <- rbind(data.frame(type = "truth", melt(Y)),
            data.frame(type = "fit", melt(Y_hat)))
colnames(mY) <- c("type", "time", "response", "value")
ggplot(mY) +
  geom_point(aes(x = time, y = value, col = type), size = .5, alpha = 0.8) +
  facet_wrap(~ response) +
  ggtitle("True vs. fitted Y_{j}'s")

## ---- example-sources ----
mHB <- rbind(data.frame(type = "truth", melt(H %*% B)),
             data.frame(type = "fit", melt(H %*% spline_fit$B)))
colnames(mHB) <- c("type", "time", "source", "value")

ggplot(mHB) +
  geom_line(aes(x = time, y = value, group = source), size = .5, alpha = 0.8) +
  facet_wrap(~type) +
  ggtitle("True vs. fitted sources (unidentifiable)")

## ---- simulation-data ----
N <- 150
P <- 20
K <- 5
L <- 6

H <- bs(1:N, df = L, degree = 1)
W <- matrix(rnorm(P * K), P, K)
B <- matrix(rnorm(L * K), L, K)
B0 <- matrix(rnorm(L * K), L, K)
W0 <- matrix(rnorm(P * K), P, K)

sigmas <- 10 ^ seq(-1, 1, length.out = 10)
na_props <- seq(0, .9, length.out = 10)
fits <- vector(length = 100, mode = "list")

ix <- 1
for(i in seq_along(sigmas)) {
  for(j in seq_along(na_props)) {
    E <- matrix(sigmas[i] * rnorm(N * P), N, P)
    Y <- H %*% B %*% t(W) + E
    Y[sample(N * P, N * P * na_props[j])] <- NA

    # fit the model
    fits[[ix]] <- spline_lf(Y, H, B0, W0, list(n_iter = 50))
    fits[[ix]]$param <- list("sigma" = sigmas[i], "na_prop" = na_props[j])
    ix <- ix + 1
  }
}

## ---- sim-funs ----
obj_line <- function(rss, w1, w2, b1, b2, lambdas, alphas) {
  rss +
    elnet_pen(lambdas[1], alphas[1], w1, w2) +
    elnet_pen(lambdas[2], alphas[2], b1, b2)
}

obj_matrix <- function(X, lambdas, alphas) {
  n <- nrow(X)
  res <- vector(length = n)
  for(i in seq_len(n)) {
    res[i] <- obj_line(X[i, "RSS"], X[i, "W1"], X[i, "W2"], X[i, "B1"],
                       X[i, "B2"], lambdas, alphas)
  }
  res
}

## ---- process-sims ----
objs <- lapply(fits, function(x) {
  obj_matrix(x$obj, c(1, 0), c(0, 1))
})

all_rss <-  lapply(fits, function(x) {
  obj <- x$obj
  obj[nrow(obj), "RSS"]
})

params <- lapply(fits, function(x) {
  colwise(function(x) { round(x,  3)})(data.frame(x$param))
})
params <- do.call(rbind, params)

params_vs_obj <- data.frame(params, obj = do.call(rbind, objs)) %>%
  melt(id.vars = c("sigma", "na_prop"))
params_vs_obj$iter <- str_extract(params_vs_obj$variable, "[0-9]+") %>%
  as.numeric()

## ---- plot-sim ----
ggplot(params_vs_obj) +
  geom_line(aes(x = iter, y = value, col = sigma, group = sigma)) +
  facet_wrap(~na_prop, scale = "free_y") +
  scale_y_log10() +
  ggtitle("Objective over iterations, across regimes")

params_vs_rss <- data.frame(params, rss = do.call(rbind, all_rss))
ggplot(params_vs_rss) +
  geom_tile(aes(x = as.factor(sigma), y = as.factor(na_prop),
                fill = log(rss))) +
  ggtitle("Reconstruction RSS across regimes")
