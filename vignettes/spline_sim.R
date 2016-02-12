
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

library("splines")
H <- bs(1:N, df = L, degree = 1)
W <- matrix(rnorm(P * K), P, K)
B <- matrix(rnorm(L * K), L, K)
B0 <- matrix(rnorm(L * K), L, K)
W0 <- matrix(rnorm(P * K), P, K)

sigmas <- 10 ^ seq(-1, 1, length.out = 10)
na_props <- seq(0, .95, length.out = 10)
fits <- vector(length = 100, mode = "list")

ix <- 1
for(i in seq_along(sigmas)) {
  for(j in seq_along(na_props)) {
    E <- matrix(sigmas[i] * rnorm(N * P), N, P)
    Y <- H %*% B %*% t(W) + E
    Y[sample(N * P, N * P * na_props[j])] <- NA # 40% missing at random

    # fit the model
    fits[[ix]] <- spline_lf(Y, H, B0, W0, list(n_iter = 50))
    ix <- ix + 1
  }
}
