
################################################################################
# Simulation to recover latent factors in binary data
################################################################################

## ---- libraries ----
library("splines")

## ---- generate-data ----
N <- 100
K <- 6
L <- 3
J <- 1000

H <- bs(1:N, df = L, degree = 1)
B <- matrix(rnorm(L * K), L, K)
Phi <- H %*% B
W <- matrix(rnorm(J * K), J, K)

P <- logit(Phi %*% t(W))
Y <- matrix(0, N, J)

for(i in 1:nrow(P)) {
  for(j in 1:ncol(P)) {
    Y[i, j] <- rbinom(1, 1, P[i, j])
  }
}

plot(Phi[, 1])
plot(P[, 1])
plot(Y[, 1])
plot(P[, 1], Y[, 1])

## ---- basic-logit ----
svd_y <- svd(Y)
W0 <- svd_y$v[, 1:K]
Phi0 <- svd_y$u[, 1:K] %*% diag(svd_y$d[1:K])

opts <- list(alpha_W = 0.7, lambda_W = 1e-2,
             alpha_Phi = 0.7, lambda_Phi = 1e-2)
lf_res <- logit_lf(Y, Phi0, W0, opts)

plot(lf_res$obj)

P_hat <- logit(lf_res$Phi %*% t(lf_res$W))
plot(P[, 1])
points(P_hat[, 1], col = 'red')

## ---- spline-logit ----
B0 <- matrix(rnorm(L * K), L, K)
spline_lf_res <- lf_spline_logit(Y, H, B0, W0, list(n_iter_B = 50, eta = 5e-5, n_iter = 30))
plot(spline_lf_res$obj[, "obj"])

P_hat <- logit(H %*% spline_lf_res$B %*% t(spline_lf_res$W))
for(j in seq_len(ncol(P))) {
  plot(P[, j], ylim = c(0, 1))
  points(P_hat[, j], col = "red")
  Sys.sleep(.2)
}

