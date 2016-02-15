
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

# ---- basic-logit ----
svd_y <- svd(Y)
W0 <- svd_y$v[, 1:K]
Phi0 <- svd_y$u[, 1:K] %*% diag(svd_y$d[1:K])

lf_res <- logit_lf(Y, Phi0, W0)

P_res <- logit(lf_res$Phi %*% t(lf_res$W))
plot(P[, 1])
points(P_res[, 1], col = 'red')
