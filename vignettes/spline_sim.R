
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
