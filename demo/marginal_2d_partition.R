library(PTT)

rmvnorm.base <- function(n, mean, sigma) {
  z <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
  sweep(z %*% chol(sigma), 2, mean, FUN = "+")
}

set.seed(12345)
n <- 650
component <- rbinom(n, size = 1, prob = 0.35)
sample1 <- rmvnorm.base(
  n, c(-1.1, 0.25), matrix(c(0.55, 0.28, 0.28, 0.45), 2)
)
sample2 <- rmvnorm.base(
  n, c(1.25, 1.15), matrix(c(0.30, -0.14, -0.14, 0.38), 2)
)
X <- (1 - component) * sample1 + component * sample2

xlim <- range(X[, 1]) + c(-0.35, 0.35)
ylim <- range(X[, 2]) + c(-0.35, 0.35)
x.grid <- seq(xlim[1], xlim[2], length.out = 80)
y.grid <- seq(ylim[1], ylim[2], length.out = 75)
prediction.grid <- expand.grid(x = x.grid, y = y.grid)
max.resol <- 7

fit <- apt(
  X = X, Xpred = prediction.grid,
  Omega.type = "standardized", max.resol = max.resol,
  rho0 = 0.2, n.grid = 4, n.s = 4
)
density.surface <- matrix(
  fit$predictive_densities,
  nrow = length(x.grid), ncol = length(y.grid)
)

old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 1),
  oma = c(0, 0, 2, 0), las = 1
)

image(
  x.grid, y.grid, density.surface,
  col = grDevices::hcl.colors(80, "Inferno"), useRaster = TRUE,
  xlab = "X1", ylab = "X2", main = "Posterior predictive density",
  asp = 1
)
contour(
  x.grid, y.grid, density.surface, add = TRUE,
  drawlabels = FALSE,
  col = grDevices::adjustcolor("white", alpha.f = 0.65), lwd = 0.8
)
points(
  X, pch = 16, cex = 0.35,
  col = grDevices::adjustcolor("white", alpha.f = 0.28)
)

terminal <- fit$part_points_hmap
terminal <- terminal[terminal[, "state"] == Inf, , drop = FALSE]
omega <- fit$Omega
scale.x <- diff(omega[1, ]) / 2^max.resol
scale.y <- diff(omega[2, ]) / 2^max.resol
xleft <- omega[1, 1] + terminal[, "X1.l"] * scale.x
xright <- omega[1, 1] + (terminal[, "X1.u"] + 1) * scale.x
ybottom <- omega[2, 1] + terminal[, "X2.l"] * scale.y
ytop <- omega[2, 1] + (terminal[, "X2.u"] + 1) * scale.y
level <- terminal[, "level"]
partition.colors <- grDevices::hcl.colors(max.resol + 1, "BluGrn")

plot(
  X, pch = 16, cex = 0.55,
  col = grDevices::adjustcolor("grey20", alpha.f = 0.42),
  xlim = xlim, ylim = ylim, xlab = "X1", ylab = "X2",
  main = "hMAP terminal partition", asp = 1
)
rect(
  xleft, ybottom, xright, ytop,
  col = grDevices::adjustcolor(
    partition.colors[level + 1], alpha.f = 0.12
  ),
  border = grDevices::adjustcolor(
    partition.colors[level + 1], alpha.f = 0.75
  ),
  lwd = 0.9
)
points(
  X, pch = 16, cex = 0.42,
  col = grDevices::adjustcolor("grey20", alpha.f = 0.35)
)
mtext("Bivariate APT with a data-adaptive sample space", outer = TRUE)

par(old.par)
