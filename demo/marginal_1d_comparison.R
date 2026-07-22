library(PTT)

set.seed(12345)
n <- 450
component <- rbinom(n, size = 1, prob = 0.35)
x <- numeric(n)
x[component == 0] <- rbeta(sum(component == 0), 3, 9)
x[component == 1] <- rbeta(sum(component == 1), 12, 4)

x.grid <- seq(0.001, 0.999, length.out = 500)
true.density <- 0.65 * dbeta(x.grid, 3, 9) +
  0.35 * dbeta(x.grid, 12, 4)
max.resol <- 8

apt.fit <- apt(
  x, Xpred = x.grid, max.resol = max.resol,
  rho0 = 0.2, n.grid = 5, n.s = 5
)
opt.fit <- opt(
  x, Xpred = x.grid, max.resol = max.resol, rho0 = 0.5
)

series.colors <- grDevices::hcl.colors(3, "Dark 3")
old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(2, 1), mar = c(4, 4.2, 2.5, 1),
  oma = c(0, 0, 2, 0), las = 1
)

hist(
  x, breaks = 28, probability = TRUE,
  col = grDevices::adjustcolor(series.colors[1], alpha.f = 0.18),
  border = "white", xlim = c(0, 1),
  xlab = "x", ylab = "Density", main = "Posterior predictive densities"
)
lines(x.grid, true.density, col = series.colors[1], lwd = 3)
lines(
  x.grid, apt.fit$predictive_densities,
  col = series.colors[2], lwd = 2.5
)
lines(
  x.grid, opt.fit$predictive_densities,
  col = series.colors[3], lwd = 2.5, lty = 2
)
rug(x, col = grDevices::adjustcolor("grey20", alpha.f = 0.22))
legend(
  "topright", c("True density", "APT", "OPT"),
  col = series.colors, lty = c(1, 1, 2), lwd = c(3, 2.5, 2.5),
  bty = "n"
)

terminal <- apt.fit$part_points_hmap
terminal <- terminal[terminal[, "state"] == Inf, , drop = FALSE]
left <- terminal[, "X1.l"] / 2^max.resol
right <- (terminal[, "X1.u"] + 1) / 2^max.resol
level <- terminal[, "level"]
level.colors <- grDevices::hcl.colors(max.resol + 1, "BluGrn")

plot(
  c(0, 1), range(c(1, level)), type = "n", bty = "n",
  xlab = "Sample space", ylab = "Terminal-cell resolution",
  main = "Adaptive terminal partition", yaxt = "n"
)
axis(2, at = sort(unique(level)))
rect(
  left, level - 0.32, right, level + 0.32,
  col = level.colors[level + 1], border = "white"
)
rug(x, col = grDevices::adjustcolor("grey20", alpha.f = 0.18))
mtext("Univariate density estimation and multiscale adaptation", outer = TRUE)

par(old.par)
