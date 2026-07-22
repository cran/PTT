library(PTT)

set.seed(12345)
nobs <- 5000
max.resol <- 11
p <- c(0.1, 0.3, 0.4, 0.2)
sample.indicator <- runif(nobs)

X <- as.matrix(
  (sample.indicator <= p[1]) * runif(nobs) +
    (p[1] < sample.indicator & sample.indicator <= sum(p[1:2])) *
      (0.25 + rbeta(nobs, 1, 1) * 0.25) +
    (sum(p[1:2]) < sample.indicator &
       sample.indicator <= sum(p[1:3])) *
      (0.25 + rbeta(nobs, 2, 2) * 0.25) +
    (sample.indicator > sum(p[1:3])) * rbeta(nobs, 5000, 2000)
)

x.grid <- seq(0.0001, 0.9999, by = 0.002)
true.density <- p[1] * dunif(x.grid) +
  p[2] * dunif(x.grid, 0.25, 0.5) +
  p[3] * 4 * dbeta((x.grid - 0.25) * 4, 2, 2) +
  p[4] * dbeta(x.grid, 5000, 2000)

apt.fit <- apt(
  X = X, Xpred = x.grid, max.resol = max.resol, rho0 = 0.2
)
opt.fit <- opt(
  X = X, Xpred = x.grid, max.resol = max.resol, rho0 = 0.2
)

series.colors <- grDevices::hcl.colors(3, "Dark 3")
ymax <- 1.05 * max(
  true.density,
  apt.fit$predictive_densities,
  opt.fit$predictive_densities
)
draw.fit <- function(fitted.density, fitted.color, title) {
  hist(
    X[, 1], breaks = 60, probability = TRUE,
    col = grDevices::adjustcolor(series.colors[1], alpha.f = 0.15),
    border = "white", xlim = c(0, 1), ylim = c(0, ymax),
    xlab = "x", ylab = "Density", main = title
  )
  lines(x.grid, true.density, col = series.colors[1], lwd = 2.8)
  lines(x.grid, fitted.density, col = fitted.color, lwd = 2.4)
  rug(X[, 1], col = grDevices::adjustcolor("grey20", alpha.f = 0.12))
  legend(
    "topleft", c("True density", "PTT estimate"),
    col = c(series.colors[1], fitted.color), lty = 1,
    lwd = c(2.8, 2.4), bty = "n"
  )
}

old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(2, 1), mar = c(4, 4.2, 2.4, 1),
  oma = c(0, 0, 2, 0), las = 1
)
draw.fit(
  apt.fit$predictive_densities, series.colors[2],
  "Markov APT posterior predictive density"
)
draw.fit(
  opt.fit$predictive_densities, series.colors[3],
  "OPT posterior predictive density"
)
mtext("Historical univariate mixture example", outer = TRUE)
par(old.par)
