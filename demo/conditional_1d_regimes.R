library(PTT)

set.seed(12345)
size <- 500
size.pred <- 100
x <- rbeta(size + size.pred, 2, 2)
y <- rbeta(size + size.pred, 10, 30)
y[x < 0.25] <- rbeta(sum(x < 0.25), 30, 20)
y[x > 0.5] <- rbeta(sum(x > 0.5), 0.5, 0.5)

X <- x[seq_len(size)]
Y <- y[seq_len(size)]
Xtest <- x[size + seq_len(size.pred)]
Ytest <- y[size + seq_len(size.pred)]

slice.locations <- c(0.12, 0.38, 0.75)
y.grid <- seq(0.001, 0.999, length.out = 300)
prediction.x <- c(
  Xtest,
  rep(slice.locations, each = length(y.grid))
)
prediction.y <- c(
  Ytest,
  rep(y.grid, times = length(slice.locations))
)

fit <- cond.opt(
  X = X, Y = Y,
  Xpred = prediction.x, Ypred = prediction.y,
  max.resX = 10, max.resY = 10,
  rho0.X = 0.5, rho0.Y = 0.5
)
test.density <- fit$predictive_densities[seq_len(size.pred)]
log.predictive.score <- sum(log(pmax(test.density, .Machine$double.xmin)))
slice.density <- matrix(
  fit$predictive_densities[-seq_len(size.pred)],
  nrow = length(y.grid), ncol = length(slice.locations)
)

print(fit$part_points_hmap)
print(log.predictive.score)

regime.colors <- grDevices::hcl.colors(3, "Dark 3")
old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(2, 2), mar = c(4, 4.2, 2.6, 1),
  oma = c(0, 0, 2, 0), las = 1
)

plot(
  c(0, 1), c(0, 1), type = "n",
  xlab = "Predictor x", ylab = "Response y",
  main = "Piecewise conditional regimes"
)
rect(
  c(0, 0.25, 0.5), 0, c(0.25, 0.5, 1), 1,
  col = grDevices::adjustcolor(regime.colors, alpha.f = 0.09),
  border = NA
)
abline(v = c(0.25, 0.5), col = "grey60", lty = 3)
points(
  X, Y, pch = 16, cex = 0.55,
  col = grDevices::adjustcolor("grey20", alpha.f = 0.45)
)

for (i in seq_along(slice.locations)) {
  x0 <- slice.locations[i]
  true.density <- if (x0 < 0.25) {
    dbeta(y.grid, 30, 20)
  } else if (x0 > 0.5) {
    dbeta(y.grid, 0.5, 0.5)
  } else {
    dbeta(y.grid, 10, 30)
  }
  fitted.density <- slice.density[, i]
  ymax <- 1.05 * max(true.density, fitted.density)

  plot(
    y.grid, fitted.density, type = "n", ylim = c(0, ymax),
    xlab = "Response y", ylab = "Density",
    main = sprintf("Conditional density at x = %.2f", x0)
  )
  polygon(
    c(y.grid, rev(y.grid)),
    c(rep(0, length(y.grid)), rev(fitted.density)),
    col = grDevices::adjustcolor(regime.colors[i], alpha.f = 0.16),
    border = NA
  )
  lines(y.grid, fitted.density, col = regime.colors[i], lwd = 2.5)
  lines(y.grid, true.density, col = "grey25", lwd = 2, lty = 2)
  rug(
    Y[abs(X - x0) < 0.055],
    col = grDevices::adjustcolor("grey20", alpha.f = 0.28)
  )
  if (i == 1) {
    legend(
      "topright", c("PTT estimate", "True density"),
      col = c(regime.colors[i], "grey25"), lty = c(1, 2),
      lwd = c(2.5, 2), bty = "n"
    )
  }
}
mtext(
  sprintf(
    "Historical conditional OPT example - held-out log score %.1f",
    log.predictive.score
  ),
  outer = TRUE
)
par(old.par)
