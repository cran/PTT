library(PTT)

set.seed(12345)
n <- 600
x <- runif(n)
shape1 <- 2 + 9 * x
shape2 <- 11 - 7 * x
y <- rbeta(n, shape1, shape2)

x.grid <- seq(0.01, 0.99, length.out = 75)
y.grid <- seq(0.01, 0.99, length.out = 100)
prediction.grid <- expand.grid(x = x.grid, y = y.grid)

fit <- cond.apt(
  X = x, Y = y,
  Xpred = prediction.grid$x, Ypred = prediction.grid$y,
  max.resX = 5, max.resY = 7,
  rho0.X = 0.25, rho0.Y = 0.2,
  n.grid = 4, n.s = 4
)
density.surface <- matrix(
  fit$predictive_densities,
  nrow = length(x.grid), ncol = length(y.grid)
)

old.par <- par(no.readonly = TRUE)
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))
par(oma = c(0, 0, 2, 0), las = 1)

par(mar = c(4.2, 4.2, 2.4, 1))
image(
  x.grid, y.grid, density.surface,
  col = grDevices::hcl.colors(80, "Inferno"), useRaster = TRUE,
  xlab = "Predictor x", ylab = "Response y",
  main = "Estimated conditional density"
)
contour(
  x.grid, y.grid, density.surface, add = TRUE,
  drawlabels = FALSE,
  col = grDevices::adjustcolor("white", alpha.f = 0.55), lwd = 0.7
)
points(
  x, y, pch = 16, cex = 0.35,
  col = grDevices::adjustcolor("white", alpha.f = 0.38)
)

slice.locations <- c(0.2, 0.8)
slice.colors <- grDevices::hcl.colors(2, "Dark 3")
for (i in seq_along(slice.locations)) {
  x0 <- slice.locations[i]
  index <- which.min(abs(x.grid - x0))
  fitted.slice <- density.surface[index, ]
  true.slice <- dbeta(y.grid, 2 + 9 * x0, 11 - 7 * x0)
  ymax <- 1.08 * max(fitted.slice, true.slice)

  par(mar = c(4.2, 4.2, 2.4, 1))
  plot(
    y.grid, fitted.slice, type = "n", ylim = c(0, ymax),
    xlab = "Response y", ylab = "Density",
    main = sprintf("Slice at x = %.1f", x0)
  )
  polygon(
    c(y.grid, rev(y.grid)),
    c(rep(0, length(y.grid)), rev(fitted.slice)),
    col = grDevices::adjustcolor(slice.colors[i], alpha.f = 0.16),
    border = NA
  )
  lines(y.grid, fitted.slice, col = slice.colors[i], lwd = 2.5)
  lines(y.grid, true.slice, col = "grey25", lwd = 2, lty = 2)
  nearby <- abs(x - x0) < 0.05
  rug(y[nearby], col = grDevices::adjustcolor("grey20", alpha.f = 0.3))
  if (i == 1) {
    legend(
      "topright", c("PTT estimate", "True density"),
      col = c(slice.colors[i], "grey25"), lty = c(1, 2),
      lwd = c(2.5, 2), bty = "n"
    )
  }
}
mtext("Conditional density: surface and representative slices", outer = TRUE)

layout(1)
par(old.par)
