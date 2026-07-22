library(PTT)

rmvnorm.base <- function(n, mean, sigma) {
  z <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
  sweep(z %*% chol(sigma), 2, mean, FUN = "+")
}

dmvnorm.base <- function(x, mean, sigma) {
  centered <- sweep(as.matrix(x), 2, mean, FUN = "-")
  quadratic <- rowSums((centered %*% solve(sigma)) * centered)
  exp(-0.5 * quadratic) / (2 * pi * sqrt(det(sigma)))
}

set.seed(12345)
nobs <- 2000
p.global <- 0.85
mean.global <- c(0.5, 0.5)
sigma.global <- diag(c(0.01, 0.0064))
mean.local <- c(0.8, 0.2)
sigma.local <- sigma.global / 25
local.component <- runif(nobs) >= p.global

X <- cbind(rbeta(nobs, 20, 10), runif(nobs))
X[local.component, 1] <- rbeta(sum(local.component), 10, 20)
Y <- rmvnorm.base(nobs, mean.global, sigma.global)
local.response <- rmvnorm.base(nobs, mean.local, sigma.local)
Y[local.component, ] <- local.response[local.component, ]

y1.grid <- seq(0.0000001, 0.9999999, by = 0.01)
y2.grid <- seq(0.0000001, 0.9999999, by = 0.01)
y.grid <- expand.grid(y1 = y1.grid, y2 = y2.grid)
xpred <- c(0.8, 0.5)
x.grid <- matrix(rep(xpred, nrow(y.grid)), byrow = TRUE, ncol = 2)
max.resX <- 7

fit <- cond.opt(
  X = X, Y = Y, Xpred = x.grid, Ypred = y.grid,
  rho0.X = 0.5, rho0.Y = 0.5,
  max.resX = max.resX, max.resY = 7
)

global.mass <- p.global * dbeta(xpred[1], 20, 10)
local.mass <- (1 - p.global) * dbeta(xpred[1], 10, 20)
global.weight <- global.mass / (global.mass + local.mass)
true.density <- global.weight *
  dmvnorm.base(y.grid, mean.global, sigma.global) +
  (1 - global.weight) *
    dmvnorm.base(y.grid, mean.local, sigma.local)
true.surface <- matrix(
  true.density, nrow = length(y1.grid), ncol = length(y2.grid)
)
fitted.surface <- matrix(
  fit$predictive_densities,
  nrow = length(y1.grid), ncol = length(y2.grid)
)
surface.max <- max(true.surface, fitted.surface)
density.colors <- grDevices::hcl.colors(80, "Inferno")

draw.surface <- function(surface, title) {
  image(
    y1.grid, y2.grid, surface,
    col = density.colors, zlim = c(0, surface.max), useRaster = TRUE,
    xlab = "Y1", ylab = "Y2", main = title, asp = 1
  )
  contour(
    y1.grid, y2.grid, surface, add = TRUE, drawlabels = FALSE,
    col = grDevices::adjustcolor("white", alpha.f = 0.62), lwd = 0.7
  )
}

terminal <- fit$part_points_hmap
terminal <- terminal[terminal[, "state"] == Inf, , drop = FALSE]
xleft <- terminal[, "X1.l"] / 2^max.resX
xright <- (terminal[, "X1.u"] + 1) / 2^max.resX
ybottom <- terminal[, "X2.l"] / 2^max.resX
ytop <- (terminal[, "X2.u"] + 1) / 2^max.resX
level <- terminal[, "level"]
partition.colors <- grDevices::hcl.colors(max.resX + 1, "BluGrn")
component.colors <- grDevices::hcl.colors(2, "Dark 3")

old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(1, 3), mar = c(4, 4.1, 2.7, 1),
  oma = c(0, 0, 2, 0), las = 1
)
draw.surface(true.surface, "True conditional density")
draw.surface(fitted.surface, "PTT conditional density")

plot(
  c(0, 1), c(0, 1), type = "n",
  xlab = "X1", ylab = "X2", main = "hMAP predictor partition", asp = 1
)
rect(
  xleft, ybottom, xright, ytop,
  col = grDevices::adjustcolor(
    partition.colors[level + 1], alpha.f = 0.12
  ),
  border = grDevices::adjustcolor(
    partition.colors[level + 1], alpha.f = 0.72
  )
)
points(
  X, pch = ifelse(local.component, 17, 16), cex = 0.42,
  col = grDevices::adjustcolor(
    component.colors[1 + local.component], alpha.f = 0.36
  )
)
legend(
  "topleft", c("Global component", "Local component"),
  pch = c(16, 17), col = component.colors, bty = "n", cex = 0.78
)
mtext(
  sprintf(
    "Historical bivariate conditional OPT example at X = (%.1f, %.1f)",
    xpred[1], xpred[2]
  ),
  outer = TRUE
)
par(old.par)
