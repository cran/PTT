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
nobs <- 3000
p.global <- 0.85
mean.global <- c(0.5, 0.5)
sigma.global <- diag(c(0.01, 0.0064))
mean.local <- c(0.8, 0.2)
sigma.local <- sigma.global / 25
local.component <- runif(nobs) >= p.global

global.sample <- rmvnorm.base(nobs, mean.global, sigma.global)
local.sample <- rmvnorm.base(nobs, mean.local, sigma.local)
X <- global.sample
X[local.component, ] <- local.sample[local.component, ]

x.grid <- seq(0.0000001, 0.9999999, by = 0.01)
y.grid <- seq(0.0000001, 0.9999999, by = 0.01)
prediction.grid <- expand.grid(x = x.grid, y = y.grid)
true.density <- p.global *
  dmvnorm.base(prediction.grid, mean.global, sigma.global) +
  (1 - p.global) *
    dmvnorm.base(prediction.grid, mean.local, sigma.local)

max.resol <- 10
n.post.samples <- 100
fit <- apt(
  X = X, Xpred = prediction.grid,
  max.resol = max.resol, rho0 = 0.2,
  tran.mode = 2, beta = 0,
  n.post.samples = n.post.samples
)

true.surface <- matrix(
  true.density, nrow = length(x.grid), ncol = length(y.grid)
)
fitted.surface <- matrix(
  fit$predictive_densities,
  nrow = length(x.grid), ncol = length(y.grid)
)
surface.max <- max(true.surface, fitted.surface)
density.colors <- grDevices::hcl.colors(80, "Inferno")
shown.observations <- sample(seq_len(nobs), 400)

draw.surface <- function(surface, title) {
  image(
    x.grid, y.grid, surface,
    col = density.colors, zlim = c(0, surface.max), useRaster = TRUE,
    xlab = "X1", ylab = "X2", main = title, asp = 1
  )
  contour(
    x.grid, y.grid, surface, add = TRUE, drawlabels = FALSE,
    col = grDevices::adjustcolor("white", alpha.f = 0.62), lwd = 0.7
  )
  points(
    X[shown.observations, ], pch = 16, cex = 0.3,
    col = grDevices::adjustcolor("white", alpha.f = 0.28)
  )
}

draw.partition <- function(sample.index) {
  partition <- fit$part_points_post_samples[[sample.index]]
  terminal <- partition[partition[, "nu"] == Inf, , drop = FALSE]
  density <- exp(terminal[, "logp"]) * 2^terminal[, "level"]
  shade <- log1p(density)
  shade.range <- range(shade)
  if (diff(shade.range) == 0) {
    color.index <- rep(40L, length(shade))
  } else {
    color.index <- 1L + floor(
      (shade - shade.range[1]) / diff(shade.range) * 79
    )
  }

  xleft <- terminal[, "X1.l"] / 2^max.resol
  xright <- (terminal[, "X1.u"] + 1) / 2^max.resol
  ybottom <- terminal[, "X2.l"] / 2^max.resol
  ytop <- (terminal[, "X2.u"] + 1) / 2^max.resol
  plot(
    c(0, 1), c(0, 1), type = "n", xlab = "X1", ylab = "X2",
    main = paste("Posterior partition sample", sample.index), asp = 1
  )
  rect(
    xleft, ybottom, xright, ytop,
    col = density.colors[color.index],
    border = grDevices::adjustcolor("white", alpha.f = 0.3)
  )
  points(
    X[shown.observations, ], pch = 16, cex = 0.28,
    col = grDevices::adjustcolor("white", alpha.f = 0.28)
  )
}

old.par <- par(no.readonly = TRUE)
par(
  mfrow = c(2, 2), mar = c(4, 4.2, 2.7, 1),
  oma = c(0, 0, 2, 0), las = 1
)
draw.surface(true.surface, "True mixture density")
draw.surface(fitted.surface, "Posterior predictive density")
draw.partition(25)
draw.partition(75)
mtext("Historical bivariate Markov APT example", outer = TRUE)
par(old.par)
