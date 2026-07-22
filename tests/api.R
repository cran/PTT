library(PTT)

set.seed(42)
x <- matrix(c(0.10, 0.20, 0.35, 0.65, 0.80, 0.90), ncol = 1)
pred <- matrix(c(0, 0.25, 0.50, 0.75, 1), ncol = 1)

fit <- apt(
  x, Xpred = pred, max.resol = 3,
  n.grid = 2, n.s = 2, n.post.samples = 1
)
stopifnot(
  length(fit$predictive_densities) == nrow(pred),
  all(is.finite(fit$predictive_densities)),
  all(fit$predictive_densities >= 0),
  identical(colnames(fit$part_points_hmap), c("X1.l", "X1.u", "level", "state")),
  identical(dim(fit$Omega), c(1L, 2L)),
  length(fit$part_points_post_samples) == 1L
)

set.seed(99)
sample.a <- apt(x, max.resol = 3, n.grid = 2, n.s = 2, n.post.samples = 1)
set.seed(99)
sample.b <- apt(x, max.resol = 3, n.grid = 2, n.s = 2, n.post.samples = 1)
stopifnot(identical(sample.a$part_points_post_samples, sample.b$part_points_post_samples))

kernel.a <- apt(x, max.resol = 2, tran.mode = 2, rho0 = 0.1, n.grid = 2, n.s = 2)
kernel.b <- apt(x, max.resol = 2, tran.mode = 2, rho0 = 0.9, n.grid = 2, n.s = 2)
stopifnot(
  identical(kernel.a$logrho, kernel.b$logrho),
  identical(kernel.a$logphi, kernel.b$logphi)
)

standardized <- apt(
  c(-3, -2, -1), Xpred = c(-3, -1),
  Omega.type = "standardized", max.resol = 2, n.grid = 2, n.s = 2
)
stopifnot(
  standardized$Omega[1, 1] < -3,
  standardized$Omega[1, 2] > -1,
  all(is.finite(standardized$predictive_densities))
)

conditional <- cond.apt(
  X = x, Y = rev(x), Xpred = pred[2:4, ], Ypred = pred[2:4, ],
  max.resX = 2, max.resY = 2, n.grid = 2, n.s = 3
)
stopifnot(
  length(conditional$predictive_densities) == 3L,
  all(is.finite(conditional$predictive_densities)),
  identical(dim(conditional$OmegaX), c(1L, 2L)),
  identical(dim(conditional$OmegaY), c(1L, 2L))
)

stopifnot(
  inherits(try(apt(c(-1, 0), max.resol = 2), silent = TRUE), "try-error"),
  inherits(try(apt(x, max.resol = 15), silent = TRUE), "try-error"),
  inherits(try(cond.apt(x, x, Xpred = 0.5), silent = TRUE), "try-error")
)
