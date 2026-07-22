library(PTT)

# These deterministic reference values guard the core marginal and conditional
# calculations against numerical drift while allowing platform-level
# floating-point noise.
tolerance <- 1e-9

x <- c(seq(0.05, 0.45, length.out = 24), seq(0.6, 0.95, length.out = 16))
grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)

apt.fit <- apt(
  x, Xpred = grid, max.resol = 4, rho0 = 0.2,
  n.grid = 3, n.s = 3
)
apt.expected <- c(
  -1.22884277881542, -0.380595133618687,
  1.04073167934391, 1.04446417461813, 0.946856014560911,
  0.960353263602841, 0.965065744499516
)
stopifnot(
  isTRUE(all.equal(
    unname(c(apt.fit$logrho, apt.fit$logphi, apt.fit$predictive_densities)),
    apt.expected,
    tolerance = tolerance
  )),
  identical(unname(apt.fit$part_points_hmap), matrix(c(0, 15, 0, Inf), nrow = 1L))
)

opt.fit <- opt(x, Xpred = grid, max.resol = 4, rho0 = 0.5)
opt.expected <- c(
  -0.0778733915972794, -0.615273788962671,
  1.01464670747005, 1.01474115491321, 0.983823971413429,
  0.98542445662195, 0.986027934553703
)
stopifnot(
  isTRUE(all.equal(
    unname(c(opt.fit$logrho, opt.fit$logphi, opt.fit$predictive_densities)),
    opt.expected,
    tolerance = tolerance
  )),
  identical(unname(opt.fit$part_points_hmap), matrix(c(0, 15, 0, Inf), nrow = 1L))
)

conditional.x <- seq(0.02, 0.98, length.out = 40)
conditional.y <- pmin(
  0.98,
  pmax(0.02, 0.15 + 0.7 * conditional.x + 0.08 * sin(8 * pi * conditional.x))
)
conditional.fit <- cond.opt(
  conditional.x, conditional.y,
  Xpred = c(0.2, 0.5, 0.8), Ypred = c(0.25, 0.5, 0.75),
  max.resX = 3, max.resY = 4
)
conditional.expected <- c(
  -28.8015601776622, 30.7668680880074,
  4.21400432146923, 4.44293777151951, 7.18750082136846
)
conditional.hmap.expected <- matrix(
  c(
    0, 0, 0, 2, 4, 4, 6,
    7, 3, 1, 3, 7, 5, 7,
    0, 1, 2, 2, 1, 2, 2,
    1, 1, Inf, Inf, 1, Inf, Inf,
    3.10199679410356e-13, 0.0718176406024455,
    0.985628998307472, 0.976837618064303,
    0.0718176406024452, 0.976837618064303, 0.985628998307472
  ),
  nrow = 7L
)
stopifnot(
  isTRUE(all.equal(
    unname(c(
      conditional.fit$logrho,
      conditional.fit$logphi,
      conditional.fit$predictive_densities
    )),
    conditional.expected,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    unname(conditional.fit$part_points_hmap),
    conditional.hmap.expected,
    tolerance = tolerance
  ))
)
