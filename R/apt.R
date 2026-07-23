#' Fit an Adaptive Pólya Tree
#'
#' Fits an adaptive Pólya tree (APT) or Markov APT to univariate or
#' multivariate observations. The partition recursively bisects the sample
#' space, while latent states control local shrinkage toward the uniform
#' centering distribution.
#'
#' @param X A numeric vector, matrix, or data frame. Rows are observations and
#'   columns are dimensions.
#' @param Xpred Optional prediction locations. A vector is interpreted as many
#'   locations for univariate data and as one location for multivariate data.
#'   If `NULL`, no predictive densities are computed.
#' @param Omega.type Sample-space construction: `"unit"` uses the unit
#'   hypercube and `"standardized"` uses a slightly padded empirical range of
#'   `X` and `Xpred` in each dimension.
#' @param max.resol Maximum partition depth. Must be an integer from 1 to 14.
#' @param rho0 Baseline prior probability of complete shrinkage. Must lie in
#'   `[0, 1]`.
#' @param rho0.mode How the complete-shrinkage probability changes with level:
#'   `0` keeps it constant; `1` uses
#'   `1 - (1 - rho0) / (level + 1)^2`; and `2` uses
#'   `1 - (1 - rho0) / 2^level`.
#' @param tran.mode Shrinkage-state transition model. `0` makes non-stopping
#'   states independent across nodes; `1` imposes stochastically increasing
#'   shrinkage and treats complete shrinkage separately; `2` applies the
#'   increasing transition kernel to all states, including complete shrinkage.
#' @param lognu.lb,lognu.ub Lower and upper support bounds for `log10(nu)`,
#'   where `nu` is the beta precision (shrinkage) parameter.
#' @param n.grid Number of midpoint quadrature points within each shrinkage
#'   state.
#' @param n.s Number of non-stopping shrinkage states. Must be an integer from
#'   1 to 65,534.
#' @param beta Nonnegative stickiness parameter in the exponential transition
#'   kernel. `beta = 0` gives uniform probability over admissible higher states.
#' @param n.post.samples Number of posterior partition samples to draw.
#'
#' @returns A list with `logrho` (the log posterior complete-shrinkage
#'   probability at the root), `logphi` (the log marginal likelihood),
#'   `part_points_hmap` (the hierarchical MAP partition),
#'   `predictive_densities`, and `Omega`. If posterior samples are requested,
#'   `part_points_post_samples` is also returned. Partition bounds are integer
#'   cell indices on the depth-`max.resol` grid; `Omega` records the physical
#'   sample-space bounds.
#'
#' @references
#' Ma, L. (2017). Adaptive shrinkage in Pólya tree type models.
#' *Bayesian Analysis*, 12(3), 779--805. \doi{10.1214/16-BA1021}.
#'
#' @seealso [opt()], [cond.apt()]
#' @export
#'
#' @examples
#' set.seed(12345)
#' x <- c(rbeta(40, 3, 8), rbeta(40, 9, 3))
#' grid <- seq(0.05, 0.95, length.out = 25)
#' fit <- apt(x, Xpred = grid, max.resol = 4, n.grid = 3, n.s = 3)
#' head(fit$predictive_densities)
apt <- function(X, Xpred = NULL, Omega.type = "unit", max.resol = 10,
                rho0 = 0.2, rho0.mode = 0, tran.mode = 1,
                lognu.lb = -1, lognu.ub = 4, n.grid = 5, n.s = 5,
                beta = 0.1, n.post.samples = 0) {
  X <- .as_ptt_matrix(X, "X")
  p <- ncol(X)
  Xpred <- if (is.null(Xpred)) {
    .empty_prediction(p)
  } else {
    .as_ptt_matrix(Xpred, "Xpred", p, allow_empty = TRUE)
  }

  Omega.type <- match.arg(Omega.type, c("unit", "standardized"))
  max.resol <- .integer_arg(max.resol, "max.resol", 1L, 14L)
  .validate_tree_capacity(p, max.resol, "X")
  rho0 <- .numeric_arg(rho0, "rho0", 0, 1)
  rho0.mode <- .integer_arg(rho0.mode, "rho0.mode", 0L, 2L)
  tran.mode <- .integer_arg(tran.mode, "tran.mode", 0L, 2L)
  lognu.lb <- .numeric_arg(lognu.lb, "lognu.lb")
  lognu.ub <- .numeric_arg(lognu.ub, "lognu.ub")
  if (lognu.lb > lognu.ub) {
    stop("`lognu.lb` must not exceed `lognu.ub`.", call. = FALSE)
  }
  n.grid <- .integer_arg(n.grid, "n.grid", 1L)
  n.s <- .integer_arg(n.s, "n.s", 1L, 65534L)
  beta <- .numeric_arg(beta, "beta", 0)
  n.post.samples <- .integer_arg(n.post.samples, "n.post.samples", 0L)

  Omega <- .make_omega(X, Xpred, Omega.type, "X")
  ans <- fitPTTcpp(
    X, Xpred, Omega, max.resol, rho0, rho0.mode, tran.mode,
    lognu.lb, lognu.ub, n.grid, n.s, beta, n.post.samples
  )

  .format_marginal_fit(ans, p, n.s, n.post.samples, Omega)
}

#' Fit an Optional Pólya Tree
#'
#' A convenience wrapper around [apt()] that uses one non-stopping state with
#' zero beta precision, reproducing the optional Pólya tree (OPT).
#'
#' @inheritParams apt
#' @returns A list with the same structure as [apt()].
#'
#' @references
#' Wong, W. H. and Ma, L. (2010). Optional Pólya tree and Bayesian inference.
#' *The Annals of Statistics*, 38(3), 1433--1459.
#' \doi{10.1214/09-AOS755}.
#'
#' @seealso [apt()], [cond.opt()]
#' @export
#'
#' @examples
#' set.seed(12345)
#' x <- rbeta(60, 4, 7)
#' fit <- opt(x, Xpred = c(0.25, 0.5, 0.75), max.resol = 4)
#' fit$predictive_densities
opt <- function(X, Xpred = NULL, Omega.type = "unit", max.resol = 10,
                rho0 = 0.5, rho0.mode = 0, n.post.samples = 0) {
  apt(
    X = X, Xpred = Xpred, Omega.type = Omega.type,
    max.resol = max.resol, rho0 = rho0, rho0.mode = rho0.mode,
    tran.mode = 1, lognu.lb = 0, lognu.ub = 0,
    n.grid = 1, n.s = 1, beta = 0,
    n.post.samples = n.post.samples
  )
}
