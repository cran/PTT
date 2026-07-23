#' Fit a Conditional Adaptive Pólya Tree
#'
#' Fits the two-stage conditional-density model of Ma (2017): the predictor
#' space is recursively partitioned and an adaptive Pólya tree models the
#' response distribution within each predictor block.
#'
#' @param X,Y Numeric vectors, matrices, or data frames containing paired
#'   predictor and response observations. They must have the same number of
#'   rows.
#' @param Xpred,Ypred Optional paired predictor and response locations at which
#'   to evaluate the conditional predictive density. Supply both or neither;
#'   they must have the same number of rows.
#' @param OmegaX.type,OmegaY.type Sample-space construction for predictors and
#'   responses. See `Omega.type` in [apt()].
#' @param max.resX,max.resY Maximum partition depths for predictor and response
#'   spaces, each between 1 and 14.
#' @param rho0.X,rho0.Y Baseline prior complete-shrinkage probabilities for the
#'   predictor partition and response model.
#' @param rho0.mode.X,rho0.mode.Y Level-dependent stopping modes for the
#'   predictor partition and response model. See `rho0.mode` in [apt()].
#' @inheritParams apt
#'
#' @returns A list containing `logrho`, `logphi`, the predictor-space
#'   `part_points_hmap`, `predictive_densities`, `OmegaX`, and `OmegaY`. If
#'   requested, `part_points_post_samples` contains sampled predictor
#'   partitions. Partition bounds use integer cell indices at depth `max.resX`.
#'
#' @references
#' Ma, L. (2017). Recursive partitioning and multi-scale modeling on
#' conditional densities. *Electronic Journal of Statistics*, 11(1),
#' 1297--1325. \doi{10.1214/17-EJS1254}.
#'
#' @seealso [cond.opt()], [apt()]
#' @export
#'
#' @examples
#' set.seed(12345)
#' x <- runif(80)
#' y <- rbeta(80, 2 + 6 * x, 8 - 5 * x)
#' fit <- cond.apt(
#'   x, y, Xpred = c(0.25, 0.75), Ypred = c(0.3, 0.7),
#'   max.resX = 3, max.resY = 4, n.grid = 3, n.s = 3
#' )
#' fit$predictive_densities
cond.apt <- function(X, Y, Xpred = NULL, Ypred = NULL,
                     OmegaX.type = "unit", OmegaY.type = "unit",
                     max.resX = 5, max.resY = 8,
                     rho0.X = 0.2, rho0.mode.X = 0,
                     rho0.Y = 0.2, rho0.mode.Y = 0,
                     tran.mode = 1, lognu.lb = -1, lognu.ub = 4,
                     n.grid = 5, n.s = 5, beta = 0.1,
                     n.post.samples = 0) {
  X <- .as_ptt_matrix(X, "X")
  Y <- .as_ptt_matrix(Y, "Y")
  if (nrow(X) != nrow(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  p.X <- ncol(X)
  p.Y <- ncol(Y)

  if (xor(is.null(Xpred), is.null(Ypred))) {
    stop("Supply both `Xpred` and `Ypred`, or leave both as `NULL`.", call. = FALSE)
  }
  if (is.null(Xpred)) {
    Xpred <- .empty_prediction(p.X)
    Ypred <- .empty_prediction(p.Y)
  } else {
    Xpred <- .as_ptt_matrix(Xpred, "Xpred", p.X, allow_empty = TRUE)
    Ypred <- .as_ptt_matrix(Ypred, "Ypred", p.Y, allow_empty = TRUE)
    if (nrow(Xpred) != nrow(Ypred)) {
      stop("`Xpred` and `Ypred` must have the same number of rows.", call. = FALSE)
    }
  }

  OmegaX.type <- match.arg(OmegaX.type, c("unit", "standardized"))
  OmegaY.type <- match.arg(OmegaY.type, c("unit", "standardized"))
  max.resX <- .integer_arg(max.resX, "max.resX", 1L, 14L)
  max.resY <- .integer_arg(max.resY, "max.resY", 1L, 14L)
  .validate_tree_capacity(p.X, max.resX, "X")
  .validate_tree_capacity(p.Y, max.resY, "Y")
  rho0.X <- .numeric_arg(rho0.X, "rho0.X", 0, 1)
  rho0.Y <- .numeric_arg(rho0.Y, "rho0.Y", 0, 1)
  rho0.mode.X <- .integer_arg(rho0.mode.X, "rho0.mode.X", 0L, 2L)
  rho0.mode.Y <- .integer_arg(rho0.mode.Y, "rho0.mode.Y", 0L, 2L)
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

  OmegaX <- .make_omega(X, Xpred, OmegaX.type, "X")
  OmegaY <- .make_omega(Y, Ypred, OmegaY.type, "Y")
  ans <- fitCondPTTcpp(
    X, Y, Xpred, Ypred, OmegaX, OmegaY, max.resX, max.resY,
    rho0.X, rho0.mode.X, rho0.Y, rho0.mode.Y, tran.mode,
    lognu.lb, lognu.ub, n.grid, n.s, beta, n.post.samples
  )

  .format_conditional_fit(ans, p.X, n.post.samples, OmegaX, OmegaY)
}

#' Fit a Conditional Optional Pólya Tree
#'
#' A convenience wrapper around [cond.apt()] that uses an optional Pólya tree
#' for the response distribution in each predictor block.
#'
#' @inheritParams cond.apt
#' @returns A list with the same structure as [cond.apt()].
#' @seealso [cond.apt()], [opt()]
#' @export
#'
#' @examples
#' set.seed(12345)
#' x <- runif(60)
#' y <- ifelse(x < 0.5, rbeta(60, 3, 8), rbeta(60, 8, 3))
#' fit <- cond.opt(
#'   x, y, Xpred = c(0.25, 0.75), Ypred = c(0.3, 0.7),
#'   max.resX = 3, max.resY = 4
#' )
#' fit$predictive_densities
cond.opt <- function(X, Y, Xpred = NULL, Ypred = NULL,
                     OmegaX.type = "unit", OmegaY.type = "unit",
                     max.resX = 7, max.resY = 7,
                     rho0.X = 0.5, rho0.mode.X = 0,
                     rho0.Y = 0.5, rho0.mode.Y = 0,
                     n.post.samples = 0) {
  cond.apt(
    X = X, Y = Y, Xpred = Xpred, Ypred = Ypred,
    OmegaX.type = OmegaX.type, OmegaY.type = OmegaY.type,
    max.resX = max.resX, max.resY = max.resY,
    rho0.X = rho0.X, rho0.mode.X = rho0.mode.X,
    rho0.Y = rho0.Y, rho0.mode.Y = rho0.mode.Y,
    tran.mode = 1, lognu.lb = 0, lognu.ub = 0,
    n.grid = 1, n.s = 1, beta = 0,
    n.post.samples = n.post.samples
  )
}
