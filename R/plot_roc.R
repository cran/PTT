#' Plot an Empirical ROC Curve
#'
#' Uses empirical null-score quantiles as thresholds. The curve treats smaller
#' scores as stronger evidence for the positive class.
#'
#' @param p Numeric scores under the alternative or positive class.
#' @param p0 Numeric scores under the null or negative class.
#' @param fpr False-positive-rate values in `[0, 1]`.
#' @param col,lty,lwd Graphical parameters passed to [graphics::plot()].
#' @param xlim,ylim Axis limits.
#' @param main Plot title.
#' @param x Scores for the legacy `plot.roc()` S3 method.
#' @param ... Additional arguments passed to [plot_roc_curve()].
#'
#' @returns Invisibly, a data frame containing `fpr`, `threshold`, and `tpr`.
#' @export
#'
#' @examples
#' set.seed(12345)
#' null <- rnorm(100)
#' alternative <- rnorm(100, mean = -1)
#' roc <- plot_roc_curve(alternative, null)
#' head(roc)
plot_roc_curve <- function(p, p0, fpr = seq(0, 1, by = 0.01), col = "black",
                           xlim = c(0, 1), ylim = c(0, 1), main = "",
                           lty = 1, lwd = 2) {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) ||
      !is.numeric(p0) || !length(p0) || any(!is.finite(p0))) {
    stop("`p` and `p0` must be non-empty finite numeric vectors.", call. = FALSE)
  }
  if (!is.numeric(fpr) || !length(fpr) || any(!is.finite(fpr)) ||
      any(fpr < 0 | fpr > 1)) {
    stop("`fpr` must contain finite probabilities in [0, 1].", call. = FALSE)
  }

  threshold <- stats::quantile(p0, probs = fpr, names = FALSE)
  tpr <- stats::ecdf(p)(threshold)
  graphics::plot(
    fpr, tpr, type = "l", lty = lty,
    xlab = "False rejection rate", ylab = "True rejection rate",
    col = col, xlim = xlim, ylim = ylim, main = main, lwd = lwd
  )
  invisible(data.frame(fpr = fpr, threshold = threshold, tpr = tpr))
}

#' @rdname plot_roc_curve
#' @export
#' @rawNamespace export(plot.roc)
plot.roc <- function(x, ...) {
  dots <- list(...)
  if (missing(x)) {
    if (!"p" %in% names(dots)) {
      stop("Scores must be supplied as `x` or `p`.", call. = FALSE)
    }
    x <- dots$p
    dots$p <- NULL
  }
  do.call(plot_roc_curve, c(list(p = x), dots))
}
