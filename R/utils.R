# Internal validation and formatting helpers ---------------------------------

.as_ptt_matrix <- function(x, name, ncol_expected = NULL, allow_empty = FALSE) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("`%s` must contain only numeric columns.", name), call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (is.numeric(x) && is.null(dim(x))) {
    if (is.null(ncol_expected) || identical(ncol_expected, 1L)) {
      x <- matrix(x, ncol = 1L)
    } else if (length(x) == ncol_expected) {
      x <- matrix(x, nrow = 1L)
    } else {
      stop(
        sprintf("`%s` must have %d columns.", name, ncol_expected),
        call. = FALSE
      )
    }
  } else if (!is.matrix(x)) {
    stop(sprintf("`%s` must be a numeric vector, matrix, or data frame.", name), call. = FALSE)
  }

  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"

  if (ncol(x) < 1L || (!allow_empty && nrow(x) < 1L)) {
    stop(sprintf("`%s` must contain at least one observation and one column.", name), call. = FALSE)
  }
  if (!is.null(ncol_expected) && ncol(x) != ncol_expected) {
    stop(sprintf("`%s` must have %d columns.", name, ncol_expected), call. = FALSE)
  }
  if (length(x) && any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }

  x
}

.empty_prediction <- function(p) {
  matrix(numeric(), nrow = 0L, ncol = p)
}

.integer_arg <- function(x, name, lower, upper = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) ||
      x < lower || x > upper) {
    stop(
      sprintf("`%s` must be a single integer between %d and %d.", name, lower, upper),
      call. = FALSE
    )
  }
  as.integer(x)
}

.validate_tree_capacity <- function(dimensions, resolution, name) {
  max_dimensions <- 65535L - resolution + 1L
  if (dimensions > max_dimensions) {
    stop(
      sprintf(
        "`%s` has %d dimensions, but resolution %d supports at most %d.",
        name, dimensions, resolution, max_dimensions
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.numeric_arg <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < lower || x > upper) {
    stop(
      sprintf("`%s` must be a finite number between %s and %s.", name, lower, upper),
      call. = FALSE
    )
  }
  as.double(x)
}

.make_omega <- function(x, xpred, type, name) {
  values <- if (nrow(xpred)) rbind(x, xpred) else x
  p <- ncol(x)

  if (identical(type, "unit")) {
    if (any(values < 0 | values > 1)) {
      stop(
        sprintf("`%s` and its prediction points must lie in [0, 1] when the sample space is `unit`.", name),
        call. = FALSE
      )
    }
    return(matrix(rep(c(0, 1), p), nrow = p, ncol = 2L, byrow = TRUE))
  }

  bounds <- apply(values, 2L, range)
  if (p == 1L) {
    bounds <- matrix(bounds, nrow = 2L)
  }
  lower <- bounds[1L, ]
  upper <- bounds[2L, ]

  # Pad both ends. Unlike multiplying the upper endpoint, this also works for
  # negative and constant-valued data.
  scale <- pmax(1, abs(lower), abs(upper))
  padding <- pmax((upper - lower) * 1e-7, sqrt(.Machine$double.eps) * scale)
  cbind(lower - padding, upper + padding)
}

.partition_names <- function(p, suffix = character()) {
  bounds <- as.vector(rbind(paste0("X", seq_len(p), ".l"), paste0("X", seq_len(p), ".u")))
  c(bounds, suffix)
}

.format_marginal_fit <- function(ans, p, n.s, n.post.samples, omega) {
  hmap <- matrix(
    unlist(ans$part_points_hmap, use.names = FALSE),
    ncol = 2L * p + 2L,
    byrow = TRUE
  )
  colnames(hmap) <- .partition_names(p, c("level", "state"))
  hmap[hmap[, "state"] == n.s + 1L, "state"] <- Inf
  ans$part_points_hmap <- hmap

  if (n.post.samples > 0L) {
    for (i in seq_len(n.post.samples)) {
      part <- matrix(
        unlist(ans$part_points_post_samples[[i]], use.names = FALSE),
        ncol = 2L * p + 1L,
        byrow = TRUE
      )
      colnames(part) <- .partition_names(p, "level")
      parameters <- matrix(
        unlist(ans$nu_and_prob_post_samples[[i]], use.names = FALSE),
        ncol = 2L,
        byrow = TRUE,
        dimnames = list(NULL, c("nu", "logp"))
      )
      ans$part_points_post_samples[[i]] <- cbind(part, parameters)
    }
  } else {
    ans$part_points_post_samples <- NULL
  }

  ans$nu_and_prob_post_samples <- NULL
  ans$Omega <- omega
  ans
}

.format_conditional_fit <- function(ans, p, n.post.samples, omega.x, omega.y) {
  hmap <- matrix(
    unlist(ans$part_points_hmap, use.names = FALSE),
    ncol = 2L * p + 2L,
    byrow = TRUE
  )
  colnames(hmap) <- .partition_names(p, c("level", "state"))
  hmap[hmap[, "state"] == 2L, "state"] <- Inf

  rho <- matrix(
    unlist(ans$rhos_hmap, use.names = FALSE),
    ncol = 1L,
    dimnames = list(NULL, "rho")
  )
  ans$part_points_hmap <- cbind(hmap, rho)
  ans$rhos_hmap <- NULL

  if (n.post.samples > 0L) {
    for (i in seq_len(n.post.samples)) {
      part <- matrix(
        unlist(ans$part_points_post_samples[[i]], use.names = FALSE),
        ncol = 2L * p + 1L,
        byrow = TRUE
      )
      colnames(part) <- .partition_names(p, "level")
      parameters <- matrix(
        unlist(ans$nu_and_prob_post_samples[[i]], use.names = FALSE),
        ncol = 2L,
        byrow = TRUE,
        dimnames = list(NULL, c("nu", "logp"))
      )
      ans$part_points_post_samples[[i]] <- cbind(part, parameters)
    }
  } else {
    ans$part_points_post_samples <- NULL
  }

  ans$nu_and_prob_post_samples <- NULL
  ans$OmegaX <- omega.x
  ans$OmegaY <- omega.y
  ans
}
