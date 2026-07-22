#' Plot a Two-Dimensional Partition
#'
#' Draws rectangular partition cells with optional density-based fill colors.
#' This helper uses grid graphics and is intended for two-dimensional HMAP or
#' posterior partitions after their integer bounds have been converted to the
#' desired plotting coordinates.
#'
#' @param xy.part A matrix or data frame with columns `xmin`, `xmax`, `ymin`,
#'   `ymax`, and `den`.
#' @param border `FALSE` uses the fill colors as borders, `TRUE` uses black, and
#'   a single color string uses that border color.
#' @param xlim,ylim Horizontal and vertical plotting limits.
#' @param zlim Density limits used to map `den` to colors.
#' @param main Plot title.
#' @param nlevels Number of colors in the palette; at least two.
#' @param plot.den Whether to color cells by `den`.
#' @param color.fun A function that returns `nlevels` colors.
#' @param plot.scale Whether to draw a density color scale.
#' @param color Optional fixed fill color. The default `NA` uses `color.fun`.
#' @param newpage Whether to start a new grid graphics page.
#' @param x A partition object for the legacy `plot.part()` S3 method.
#' @param ... Additional arguments passed to [plot_partition()].
#'
#' @returns Invisibly, the clipped partition data with an added `fill` column.
#' @export
#'
#' @examples
#' cells <- data.frame(
#'   xmin = c(0, 0.5), xmax = c(0.5, 1),
#'   ymin = c(0, 0), ymax = c(1, 1), den = c(0.6, 1.4)
#' )
#' plot_partition(cells, plot.scale = FALSE)
plot_partition <- function(xy.part, border = FALSE, xlim = c(0, 1), ylim = c(0, 1),
                           zlim = grDevices::extendrange(xy.part[, "den"]), main = "",
                           nlevels = 100, plot.den = TRUE,
                           color.fun = grDevices::terrain.colors, plot.scale = TRUE,
                           color = NA, newpage = TRUE) {
  xy.part <- as.data.frame(xy.part)
  required <- c("xmin", "xmax", "ymin", "ymax", "den")
  if (!all(required %in% names(xy.part))) {
    stop("`xy.part` must contain columns xmin, xmax, ymin, ymax, and den.", call. = FALSE)
  }
  if (!all(vapply(xy.part[required], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(xy.part[required])))) {
    stop("The required `xy.part` columns must contain finite numeric values.", call. = FALSE)
  }
  if (nrow(xy.part) < 1L) {
    stop("`xy.part` must contain at least one partition cell.", call. = FALSE)
  }

  nlevels <- .integer_arg(nlevels, "nlevels", 2L)
  xlim <- as.double(xlim)
  ylim <- as.double(ylim)
  zlim <- as.double(zlim)
  if (length(xlim) != 2L || any(!is.finite(xlim)) || xlim[1L] >= xlim[2L] ||
      length(ylim) != 2L || any(!is.finite(ylim)) || ylim[1L] >= ylim[2L]) {
    stop("`xlim` and `ylim` must be finite, increasing two-element vectors.", call. = FALSE)
  }
  if (length(zlim) != 2L || any(!is.finite(zlim)) || zlim[1L] > zlim[2L]) {
    stop("`zlim` must be a finite, increasing two-element vector.", call. = FALSE)
  }
  zlim[1L] <- max(0, zlim[1L])
  if (zlim[1L] == zlim[2L]) {
    delta <- max(abs(zlim[1L]), 1) * sqrt(.Machine$double.eps)
    zlim <- c(max(0, zlim[1L] - delta), zlim[2L] + delta)
  }

  palette <- color.fun(nlevels)
  if (length(palette) != nlevels) {
    stop("`color.fun` must return exactly `nlevels` colors.", call. = FALSE)
  }
  density <- pmin(pmax(xy.part$den, zlim[1L]), zlim[2L])
  color.index <- 1L + floor(
    (density - zlim[1L]) / diff(zlim) * (nlevels - 1L)
  )

  if (!all(is.na(color))) {
    fill.colors <- rep(color, length.out = nrow(xy.part))
  } else if (isTRUE(plot.den)) {
    fill.colors <- palette[color.index]
  } else {
    fill.colors <- rep("transparent", nrow(xy.part))
  }

  if (isTRUE(border)) {
    border.colors <- "black"
  } else if (identical(border, FALSE)) {
    border.colors <- fill.colors
  } else if (is.character(border) && length(border) == 1L && !is.na(border)) {
    border.colors <- border
  } else {
    stop("`border` must be TRUE, FALSE, or a single color string.", call. = FALSE)
  }

  xy.part$xmin <- pmin(pmax(xy.part$xmin, xlim[1L]), xlim[2L])
  xy.part$xmax <- pmin(pmax(xy.part$xmax, xlim[1L]), xlim[2L])
  xy.part$ymin <- pmin(pmax(xy.part$ymin, ylim[1L]), ylim[2L])
  xy.part$ymax <- pmin(pmax(xy.part$ymax, ylim[1L]), ylim[2L])
  keep <- xy.part$xmax > xy.part$xmin & xy.part$ymax > xy.part$ymin
  xy.part <- xy.part[keep, , drop = FALSE]
  fill.colors <- fill.colors[keep]
  border.colors <- rep(border.colors, length.out = length(keep))[keep]

  if (isTRUE(newpage)) {
    grid::grid.newpage()
  }
  if (isTRUE(plot.scale)) {
    vp.plot <- grid::viewport(
      x = grid::unit(4, "lines"),
      y = grid::unit(0.2, "npc") + grid::unit(2, "lines"),
      just = c("left", "bottom"),
      width = grid::unit(1, "npc") - grid::unit(6, "lines"),
      height = grid::unit(0.70, "npc"),
      xscale = xlim, yscale = ylim, name = "plotRegion"
    )
    vp.scale <- grid::viewport(
      x = grid::unit(4, "lines"), y = grid::unit(4, "lines"),
      just = c("left", "bottom"),
      width = grid::unit(1, "npc") - grid::unit(6, "lines"),
      height = grid::unit(0.2, "npc"), name = "scaleRegion"
    )
    grid::pushViewport(vp.plot)
  } else {
    grid::pushViewport(grid::plotViewport(c(5, 4, 2, 2)))
    grid::pushViewport(grid::viewport(
      xscale = xlim, yscale = ylim, name = "plotRegion"
    ))
  }

  grid::grid.rect(
    x = grid::unit((xy.part$xmin + xy.part$xmax) / 2, "native"),
    y = grid::unit((xy.part$ymin + xy.part$ymax) / 2, "native"),
    width = grid::unit(xy.part$xmax - xy.part$xmin, "native"),
    height = grid::unit(xy.part$ymax - xy.part$ymin, "native"),
    gp = grid::gpar(col = border.colors, fill = fill.colors)
  )
  grid::popViewport(if (isTRUE(plot.scale)) 1L else 2L)

  if (isTRUE(plot.scale)) {
    grid::pushViewport(vp.scale)
    width <- 1 / nlevels
    grid::grid.rect(
      x = grid::unit((seq_len(nlevels) - 0.5) * width, "npc"),
      y = grid::unit(0, "npc"), just = c("centre", "bottom"),
      height = grid::unit(0.2, "npc"), width = grid::unit(width, "npc"),
      gp = grid::gpar(col = NA, fill = palette)
    )
    grid::grid.xaxis(
      at = c(0, 0.5, 1),
      label = signif(c(zlim[1L], mean(zlim), zlim[2L]), 3)
    )
    grid::grid.text("Density scale", y = grid::unit(0.1, "npc"))
    grid::popViewport()
  }

  grid::grid.text(main, y = grid::unit(1, "npc") - grid::unit(1, "lines"))
  xy.part$fill <- fill.colors
  invisible(xy.part)
}

#' @rdname plot_partition
#' @export
#' @rawNamespace export(plot.part)
plot.part <- function(x, ...) {
  dots <- list(...)
  if (missing(x)) {
    if (!"xy.part" %in% names(dots)) {
      stop("A partition must be supplied as `x` or `xy.part`.", call. = FALSE)
    }
    x <- dots$xy.part
    dots$xy.part <- NULL
  }
  do.call(plot_partition, c(list(xy.part = x), dots))
}
