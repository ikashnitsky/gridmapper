#' Approximate adjacency graph of features
#'
#' Treats two features as adjacent when their bounding boxes touch or
#' overlap, allowing for a small buffer. Edge weights approximate the length
#' of the shared boundary by the larger overlap of the two boxes along the x
#' or y axis. This is fast and needs no topology, but can link features whose
#' boxes overlap without the features touching.
#'
#' @param x An [sf][sf::st_sf] or [sfc][sf::st_sfc] object.
#' @param buffer_fraction Tolerance, as a fraction of the mean feature width.
#'
#' @returns A list of class `gridmapper_adjacency` with elements
#'   * `nodes`: a data frame with an `id` per feature (the `id` column of
#'     `x` if present, otherwise the row number);
#'   * `edges`: a data frame with 1-based `from` and `to` node indices and a
#'     `weight`;
#'   * `bboxes`: a matrix of feature bounding boxes.
#' @export
#' @examples
#' london <- sf::st_read(
#'   system.file("extdata", "london.geojson", package = "gridmapper"),
#'   quiet = TRUE
#' )
#' graph <- adjacency_graph(london)
#' head(graph$edges)
adjacency_graph <- function(x, buffer_fraction = 0.02) {
  check_number(buffer_fraction, min = 0)
  geoms <- as_geometry_list(x)
  ids <- if (inherits(x, "sf") && "id" %in% names(x)) {
    as.character(x$id)
  } else {
    as.character(seq_along(geoms))
  }
  bboxes <- geoms_bbox(geoms)
  edges <- bbox_adjacency(bboxes, buffer_fraction)
  structure(
    list(
      nodes = data.frame(id = ids, index = seq_along(geoms)),
      edges = edges,
      bboxes = bboxes
    ),
    class = "gridmapper_adjacency"
  )
}

#' @export
print.gridmapper_adjacency <- function(x, ...) {
  cli::cli_text(
    "<gridmapper_adjacency> {nrow(x$nodes)} node{?s}, {nrow(x$edges)} edge{?s}"
  )
  invisible(x)
}

bbox_adjacency <- function(bboxes, buffer_fraction = 0.02) {
  n <- nrow(bboxes)
  avg_width <- if (n > 0) sum(bboxes[, "width"]) / n else 0
  buffer <- max(1e-9, avg_width * buffer_fraction)

  min_x <- bboxes[, "min_x"]
  max_x <- bboxes[, "max_x"]
  min_y <- bboxes[, "min_y"]
  max_y <- bboxes[, "max_y"]
  separated <- outer(max_x + buffer, min_x - buffer, "<") |
    outer(min_x - buffer, max_x + buffer, ">") |
    outer(max_y + buffer, min_y - buffer, "<") |
    outer(min_y - buffer, max_y + buffer, ">")

  pairs <- which(!separated & upper.tri(separated), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  i <- pairs[, 1]
  j <- pairs[, 2]
  overlap_x <- pmax(0, pmin(max_x[i], max_x[j]) - pmax(min_x[i], min_x[j]))
  overlap_y <- pmax(0, pmin(max_y[i], max_y[j]) - pmax(min_y[i], min_y[j]))
  shared <- pmax(overlap_x, overlap_y)
  shared[shared <= 0] <- buffer

  data.frame(from = unname(i), to = unname(j), weight = unname(shared))
}

#' Force-directed embedding of an adjacency graph
#'
#' Lays out the nodes of a graph in two dimensions with a
#' Fruchterman-Reingold style simulation: all nodes repel each other and
#' edges pull their endpoints together. The result can serve as
#' `embedding_targets` in [allocate_grid()]. Starting positions are random,
#' so call [set.seed()] for reproducible results.
#'
#' @param graph An adjacency graph from [adjacency_graph()].
#' @param iterations Number of simulation steps.
#' @param width,height Size of the layout area.
#' @param jitter Minimum step size per iteration.
#'
#' @returns A numeric matrix with columns `x` and `y` and one row per node.
#' @export
#' @examples
#' london <- sf::st_read(
#'   system.file("extdata", "london.geojson", package = "gridmapper"),
#'   quiet = TRUE
#' )
#' set.seed(1)
#' head(embed_graph(adjacency_graph(london)))
embed_graph <- function(
  graph,
  iterations = 400,
  width = 1,
  height = 1,
  jitter = 0.01
) {
  if (!inherits(graph, "gridmapper_adjacency")) {
    cli::cli_abort("{.arg graph} must be created by {.fn adjacency_graph}.")
  }
  check_count(iterations)
  check_number(width, min = 0)
  check_number(height, min = 0)
  check_number(jitter, min = 0)

  n <- nrow(graph$nodes)
  x <- stats::runif(n)
  y <- stats::runif(n)

  e <- graph$edges
  w <- e$weight
  w[is.na(w) | w == 0] <- 1
  src <- c(e$from, e$to)
  dst <- c(e$to, e$from)
  ww <- c(w, w)
  k <- sqrt(1 / max(1, n))

  for (iter in seq_len(iterations) - 1) {
    dx <- outer(x, x, "-")
    dy <- outer(y, y, "-")
    dist <- sqrt(dx^2 + dy^2 + 1e-9)
    f <- k^2 / dist
    disp_x <- rowSums(dx / dist * f)
    disp_y <- rowSums(dy / dist * f)

    if (length(src)) {
      adx <- x[src] - x[dst]
      ady <- y[src] - y[dst]
      ad <- sqrt(adx^2 + ady^2 + 1e-9)
      af <- ad^2 / k / pmax(1e-6, ww)
      disp_x <- disp_x - tapply_sum(adx / ad * af, src, n)
      disp_y <- disp_y - tapply_sum(ady / ad * af, src, n)
    }

    temp <- (1 - iter / iterations) * 0.1 + jitter
    len <- sqrt(disp_x^2 + disp_y^2)
    step <- ifelse(len > 1e-9, pmin(len, temp) / len, 0)
    x <- pmin(1, pmax(0, x + disp_x * step))
    y <- pmin(1, pmax(0, y + disp_y * step))
  }

  cbind(x = x * width, y = y * height)
}

tapply_sum <- function(values, index, n) {
  out <- numeric(n)
  sums <- rowsum(values, index)
  out[as.integer(rownames(sums))] <- sums[, 1]
  out
}
