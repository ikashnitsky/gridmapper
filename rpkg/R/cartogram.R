#' Dorling and Demers cartograms
#'
#' Replaces each feature by a circle, hexagon or square whose area is
#' proportional to `weight`, centred on the feature's [polygon_centroid()],
#' then runs a force simulation that removes overlaps while anchoring shapes
#' to their original positions.
#'
#' `demers_cartogram()` is `dorling_cartogram()` with squares, adjacency
#' attraction and a strict final overlap resolver.
#'
#' @param x An [sf][sf::st_sf] object.
#' @param weight Name of a numeric column of `x` holding the weights.
#'   Features with missing, non-positive or (if given) below-`min_weight`
#'   weights are dropped.
#' @param shape Shape of the symbols.
#' @param k Size of the largest symbol, as the percentage of the bounding
#'   box area it covers.
#' @param iterations Number of force simulation iterations.
#' @param repulsion Initial strength of the overlap repulsion.
#' @param anchor Strength of the pull back to the original positions.
#' @param preserve_adjacency If `TRUE`, neighbouring features (see
#'   [adjacency_graph()]) attract each other.
#' @param adjacency_strength Strength of the adjacency attraction.
#' @param cooling Factor by which repulsion decays each iteration.
#' @param min_weight Optional minimum weight.
#' @param seed Optional seed for the simulation's own random number
#'   generator. With a seed, results are reproducible and identical to the
#'   JavaScript library. Without one, the seed is drawn from R's random
#'   number generator, so [set.seed()] also gives reproducible results.
#' @param refine_overlaps If `TRUE` and overlaps remain, run a second,
#'   unanchored simulation.
#' @param refinement_iterations Number of iterations of the refinement.
#' @param strict_no_overlap If `TRUE`, squares and rectangles that still
#'   overlap are separated by a deterministic resolver.
#' @param ... Arguments passed on to `dorling_cartogram()`.
#'
#' @returns An [sf][sf::st_sf] object with the retained features of `x`, their
#'   geometry replaced by the cartogram shapes, and columns `cartogram_x`,
#'   `cartogram_y` (shape centres), `cartogram_radius`, and
#'   `cartogram_original_x`, `cartogram_original_y` (the original centroids).
#'   Diagnostics are stored in the `"cartogram"` attribute, including
#'   `has_overlaps`.
#' @export
#' @examples
#' london <- sf::st_read(
#'   system.file("extdata", "london.geojson", package = "gridmapper"),
#'   quiet = TRUE
#' )
#' london$area <- as.numeric(sf::st_area(sf::st_transform(london, 27700)))
#' dorl <- dorling_cartogram(london, "area", k = 3, seed = 1, iterations = 500)
#' plot(sf::st_geometry(dorl))
#' attr(dorl, "cartogram")$has_overlaps
#'
#' dem <- demers_cartogram(london, "area", k = 3, seed = 1, iterations = 500)
#' plot(sf::st_geometry(dem))
dorling_cartogram <- function(
  x,
  weight,
  shape = c("circle", "hexagon", "square", "rectangle"),
  k = 10,
  iterations = 2000,
  repulsion = 1,
  anchor = 0.01,
  preserve_adjacency = FALSE,
  adjacency_strength = 0.05,
  cooling = 0.998,
  min_weight = NULL,
  seed = NULL,
  refine_overlaps = TRUE,
  refinement_iterations = 1200,
  strict_no_overlap = FALSE
) {
  if (!inherits(x, "sf")) {
    cli::cli_abort(
      "{.arg x} must be an {.cls sf} object, not {.obj_type_friendly {x}}."
    )
  }
  if (!rlang::is_string(weight) || !weight %in% names(x)) {
    cli::cli_abort("{.arg weight} must be the name of a column of {.arg x}.")
  }
  shape <- rlang::arg_match(shape)
  check_number(k, min = 0)
  check_count(iterations)
  check_number(repulsion)
  check_number(anchor)
  check_flag(preserve_adjacency)
  check_number(adjacency_strength)
  check_number(cooling, min = 0, max = 1)
  check_flag(refine_overlaps)
  check_count(refinement_iterations)
  check_flag(strict_no_overlap)
  if (!is.null(seed)) {
    check_number(seed)
  }

  w <- x[[weight]]
  if (!is.numeric(w)) {
    cli::cli_abort("Column {.field {weight}} must be numeric.")
  }
  valid <- !is.na(w) & is.finite(w) & w > 0
  if (!is.null(min_weight)) {
    valid <- valid & w >= min_weight
  }
  if (!any(valid)) {
    cli::cli_abort(
      "No features with valid positive weights found in column {.field {weight}}."
    )
  }
  x <- x[valid, ]
  w <- w[valid]

  cent <- polygon_centroid(x)
  if (anyNA(cent)) {
    cli::cli_abort(
      "Can't compute the centroid of feature{?s} {which(is.na(cent[, 1]))}."
    )
  }
  b <- calc_bounds(cent[, 1], cent[, 2])
  scale <- or_one(max(b$width, b$height))

  # Simulate in unit space
  nx <- (cent[, 1] - b$min_x) / scale
  ny <- (cent[, 2] - b$min_y) / scale
  scale_factor <- (b$width / scale) * (b$height / scale) * k / (100 * max(w))
  radius <- sqrt(abs(w) * scale_factor / pi)
  rect <- shape %in% c("square", "rectangle")
  side <- if (rect) radius * sqrt(pi) else rep(0, length(radius))

  edges <- if (preserve_adjacency) {
    adjacency_graph(x)$edges
  } else {
    data.frame(from = integer(), to = integer())
  }
  sim_seed <- function() seed %||% floor(stats::runif(1, 1, 4294967295))

  run <- function(
    px,
    py,
    iters,
    rep_strength,
    anchor_strength,
    adj,
    cool,
    min_rep,
    s,
    jitter
  ) {
    force_simulation_cpp(
      px,
      py,
      nx,
      ny,
      radius,
      side,
      side,
      rect,
      as.integer(edges$from - 1L),
      as.integer(edges$to - 1L),
      as.integer(iters),
      rep_strength,
      anchor_strength,
      adj,
      cool,
      min_rep,
      s,
      jitter
    )
  }

  s <- sim_seed()
  pos <- run(
    nx,
    ny,
    iterations,
    repulsion,
    anchor,
    adjacency_strength,
    cooling,
    0.5,
    s,
    if (is.null(seed)) 0 else 1e-6
  )

  refined <- FALSE
  if (
    refine_overlaps && has_overlaps_cpp(pos$x, pos$y, radius, side, side, rect)
  ) {
    refined <- TRUE
    pos <- run(
      pos$x,
      pos$y,
      refinement_iterations,
      max(1.2, repulsion),
      0,
      if (preserve_adjacency) adjacency_strength * 0.5 else 0,
      0.999,
      0.9,
      sim_seed(),
      0
    )
  }

  strict <- NULL
  if (
    strict_no_overlap &&
      rect &&
      has_overlaps_cpp(pos$x, pos$y, radius, side, side, rect)
  ) {
    strict <- resolve_rect_overlaps_cpp(
      pos$x,
      pos$y,
      nx,
      ny,
      side,
      side,
      as.integer(max(2000, refinement_iterations)),
      1e-4,
      1e-4,
      0.001
    )
    pos$x <- strict$x
    pos$y <- strict$y
    strict$x <- NULL
    strict$y <- NULL
  }

  # Back to map coordinates
  cx <- pos$x * scale + b$min_x
  cy <- pos$y * scale + b$min_y
  r <- radius * scale
  s_side <- side * scale
  still_overlaps <- has_overlaps_cpp(cx, cy, r, s_side, s_side, rect)

  polys <- lapply(seq_along(cx), function(i) {
    ring <- switch(
      shape,
      circle = circle_ring(cx[i], cy[i], r[i]),
      hexagon = hexagon_ring(cx[i], cy[i], r[i]),
      rectangle_ring(cx[i], cy[i], s_side[i], s_side[i])
    )
    sf::st_polygon(list(ring))
  })

  out <- x
  out$cartogram_x <- cx
  out$cartogram_y <- cy
  out$cartogram_radius <- r
  out$cartogram_original_x <- cent[, 1]
  out$cartogram_original_y <- cent[, 2]
  sf::st_geometry(out) <- sf::st_sfc(polys, crs = sf::st_crs(x))
  attr(out, "cartogram") <- list(
    weight = weight,
    shape = shape,
    k = k,
    iterations = iterations,
    seed = seed,
    refinement_applied = refined,
    refinement_iterations = refinement_iterations,
    strict_no_overlap = strict_no_overlap,
    strict_resolver = strict,
    feature_count = length(cx),
    has_overlaps = still_overlaps,
    bounds = b,
    scale_factor = scale_factor
  )
  out
}

#' @rdname dorling_cartogram
#' @export
demers_cartogram <- function(x, weight, ...) {
  defaults <- list(
    shape = "square",
    preserve_adjacency = TRUE,
    adjacency_strength = 0.08,
    repulsion = 1.25,
    anchor = 0.006,
    refine_overlaps = TRUE,
    refinement_iterations = 1800,
    strict_no_overlap = TRUE
  )
  args <- utils::modifyList(defaults, list(...))
  do.call(dorling_cartogram, c(list(x = x, weight = weight), args))
}
