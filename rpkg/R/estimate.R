#' Estimate grid parameters from geographic features
#'
#' Recommends grid dimensions, compactness and whether to use PCA rotation
#' for [allocate_grid()], using heuristics on the shape of the features:
#' elongation and diagonality of the point cloud, how much of the map the
#' features cover, holes, adjacency and connected components. See
#' `vignette("algorithms")` in the original JavaScript library for details.
#'
#' All geometry is computed on the raw coordinates, as in the JavaScript
#' library.
#'
#' @param x An [sf][sf::st_sf] object with `POINT`, `POLYGON` or
#'   `MULTIPOLYGON` geometries.
#' @param dimension_search If `TRUE`, refine the dimensions by scoring
#'   nearby candidate grids.
#' @param embedding If `TRUE`, also compute an [adjacency_graph()] and a
#'   force-directed [embed_graph()] layout scaled to the map bounds, for use
#'   as `embedding_targets` in [allocate_grid()].
#' @param embedding_options A list of arguments passed to [embed_graph()].
#'
#' @returns A list with elements `rows`, `cols`, `compactness` and
#'   `rotate_pca`, plus `adjacency_graph` and `embedding_targets` when
#'   `embedding = TRUE`.
#' @seealso [estimate_parameters_from_data()] for a simpler estimate from
#'   point coordinates.
#' @export
#' @examples
#' london <- sf::st_read(
#'   system.file("extdata", "london.geojson", package = "gridmapper"),
#'   quiet = TRUE
#' )
#' params <- estimate_parameters(london)
#' params
#' gm <- allocate_grid(
#'   london,
#'   rows = params$rows,
#'   cols = params$cols,
#'   compactness = params$compactness,
#'   rotate_pca = params$rotate_pca
#' )
estimate_parameters <- function(
  x,
  dimension_search = TRUE,
  embedding = FALSE,
  embedding_options = list()
) {
  check_flag(dimension_search)
  check_flag(embedding)
  geoms <- as_geometry_list(x)
  n <- length(geoms)
  if (n == 0) {
    cli::cli_abort("{.arg x} must contain at least one feature.")
  }

  # 1. Centroids and geometry descriptors
  cent <- t(vapply(
    geoms,
    \(g) {
      if (g$type %in% c("Point", "Polygon", "MultiPolygon")) {
        geom_centroid(g)
      } else {
        c(0, 0)
      }
    },
    numeric(2)
  ))
  px <- cent[, 1]
  py <- cent[, 2]
  fb <- geoms_bbox(geoms)
  stats <- t(vapply(geoms, geom_stats, numeric(4)))

  global <- bounds_of(
    fb[, "min_x"],
    fb[, "max_x"],
    fb[, "min_y"],
    fb[, "max_y"]
  )

  # 2. PCA rotation decision
  angle <- pca_angle(px, py)
  center <- c(mean(px), mean(py))
  rot <- pca_rotate(px, py, -angle, center)
  rot_pts <- calc_bounds(rot[, 1], rot[, 2])
  pts <- calc_bounds(px, py)

  corners_x <- c(fb[, "min_x"], fb[, "max_x"], fb[, "max_x"], fb[, "min_x"])
  corners_y <- c(fb[, "min_y"], fb[, "min_y"], fb[, "max_y"], fb[, "max_y"])
  rc <- pca_rotate(corners_x, corners_y, -angle, center)
  rot_global <- calc_bounds(rc[, 1], rc[, 2])

  directional <- directional_indicators(px, py, angle)
  geo <- geo_indicators(px, py, stats, fb)

  aspect <- function(b) max(b$width, b$height) / min(b$width, b$height)
  area_original <- global$width * global$height
  area_rotated <- rot_global$width * rot_global$height
  points_area_gain <- (pts$width * pts$height) /
    max(rot_pts$width * rot_pts$height, 1e-9)
  aspect_original <- aspect(global)
  aspect_rotated <- aspect(rot_global)
  points_aspect_original <- aspect(pts)
  points_aspect_rotated <- aspect(rot_pts)

  rotate <- isTRUE(area_original > area_rotated * 1.1) ||
    isTRUE(aspect_rotated > aspect_original * 1.5) ||
    isTRUE(aspect_original > 1.5 && area_original > area_rotated * 1.05) ||
    isTRUE(points_aspect_rotated > points_aspect_original * 1.5)

  # Keep the diagonal mental map of strongly linear chains (Japan, Chile)
  # unless PCA gives a large packing improvement.
  diagonal_chain <- isTRUE(
    directional$linearity > 0.52 &&
      directional$diagonality > 0.55 &&
      points_aspect_rotated > 2.2
  )
  if (rotate && diagonal_chain && points_area_gain < 1.1) {
    rotate <- FALSE
  }

  # 3. Grid dimensions
  active <- if (rotate) rot_pts else pts
  intrinsic_ratio <- rot_pts$height / max(rot_pts$width, 1e-9)
  raw_ratio <- active$height / max(active$width, 1e-9)
  elong <- function(b) {
    max(b$width, b$height) / max(1e-9, min(b$width, b$height))
  }
  intrinsic_elongation <- elong(rot_pts)
  elongation <- elong(active)
  if (!rotate && diagonal_chain) {
    raw_ratio <- min(raw_ratio, intrinsic_ratio * 1.05)
    elongation <- max(elongation, intrinsic_elongation)
  }

  slack_base <- if (rotate) 0.55 else 0.35
  max_slack <- if (rotate) 1.6 else 1.0
  slack <- min(max_slack, max(0, (elongation - 1.3) * slack_base))

  map_area <- max(global$width * global$height, 1e-9)
  feature_area <- feature_areas(stats, fb)
  coverage_ratio <- min(1, sum(feature_area) / map_area)
  coverage_boost <- max(0, coverage_ratio - 0.35)
  circ <- stats[, "circularity"]
  circularity_median <- js_median(circ[is.finite(circ) & circ > 0])
  hole_ratio <- geo$hole_ratio
  occupancy_ratio <- geo$occupancy_ratio
  component_count <- geo$component_count
  component_penalty <- max(0, component_count - 1) / max(1, sqrt(n))
  sparsity <- geo$sparsity_index

  raw_log_ratio <- log(max(raw_ratio, 1e-9))
  log_compression <- if (rotate) 0.85 else 0.75
  if (!rotate && diagonal_chain) {
    log_compression <- 0.96
  }
  if (elongation > 2.2) {
    log_compression <- log_compression + 0.08
  }
  if (elongation > 3.5) {
    log_compression <- log_compression + 0.05
  }
  if (coverage_ratio > 0.55 && n < 90) {
    log_compression <- log_compression - 0.08
  }
  if (sparsity > 0.45) {
    log_compression <- log_compression + 0.04
  }
  if (hole_ratio > 0.1) {
    log_compression <- log_compression + 0.03
  }
  log_compression <- max(0.65, min(0.98, log_compression))

  ratio <- exp(raw_log_ratio * log_compression)
  min_ratio <- if (elongation > 4) {
    0.18
  } else if (elongation > 2.5) {
    0.24
  } else {
    0.3
  }
  ratio <- max(min_ratio, min(1 / min_ratio, ratio))
  if (n < 50 && coverage_boost > 0.25 && elongation < 1.9) {
    ratio <- ratio * 0.7 + 0.3 * 0.85
  }
  if (!diagonal_chain && !rotate && coverage_ratio > 0.32 && elongation < 1.8) {
    ratio <- ratio * 0.55 + 0.45 * 0.9
  }

  adjacency_boost <- max(0, geo$average_degree - 2.5) * 0.1
  if (elongation > 2.2) {
    adjacency_boost <- adjacency_boost * 0.6
  }
  if (coverage_ratio < 0.3) {
    adjacency_boost <- adjacency_boost * 0.7
  }
  if (circularity_median > 0.45) {
    adjacency_boost <- adjacency_boost * 0.85
  }
  density_boost <- coverage_boost + adjacency_boost
  sparse_boost <- max(0, sparsity - 0.35)
  occupancy_boost <- max(0, 0.32 - occupancy_ratio)

  target_cells <- ceiling(n * (1 + slack + density_boost * 0.45))
  target_cells <- ceiling(
    target_cells *
      (1 +
        sparse_boost * 0.35 +
        occupancy_boost * 0.25 +
        component_penalty * 0.4 +
        hole_ratio * 0.2)
  )
  if (coverage_ratio > 0.45 && elongation < 1.8 && circularity_median > 0.35) {
    target_cells <- max(n, ceiling(target_cells * 0.96))
  }
  if (
    !diagonal_chain &&
      !rotate &&
      directional$diagonality < 0.35 &&
      coverage_ratio > 0.22 &&
      elongation < 1.9
  ) {
    target_cells <- max(n, ceiling(target_cells * 0.92))
  }

  elongation_boost <- max(0, log2(max(elongation, 1)))
  sqrt_n <- sqrt(n)
  min_rows <- max(
    1,
    ceiling(sqrt_n * (1 + 0.10 * elongation_boost + density_boost * 0.18))
  )
  min_cols <- max(
    1,
    ceiling(
      sqrt_n *
        (1 +
          (if (rotate) 0.35 else 0.18) * elongation_boost +
          density_boost * 0.2)
    )
  )
  min_cols <- min(min_cols, max(3, ceiling(sqrt_n * 3)))

  dims <- grid_dimensions(
    n,
    aspect_ratio = ratio,
    target_cells = target_cells,
    min_rows = min_rows,
    min_cols = min_cols
  )
  rows <- dims[["rows"]]
  cols <- dims[["cols"]]

  if (!diagonal_chain && !rotate && coverage_ratio > 0.5 && elongation < 2.5) {
    while (cols - rows > 1 && (rows + 1) * (cols - 1) >= n) {
      rows <- rows + 1L
      cols <- cols - 1L
    }
  }

  if (dimension_search) {
    refined <- refine_dimensions(
      base_rows = rows,
      base_cols = cols,
      n = n,
      target_cells = target_cells,
      target_aspect = ratio,
      min_rows = min_rows,
      min_cols = min_cols,
      diagonal_chain = diagonal_chain,
      coverage_ratio = coverage_ratio,
      sparsity_index = sparsity,
      elongation = elongation
    )
    rows <- refined[["rows"]]
    cols <- refined[["cols"]]
  }

  # 4. Compactness
  compactness <- estimate_compactness_advanced(
    px,
    py,
    feature_area,
    global,
    average_degree = geo$average_degree,
    hole_ratio = hole_ratio,
    sparsity_index = sparsity,
    occupancy_ratio = occupancy_ratio,
    component_count = component_count
  )
  if (rotate && diagonal_chain) {
    compactness <- compactness - 0.10
  }
  if (sparsity > 0.25 && elongation > 2.4) {
    compactness <- compactness - 0.03
  }
  if (coverage_ratio < 0.2) {
    compactness <- compactness - 0.02
  }
  compactness <- max(0.4, min(0.8, compactness))

  out <- list(
    rows = as.integer(rows),
    cols = as.integer(cols),
    compactness = compactness,
    rotate_pca = rotate
  )

  if (embedding) {
    graph <- adjacency_graph(x)
    region <- if (rotate) rot_global else global
    emb_args <- utils::modifyList(
      list(
        iterations = 400,
        width = max(1e-9, region$width),
        height = max(1e-9, region$height)
      ),
      embedding_options
    )
    emb <- do.call(embed_graph, c(list(graph), emb_args))
    targets <- cbind(x = region$min_x + emb[, 1], y = region$min_y + emb[, 2])
    if (rotate) {
      # Express the targets in the original coordinate space.
      targets <- pca_rotate(targets[, 1], targets[, 2], angle, center)
    }
    out$adjacency_graph <- graph
    out$embedding_targets <- targets
  }

  out
}

#' Estimate grid parameters from point coordinates
#'
#' A simpler alternative to [estimate_parameters()] for plain points: grid
#' dimensions follow the aspect ratio of the points, compactness their
#' spread, and PCA rotation is recommended when it gives a tighter or more
#' elongated bounding box.
#'
#' @inheritParams allocate_grid
#'
#' @returns A list with elements `rows`, `cols`, `compactness` and
#'   `rotate_pca`.
#' @export
#' @examples
#' estimate_parameters_from_data(us_states, x = "lon", y = "lat")
estimate_parameters_from_data <- function(data, x = "x", y = "y") {
  xy <- extract_xy(data, x, y)
  px <- xy[, 1]
  py <- xy[, 2]
  b <- calc_bounds(px, py)
  dims <- grid_dimensions(nrow(xy), aspect_ratio = b$height / or_one(b$width))

  # Compactness from the spread of points around their mean
  dist <- sqrt(
    ((px - mean(px)) / or_one(b$width))^2 +
      ((py - mean(py)) / or_one(b$height))^2
  )
  spread <- min(1, sqrt(mean((dist - mean(dist))^2)) / 0.4)
  compactness <- max(0.2, min(0.8, 0.3 + (1 - spread) * 0.4))

  rotate <- FALSE
  if (length(px) >= 2) {
    angle <- pca_angle(px, py)
    rot <- pca_rotate(px, py, -angle)
    rb <- calc_bounds(rot[, 1], rot[, 2])
    area_o <- b$width * b$height
    area_r <- rb$width * rb$height
    aspect_o <- max(b$width, b$height) / min(b$width, b$height)
    aspect_r <- max(rb$width, rb$height) / min(rb$width, rb$height)
    rotate <- isTRUE(area_o > area_r * 1.1) ||
      isTRUE(aspect_r > aspect_o * 1.5) ||
      isTRUE(aspect_o > 1.5 && area_o > area_r * 1.05)
  }

  list(
    rows = dims[["rows"]],
    cols = dims[["cols"]],
    compactness = compactness,
    rotate_pca = rotate
  )
}

# Helpers -----------------------------------------------------------------------

bounds_of <- function(min_x, max_x, min_y, max_y) {
  b <- list(
    min_x = min(min_x),
    max_x = max(max_x),
    min_y = min(min_y),
    max_y = max(max_y)
  )
  b$width <- b$max_x - b$min_x
  b$height <- b$max_y - b$min_y
  b
}

# Feature area, falling back to the bounding box area (JS `area || bbox`).
feature_areas <- function(stats, fb) {
  a <- stats[, "area"]
  ifelse(is.na(a) | a == 0, fb[, "area"], a)
}

js_median <- function(x) {
  if (length(x) == 0) 0 else stats::median(x)
}

directional_indicators <- function(px, py, angle) {
  n <- length(px)
  if (n < 2) {
    return(list(linearity = 0, diagonality = 0))
  }
  dx <- px - mean(px)
  dy <- py - mean(py)
  cxx <- sum(dx * dx) / n
  cyy <- sum(dy * dy) / n
  cxy <- sum(dx * dy) / n
  tr <- cxx + cyy
  det <- cxx * cyy - cxy * cxy
  disc <- sqrt(max(0, tr * tr * 0.25 - det))
  l1 <- tr * 0.5 + disc
  l2 <- tr * 0.5 - disc
  linearity <- if (l1 + l2 > 1e-9) (l1 - l2) / (l1 + l2) else 0

  a <- abs(angle) %% pi
  if (a > pi / 2) {
    a <- pi - a
  }
  axis_dist <- min(a, abs(pi / 2 - a))
  diagonality <- sin(axis_dist * 2)

  list(
    linearity = max(0, min(1, linearity)),
    diagonality = max(0, min(1, diagonality))
  )
}

geo_indicators <- function(px, py, stats, fb) {
  n <- length(px)
  adj <- hybrid_adjacency_stats(px, py, fb)
  total_area <- sum(stats[, "area"])
  total_holes <- sum(stats[, "hole_area"])
  hull_area <- convex_hull_area(cbind(px, py))

  occupancy_ratio <- if (hull_area > 1e-9) min(1, total_area / hull_area) else 1
  hole_ratio <- if (total_area > 1e-9) {
    min(0.95, total_holes / total_area)
  } else {
    0
  }
  component_count <- max(1, adj$components)
  sparsity <- max(
    0,
    min(
      1,
      0.5 *
        (1 - occupancy_ratio) +
        0.25 * min(1, hole_ratio * 3) +
        0.25 * min(1, (component_count - 1) / max(1, sqrt(n)))
    )
  )

  list(
    average_degree = adj$average_degree,
    hole_ratio = hole_ratio,
    occupancy_ratio = occupancy_ratio,
    component_count = component_count,
    sparsity_index = sparsity
  )
}

# Degrees from bounding-box topology plus pruned k-nearest-neighbour
# proximity edges; components from topology only.
hybrid_adjacency_stats <- function(px, py, fb) {
  n <- length(px)
  if (n < 2) {
    return(list(average_degree = 0, components = 1))
  }
  topo <- bbox_adjacency(fb, buffer_fraction = 0.015)
  i <- topo$from
  j <- topo$to

  if (n >= 3) {
    prox <- proximity_edges(px, py, k = 4)
    keep <- prox$dist <= max(1e-9, js_median(prox$dist) * 2.2)
    i <- c(i, prox$i[keep])
    j <- c(j, prox$j[keep])
  }
  key <- unique(paste(pmin(i, j), pmax(i, j)))
  ends <- as.integer(unlist(strsplit(key, " ", fixed = TRUE)))
  degrees <- tabulate(ends, nbins = n)

  list(
    average_degree = sum(degrees) / n,
    components = count_components(n, topo$from, topo$to)
  )
}

proximity_edges <- function(px, py, k = 4) {
  n <- length(px)
  d <- sqrt(outer(px, px, "-")^2 + outer(py, py, "-")^2)
  kk <- max(1, min(k, n - 1))
  i <- integer()
  j <- integer()
  for (a in seq_len(n)) {
    others <- setdiff(seq_len(n), a)
    nb <- others[order(d[a, others])][seq_len(kk)]
    i <- c(i, rep(a, kk))
    j <- c(j, nb)
  }
  lo <- pmin(i, j)
  hi <- pmax(i, j)
  dup <- duplicated(paste(lo, hi))
  data.frame(i = lo[!dup], j = hi[!dup], dist = d[cbind(lo[!dup], hi[!dup])])
}

count_components <- function(n, from, to) {
  parent <- seq_len(n)
  find <- function(v) {
    while (parent[v] != v) {
      v <- parent[v]
    }
    v
  }
  for (e in seq_along(from)) {
    a <- find(from[e])
    b <- find(to[e])
    if (a != b) parent[a] <- b
  }
  length(unique(vapply(seq_len(n), find, integer(1))))
}

refine_dimensions <- function(
  base_rows,
  base_cols,
  n,
  target_cells,
  target_aspect,
  min_rows,
  min_cols,
  diagonal_chain,
  coverage_ratio,
  sparsity_index,
  elongation
) {
  cand <- list()
  add <- function(r, c) {
    r <- max(min_rows, js_round(r))
    c <- max(min_cols, js_round(c))
    if (r < 1 || c < 1 || r * c < n) {
      return()
    }
    cand[[paste0(r, "x", c)]] <<- c(rows = r, cols = c)
  }
  add(base_rows, base_cols)
  for (dr in -2:2) {
    for (dc in -2:2) {
      if (dr != 0 || dc != 0) add(base_rows + dr, base_cols + dc)
    }
  }
  safe_aspect <- max(1e-6, target_aspect)
  cells <- max(target_cells, n)
  c_aspect <- max(min_cols, js_round(sqrt(cells / safe_aspect)))
  r_aspect <- max(min_rows, ceiling(cells / c_aspect))
  add(r_aspect, c_aspect)
  add(r_aspect - 1, c_aspect + 1)
  add(r_aspect + 1, c_aspect - 1)
  add(r_aspect + 1, c_aspect + 1)

  ctx <- list(
    target_cells = target_cells,
    target_aspect = target_aspect,
    n = n,
    compact_map = coverage_ratio > 0.48 && elongation < 2.4 && !diagonal_chain,
    diagonal_chain = diagonal_chain,
    sparse_map = sparsity_index > 0.32,
    elongation = elongation,
    sparsity_index = sparsity_index
  )

  best <- c(rows = base_rows, cols = base_cols)
  best_score <- score_dimensions(best, ctx)
  for (cd in cand) {
    s <- score_dimensions(cd, ctx)
    if (s < best_score) {
      best_score <- s
      best <- cd
    }
  }
  best
}

score_dimensions <- function(cand, ctx) {
  rows <- cand[["rows"]]
  cols <- cand[["cols"]]
  n <- ctx$n
  cells <- rows * cols
  ratio <- rows / max(cols, 1e-9)
  target_aspect <- max(1e-9, ctx$target_aspect)
  target_cells <- ctx$target_cells

  score <- abs(log(max(ratio, 1e-9) / target_aspect)) *
    0.6 +
    abs(cells - target_cells) / max(1, target_cells) * 0.28 +
    max(0, cells - target_cells) / max(1, target_cells) * 0.12

  if (ctx$compact_map) {
    imbalance <- abs(cols - rows) / max(1, sqrt(n))
    score <- score + max(0, imbalance - 0.35) * 0.35
  }
  if (ctx$diagonal_chain) {
    score <- score + max(0, ratio - min(0.75, target_aspect * 1.25)) * 0.8
  }
  if (ctx$sparse_map) {
    score <- score + max(0, ratio - 0.78) * 0.3
  }
  # "Chain-like" maps are the diagonal chains; the JS code also adds a
  # non-diagonal branch here that can never trigger.
  if (ctx$diagonal_chain) {
    long_axis <- max(rows, cols)
    elong_boost <- max(0, log2(max(ctx$elongation, 1)))
    spread_target <- sqrt(max(n, 1)) *
      (1.35 + 0.35 * elong_boost + 0.25 * max(0, ctx$sparsity_index - 0.25))
    score <- score +
      max(0, spread_target - long_axis) / max(1, spread_target) * 0.9
    lane_target <- 3.2 - min(1.2, elong_boost * 0.35)
    lane_density <- n / max(1, long_axis)
    score <- score +
      max(0, lane_density - lane_target) / max(1, lane_target) * 0.6
  }
  score
}

estimate_compactness_advanced <- function(
  px,
  py,
  feature_area,
  bounds,
  average_degree,
  hole_ratio,
  sparsity_index,
  occupancy_ratio,
  component_count
) {
  n <- length(px)
  if (n < 2) {
    return(0.5)
  }
  total_area <- bounds$width * bounds$height
  if (is.na(total_area) || total_area == 0) {
    total_area <- 1
  }
  pb <- calc_bounds(px, py)
  aspect <- max(pb$width, pb$height) / max(1e-9, min(pb$width, pb$height))

  # Nearest neighbour index: < 1 clustered, > 1 dispersed
  d <- sqrt(outer(px, px, "-")^2 + outer(py, py, "-")^2)
  diag(d) <- Inf
  observed <- sum(apply(d, 1, min)) / n
  expected <- 0.5 / sqrt(n / total_area)
  nni <- if (expected > 0) observed / expected else 1

  coverage <- min(1, sum(feature_area) / total_area)

  comp <- 0.55
  if (nni < 0.7) {
    comp <- comp + 0.05
  } else if (nni > 1.2) {
    comp <- comp - 0.03
  }
  if (coverage > 0.5) {
    comp <- comp + 0.03
  } else if (coverage < 0.1) {
    comp <- comp + 0.05
  }
  if (sparsity_index > 0.45) {
    comp <- comp - 0.07
  } else if (sparsity_index > 0.3) {
    comp <- comp - 0.04
  }
  if (hole_ratio > 0.12) {
    comp <- comp - 0.04
  }
  if (occupancy_ratio < 0.32) {
    comp <- comp - 0.03
  }
  if (component_count > 1) {
    comp <- comp - min(0.06, (component_count - 1) * 0.015)
  }
  if (average_degree > 4) {
    comp <- comp + 0.02
  } else if (average_degree < 1.5) {
    comp <- comp + 0.01
  }

  if (aspect > 2.3) {
    weight <- min(1, (aspect - 2.3) / 2)
    comp <- 0.5 + (comp - 0.5) * (1 - weight)
    comp <- max(0.45, min(0.65, comp))
  }
  max(0.35, min(0.85, comp))
}
