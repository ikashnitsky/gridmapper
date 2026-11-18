#' Allocate geographic features to a grid
#'
#' Assigns each feature to one cell of a rectangular or hexagonal grid by
#' solving a mixed integer programme (MIP) with [HiGHS][highs::highs_solve].
#' The objective minimises the distance between each feature's
#' (compactness-scaled) position and its cell, optionally rewarding occupied
#' neighbouring cells and adjacent feature pairs.
#'
#' Grid cells are indexed from 1. Row 1 is at the bottom (minimum `y`) and
#' column 1 at the left (minimum `x`). In hexagonal grids even rows are
#' shifted right by half a cell.
#'
#' @param data Features to allocate. One of:
#'   * a data frame with coordinate columns named by `x` and `y`;
#'   * a numeric matrix with two columns;
#'   * an [sf][sf::st_sf] object, in which case each feature is represented
#'     by its [polygon_centroid()].
#' @param x,y Names of the coordinate columns when `data` is a data frame.
#' @param rows,cols Grid dimensions. If either is `NULL`, both are chosen
#'   with [grid_dimensions()] from the aspect ratio of the data. Use
#'   [estimate_parameters()] for a more elaborate choice.
#' @param compactness Number between 0 and 1. At `0.5` features keep their
#'   relative geographic position; larger values pull them towards the
#'   centre of the grid, smaller values push them to the edges.
#' @param grid_type Grid geometry, `"rect"` or `"hex"`.
#' @param distance Distance in the objective: squared `"euclidean"` or
#'   `"manhattan"`.
#' @param rotate_pca If `TRUE`, rotate the coordinates so that their first
#'   principal axis is horizontal before allocating. Useful for diagonal
#'   geographies.
#' @param spacers Cells to keep empty: a two-column matrix or data frame of
#'   1-based `row` and `col` positions.
#' @param mask Optional polygon used to derive spacers when `spacers` is
#'   `NULL`; cells outside it become spacers. See [compute_spacers()]. The
#'   mask is rotated together with the data when `rotate_pca = TRUE`.
#' @param spacer_mode `"hard"` excludes spacer cells from the grid; `"soft"`
#'   keeps them available but adds `mask_penalty` to their cost.
#' @param mask_penalty Cost added for placing a feature in a soft spacer.
#' @param compactness_weight Multiplier of the distance cost.
#' @param adjacency_weight Reward for each pair of neighbouring occupied
#'   cells. Values above 0 make the layout more contiguous but enlarge the
#'   MIP.
#' @param adjacency_diagonal If `TRUE`, diagonal cells count as neighbours.
#' @param embedding_targets Optional two-column matrix (one row per feature)
#'   of target coordinates, in the same coordinate space as the data, e.g.
#'   the `embedding_targets` from [estimate_parameters()]. Each feature is
#'   additionally pulled towards its target by `embedding_weight`.
#' @param embedding_weight Weight of the embedding cost.
#' @param adjacency_graph Optional feature adjacency graph from
#'   [adjacency_graph()]. Needed by `pairwise_adjacency` and by
#'   `post_process_options(fix_adjacency = TRUE)`.
#' @param pairwise_adjacency Optional [pairwise_adjacency_options()]
#'   rewarding adjacent features placed in adjacent cells. Requires
#'   `adjacency_graph`.
#' @param post_process Optional heuristic refinement after solving: `TRUE`
#'   for the defaults of [post_process_options()], or a call to it.
#' @param solver_options A list of options passed to
#'   [highs::highs_control()], e.g. `list(time_limit = 10)`.
#'
#' @returns A `gridmap` object: a list with elements
#'   * `assignments`: a data frame with the input coordinates (`x`, `y`) and
#'     the assigned `row` and `col` of each feature;
#'   * `data`: the input `data`;
#'   * `rows`, `cols`, `grid_type`;
#'   * `score`: the objective value, `status`: the solver status, and
#'     `feasible`: whether an optimal solution was found;
#'   * `bounds`: the bounding box of the (unrotated) input coordinates;
#'   * `spacers`: the spacer cells used;
#'   * `post_process`: diagnostics of the post-processing step, if any.
#'
#'   Use [as.data.frame()] to bind the assignments to the input data and
#'   [grid_to_sf()] to turn them into polygons.
#' @seealso [allocate_grid_simple()] for the basic interface of the original
#'   method.
#' @export
#' @examples
#' gm <- allocate_grid(us_states, x = "lon", y = "lat", compactness = 0.6)
#' gm
#' head(as.data.frame(gm))
#' plot(gm, labels = "abbr")
allocate_grid <- function(
  data,
  x = "x",
  y = "y",
  rows = NULL,
  cols = NULL,
  compactness = 0.5,
  grid_type = c("rect", "hex"),
  distance = c("euclidean", "manhattan"),
  rotate_pca = FALSE,
  spacers = NULL,
  mask = NULL,
  spacer_mode = c("hard", "soft"),
  mask_penalty = 1000,
  compactness_weight = 1,
  adjacency_weight = 0,
  adjacency_diagonal = FALSE,
  embedding_targets = NULL,
  embedding_weight = 0,
  adjacency_graph = NULL,
  pairwise_adjacency = NULL,
  post_process = NULL,
  solver_options = list()
) {
  grid_type <- rlang::arg_match(grid_type)
  distance <- rlang::arg_match(distance)
  spacer_mode <- rlang::arg_match(spacer_mode)
  check_count(rows, allow_null = TRUE)
  check_count(cols, allow_null = TRUE)
  check_number(compactness, min = 0, max = 1)
  check_flag(rotate_pca)
  check_number(mask_penalty)
  check_number(compactness_weight)
  check_number(adjacency_weight)
  check_flag(adjacency_diagonal)
  check_number(embedding_weight)

  xy <- extract_xy(data, x, y)
  n <- nrow(xy)
  bounds <- calc_bounds(xy[, 1], xy[, 2])

  if (is.null(rows) || is.null(cols)) {
    dims <- grid_dimensions(
      n,
      aspect_ratio = bounds$height / or_one(bounds$width)
    )
    rows <- dims[["rows"]]
    cols <- dims[["cols"]]
  }

  work <- xy
  emb <- check_embedding(embedding_targets, n)
  ring <- mask_ring(mask)
  if (rotate_pca) {
    angle <- pca_angle(xy[, 1], xy[, 2])
    center <- colMeans(xy)
    work <- pca_rotate(xy[, 1], xy[, 2], -angle, center)
    if (!is.null(emb)) {
      emb <- pca_rotate(emb[, 1], emb[, 2], -angle, center)
    }
    if (!is.null(ring)) {
      ring <- pca_rotate(ring[, 1], ring[, 2], -angle, center)
    }
  }
  work_bounds <- calc_bounds(work[, 1], work[, 2])

  keys <- if (!is.null(spacers)) {
    unique(spacer_keys(spacers))
  } else if (!is.null(ring)) {
    spacer_keys_from_ring(work_bounds, rows, cols, ring, grid_type)
  } else {
    character()
  }

  cells <- grid_cells(rows, cols)
  cell_keys <- paste0(cells$row, "_", cells$col)
  is_spacer <- cell_keys %in% keys
  n_spacers <- sum(is_spacer)
  if (spacer_mode == "hard") {
    cells <- cells[!is_spacer, , drop = FALSE]
    is_spacer <- is_spacer[!is_spacer]
  }
  if (nrow(cells) < n) {
    cli::cli_abort(c(
      "Can't allocate {n} feature{?s} to {nrow(cells)} available cell{?s}.",
      i = "The {rows} x {cols} grid has {n_spacers} spacer{?s}.",
      i = "Increase {.arg rows} or {.arg cols}, or reduce the spacers."
    ))
  }

  norm <- normalize_points(
    work[, 1],
    work[, 2],
    work_bounds,
    rows,
    cols,
    compactness,
    grid_type
  )
  norm_emb <- if (!is.null(emb)) {
    normalize_points(
      emb[, 1],
      emb[, 2],
      work_bounds,
      rows,
      cols,
      compactness,
      grid_type
    )
  }

  pairwise <- NULL
  if (!is.null(pairwise_adjacency)) {
    if (!inherits(pairwise_adjacency, "gridmapper_pairwise_options")) {
      cli::cli_abort(
        "{.arg pairwise_adjacency} must be created by {.fn pairwise_adjacency_options}."
      )
    }
    check_adjacency_graph(adjacency_graph, n, "pairwise_adjacency")
    pairwise <- c(
      unclass(pairwise_adjacency),
      list(edges = adjacency_graph$edges)
    )
  }

  solution <- solve_allocation(
    norm$x,
    norm$y,
    cells,
    soft_spacer = is_spacer & spacer_mode == "soft",
    grid_type = grid_type,
    distance = distance,
    compactness_weight = compactness_weight,
    mask_penalty = mask_penalty,
    adjacency_weight = adjacency_weight,
    adjacency_diagonal = adjacency_diagonal,
    embedding = norm_emb,
    embedding_weight = embedding_weight,
    pairwise = pairwise,
    solver_options = solver_options
  )

  row0 <- cells$row[solution$cell]
  col0 <- cells$col[solution$cell]

  pp_meta <- NULL
  if (isTRUE(post_process)) {
    post_process <- post_process_options()
  }
  if (inherits(post_process, "gridmapper_post_process_options")) {
    if (post_process$fix_adjacency) {
      check_adjacency_graph(adjacency_graph, n, "fix_adjacency")
    }
    pp <- post_process_assignments(
      row0,
      col0,
      norm,
      rows,
      cols,
      opts = post_process,
      grid_type = grid_type,
      distance = distance,
      compactness_weight = compactness_weight,
      adjacency_weight = adjacency_weight,
      adjacency_diagonal = adjacency_diagonal,
      spacer_mode = spacer_mode,
      mask_penalty = mask_penalty,
      spacer_keys = keys,
      adjacency_graph = adjacency_graph
    )
    row0 <- pp$row
    col0 <- pp$col
    pp_meta <- pp$meta
  } else if (!is.null(post_process) && !isFALSE(post_process)) {
    cli::cli_abort(
      "{.arg post_process} must be `TRUE`, `FALSE`, `NULL` or created by {.fn post_process_options}."
    )
  }

  structure(
    list(
      assignments = data.frame(
        x = unname(xy[, 1]),
        y = unname(xy[, 2]),
        row = row0 + 1L,
        col = col0 + 1L
      ),
      data = data,
      rows = as.integer(rows),
      cols = as.integer(cols),
      grid_type = grid_type,
      score = solution$score,
      status = solution$status,
      feasible = solution$feasible,
      bounds = bounds,
      spacers = keys_to_spacer_matrix(keys),
      pairwise_vars = solution$pairwise_vars,
      post_process = pp_meta
    ),
    class = "gridmap"
  )
}

#' Allocate points to a grid (simple interface)
#'
#' The basic grid map allocation of Jo Wood's original method: points are
#' assigned to the cells of a rectangular grid, minimising the squared
#' distance between each point's compactness-scaled position and its cell.
#'
#' @param points A two-column numeric matrix or data frame of `x` and `y`
#'   coordinates.
#' @param n_rows,n_cols Grid dimensions.
#' @param compactness Number between 0 and 1; see [allocate_grid()].
#' @param spacers Cells to keep empty: a two-column matrix or data frame of
#'   1-based `row` and `col` positions.
#' @inheritParams allocate_grid
#'
#' @returns A list with elements `n_rows`, `n_cols` and `cells`, an integer
#'   matrix with columns `row` and `col` giving the cell of each point, in
#'   input order.
#' @export
#' @examples
#' pts <- as.matrix(us_states[c("lon", "lat")])
#' res <- allocate_grid_simple(pts, n_rows = 6, n_cols = 10, compactness = 0.6)
#' head(res$cells)
allocate_grid_simple <- function(
  points,
  n_rows,
  n_cols,
  compactness = 1,
  spacers = NULL,
  solver_options = list()
) {
  check_count(n_rows)
  check_count(n_cols)
  check_number(compactness, min = 0, max = 1)
  xy <- extract_xy(as.matrix(points))

  keys <- unique(spacer_keys(spacers))
  cells <- grid_cells(n_rows, n_cols)
  cells <- cells[!paste0(cells$row, "_", cells$col) %in% keys, , drop = FALSE]
  if (nrow(xy) > nrow(cells)) {
    cli::cli_abort(
      "Can't allocate {nrow(xy)} point{?s} to a grid with only {nrow(cells)} cell{?s}."
    )
  }

  bounds <- calc_bounds(xy[, 1], xy[, 2])
  norm <- normalize_points(
    xy[, 1],
    xy[, 2],
    bounds,
    n_rows,
    n_cols,
    compactness,
    "rect"
  )
  solution <- solve_allocation(
    norm$x,
    norm$y,
    cells,
    soft_spacer = rep(FALSE, nrow(cells)),
    solver_options = solver_options
  )

  out <- cbind(
    row = cells$row[solution$cell] + 1L,
    col = cells$col[solution$cell] + 1L
  )
  list(n_rows = as.integer(n_rows), n_cols = as.integer(n_cols), cells = out)
}

#' Options for the pairwise adjacency reward
#'
#' Configures the sparse pairwise term of [allocate_grid()]: for the
#' strongest edges of the adjacency graph, an auxiliary binary variable
#' rewards placing the two features in neighbouring cells. To keep the MIP
#' small, only cell pairs among the `nearest_cells` closest cells of each
#' feature are considered.
#'
#' @param edge_limit Maximum number of edges (strongest first) to use.
#'   `NULL` uses all edges.
#' @param nearest_cells Number of candidate cells per feature.
#' @param weight Reward multiplier; each edge contributes
#'   `weight * edge weight`.
#'
#' @returns A list of class `gridmapper_pairwise_options`.
#' @export
#' @examples
#' pairwise_adjacency_options(edge_limit = 100, nearest_cells = 6)
pairwise_adjacency_options <- function(
  edge_limit = NULL,
  nearest_cells = 6,
  weight = 1
) {
  check_count(edge_limit, allow_null = TRUE)
  check_count(nearest_cells)
  check_number(weight)
  structure(
    list(
      edge_limit = edge_limit,
      nearest_cells = nearest_cells,
      weight = weight
    ),
    class = "gridmapper_pairwise_options"
  )
}

# Methods --------------------------------------------------------------------

#' @export
print.gridmap <- function(x, ...) {
  n <- nrow(x$assignments)
  cli::cli_text(
    "<gridmap> {n} feature{?s} on a {x$rows} x {x$cols} {x$grid_type} grid"
  )
  cli::cli_text("Score: {format(x$score, digits = 6)} ({x$status})")
  if (nrow(x$spacers) > 0) {
    cli::cli_text("Spacers: {nrow(x$spacers)}")
  }
  if (!is.null(x$post_process)) {
    cli::cli_text(
      "Post-processing: {x$post_process$swaps} swap{?s} in {x$post_process$iterations} iteration{?s}"
    )
  }
  invisible(x)
}

#' Convert a gridmap to a data frame
#'
#' Binds the assigned `row` and `col` of each feature to the input data.
#'
#' @param x A `gridmap` object from [allocate_grid()].
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#'
#' @returns A data frame with the columns of the input data (without
#'   geometry for `sf` input) followed by `row` and `col`.
#' @export
#' @examples
#' gm <- allocate_grid(us_states, x = "lon", y = "lat")
#' head(as.data.frame(gm))
as.data.frame.gridmap <- function(x, row.names = NULL, optional = FALSE, ...) {
  data <- input_attributes(x$data)
  cbind(data, x$assignments[c("row", "col")])
}

input_attributes <- function(data) {
  if (inherits(data, "sf")) {
    data <- sf::st_drop_geometry(data)
  }
  if (is.data.frame(data)) {
    data <- as.data.frame(data)
    data$row <- NULL
    data$col <- NULL
    return(data)
  }
  data.frame(id = seq_len(NROW(data)))[0]
}

#' Plot a gridmap
#'
#' Draws the grid cells of an allocation, optionally labelled.
#'
#' @param x A `gridmap` object from [allocate_grid()].
#' @param labels Optional name of a column of the input data to use as cell
#'   labels.
#' @param col Fill colour of the cells.
#' @param border Border colour of the cells.
#' @param cex Size of the labels.
#' @param ... Passed on to [plot()][sf::plot.sf] for the cell geometries.
#'
#' @returns `x`, invisibly.
#' @export
#' @examples
#' gm <- allocate_grid(us_states, x = "lon", y = "lat", grid_type = "hex")
#' plot(gm, labels = "abbr")
plot.gridmap <- function(
  x,
  labels = NULL,
  col = "grey90",
  border = "grey40",
  cex = 0.7,
  ...
) {
  tiles <- grid_to_sf(x)
  plot(sf::st_geometry(tiles), col = col, border = border, ...)
  if (!is.null(labels)) {
    lab <- as.data.frame(x)[[labels]]
    if (is.null(lab)) {
      cli::cli_abort("Column {.field {labels}} not found in the input data.")
    }
    centers <- tile_centers(x)
    graphics::text(centers[, 1], centers[, 2], labels = lab, cex = cex)
  }
  invisible(x)
}

# Helpers ---------------------------------------------------------------------

extract_xy <- function(data, x = "x", y = "y", call = rlang::caller_env()) {
  if (inherits(data, c("sf", "sfc"))) {
    xy <- polygon_centroid(data)
  } else if (is.matrix(data)) {
    if (!is.numeric(data) || ncol(data) < 2) {
      cli::cli_abort(
        "A matrix {.arg data} must be numeric with two columns.",
        call = call
      )
    }
    xy <- data[, 1:2, drop = FALSE]
  } else if (is.data.frame(data)) {
    missing <- setdiff(c(x, y), names(data))
    if (length(missing)) {
      cli::cli_abort(
        "Column{?s} {.field {missing}} not found in {.arg data}.",
        call = call
      )
    }
    xy <- cbind(data[[x]], data[[y]])
  } else {
    cli::cli_abort(
      "{.arg data} must be a data frame, matrix or {.cls sf} object, not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  xy <- unname(xy)
  storage.mode(xy) <- "double"
  if (nrow(xy) == 0) {
    cli::cli_abort(
      "{.arg data} must contain at least one feature.",
      call = call
    )
  }
  if (anyNA(xy) || any(!is.finite(xy))) {
    cli::cli_abort("Coordinates must be finite and not missing.", call = call)
  }
  colnames(xy) <- c("x", "y")
  xy
}

check_embedding <- function(targets, n, call = rlang::caller_env()) {
  if (is.null(targets)) {
    return(NULL)
  }
  targets <- as.matrix(targets)
  if (!is.numeric(targets) || ncol(targets) < 2 || nrow(targets) != n) {
    cli::cli_abort(
      "{.arg embedding_targets} must be a two-column numeric matrix with one row per feature.",
      call = call
    )
  }
  targets[, 1:2, drop = FALSE]
}

check_adjacency_graph <- function(graph, n, what, call = rlang::caller_env()) {
  if (!inherits(graph, "gridmapper_adjacency")) {
    cli::cli_abort(
      "{.arg {what}} requires {.arg adjacency_graph}, created by {.fn adjacency_graph}.",
      call = call
    )
  }
  if (nrow(graph$nodes) != n) {
    cli::cli_abort(
      c(
        "{.arg adjacency_graph} must have one node per feature.",
        x = "It has {nrow(graph$nodes)} node{?s}; the data has {n} feature{?s}."
      ),
      call = call
    )
  }
  invisible(graph)
}
