HEX_Y <- 0.8660254037844386

# Centres of 0-based grid cells in grid coordinate space. Hex grids are
# "odd-r": odd (0-based) rows are shifted right by half a cell.
cell_center <- function(row, col, grid_type = "rect") {
  if (grid_type == "hex") {
    list(x = col + ifelse(row %% 2 == 1, 0.5, 0), y = row * HEX_Y)
  } else {
    list(x = col, y = row)
  }
}

# All 0-based cells of a grid in row-major order.
grid_cells <- function(rows, cols) {
  data.frame(
    row = rep(seq_len(rows) - 1L, each = cols),
    col = rep(seq_len(cols) - 1L, times = rows)
  )
}

neighbour_offsets <- function(diagonal = FALSE) {
  off <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  if (diagonal) {
    off <- rbind(off, c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  }
  off
}

#' Choose grid dimensions for a number of features
#'
#' Picks the number of rows and columns of a grid with at least
#' `target_cells` cells whose shape approximates `aspect_ratio`. This is the
#' rule [allocate_grid()] uses when `rows` or `cols` are not given.
#'
#' @param n Number of features to allocate.
#' @param aspect_ratio Desired ratio of height to width.
#' @param target_cells Desired number of cells. Never less than `n`.
#' @param min_rows,min_cols Minimum number of rows and columns.
#'
#' @returns A named integer vector with elements `rows` and `cols`.
#' @export
#' @examples
#' grid_dimensions(50, aspect_ratio = 0.6)
#' grid_dimensions(50, aspect_ratio = 0.6, target_cells = 70)
grid_dimensions <- function(
  n,
  aspect_ratio = 1,
  target_cells = n,
  min_rows = 1,
  min_cols = 1
) {
  check_count(n)
  check_number(aspect_ratio)
  check_number(target_cells)

  safe_aspect <- if (aspect_ratio > 0) aspect_ratio else 1
  cells <- max(target_cells, n, 1)

  cols <- max(1, js_round(sqrt(cells / safe_aspect)))
  rows <- ceiling(cells / cols)
  rows <- max(rows, min_rows)
  cols <- max(cols, min_cols)

  while (rows * cols < cells) {
    if (rows / safe_aspect <= cols) {
      rows <- rows + 1
    } else {
      cols <- cols + 1
    }
  }
  while (rows * cols < n) {
    if (rows <= cols) {
      rows <- rows + 1
    } else {
      cols <- cols + 1
    }
  }

  c(rows = as.integer(rows), cols = as.integer(cols))
}

# Normalisation --------------------------------------------------------------

# Scales points into the "ideal" floating-point grid position. Compactness
# controls the size of the target rectangle: 1 clusters in the centre, 0.5
# keeps the relative geographic position, 0 spreads to the edges.
normalize_points <- function(x, y, bounds, rows, cols, compactness, grid_type) {
  stretch <- 1 / (compactness + 0.001) - 1
  if (grid_type == "hex") {
    eff_cols <- (cols - 1) + 0.5
    eff_rows <- (rows - 1) * HEX_Y
    cc <- eff_cols / 2
    cr <- eff_rows / 2
    rw <- stretch * eff_cols + 1
    rh <- stretch * eff_rows + 1
  } else {
    cc <- (cols - 1) / 2
    cr <- (rows - 1) / 2
    rw <- stretch * (cols - 1) + 1
    rh <- stretch * (rows - 1) + 1
  }
  nx <- (x - bounds$min_x) / or_one(bounds$width)
  ny <- (y - bounds$min_y) / or_one(bounds$height)
  list(x = cc - rw / 2 + nx * rw, y = cr - rh / 2 + ny * rh)
}

# PCA ------------------------------------------------------------------------

#' Rotate points onto their principal axis
#'
#' `pca_angle()` returns the angle (in radians) of the first principal axis
#' of a point cloud. `pca_rotate()` rotates points by an angle around a
#' centre. [allocate_grid()] with `rotate_pca = TRUE` rotates by
#' `-pca_angle()` so that the principal axis lies along the grid's x axis.
#'
#' @param x,y Numeric vectors of coordinates.
#' @param angle Rotation angle in radians (counter-clockwise).
#' @param center Optional numeric vector `c(x, y)` giving the centre of
#'   rotation. Defaults to the mean of the points.
#'
#' @returns `pca_angle()` returns a single number. `pca_rotate()` returns a
#'   matrix with columns `x` and `y`.
#' @export
#' @examples
#' x <- c(0, 1, 2, 3)
#' y <- c(0, 1, 2, 3)
#' angle <- pca_angle(x, y)
#' angle
#' pca_rotate(x, y, -angle)
pca_angle <- function(x, y) {
  if (length(x) < 2) {
    return(0)
  }
  dx <- x - mean(x)
  dy <- y - mean(y)
  0.5 * atan2(2 * sum(dx * dy), sum(dx * dx) - sum(dy * dy))
}

#' @rdname pca_angle
#' @export
pca_rotate <- function(x, y, angle, center = NULL) {
  if (is.null(center)) {
    center <- c(mean(x), mean(y))
  }
  cos_a <- cos(angle)
  sin_a <- sin(angle)
  dx <- x - center[1]
  dy <- y - center[2]
  cbind(
    x = center[1] + dx * cos_a - dy * sin_a,
    y = center[2] + dx * sin_a + dy * cos_a
  )
}

# Spacers --------------------------------------------------------------------

#' Find grid cells that fall outside a mask
#'
#' Maps every cell of a `rows` x `cols` grid onto the bounding box of the
#' points and marks as a spacer each cell whose position lies outside `mask`.
#' Spacer cells are left empty (or penalised) by [allocate_grid()].
#'
#' @param x,y Numeric vectors of point coordinates.
#' @param rows,cols Grid dimensions.
#' @param mask Polygon outside which cells become spacers: an
#'   [sf][sf::st_sf] / [sfc][sf::st_sfc] object (the outer ring of the first
#'   polygon is used), or a two-column matrix of ring coordinates. `NULL`
#'   (the default) uses the convex hull of the points.
#' @param grid_type Grid geometry, `"rect"` or `"hex"`.
#'
#' @returns An integer matrix with columns `row` and `col` (1-based) of
#'   spacer cells.
#' @export
#' @examples
#' set.seed(1)
#' x <- runif(30)
#' y <- runif(30)
#' compute_spacers(x, y, rows = 7, cols = 7)
compute_spacers <- function(
  x,
  y,
  rows,
  cols,
  mask = NULL,
  grid_type = c("rect", "hex")
) {
  grid_type <- rlang::arg_match(grid_type)
  check_count(rows)
  check_count(cols)
  bounds <- calc_bounds(x, y)
  ring <- mask_ring(mask) %||% convex_hull(cbind(x, y))
  keys <- spacer_keys_from_ring(bounds, rows, cols, ring, grid_type)
  keys_to_spacer_matrix(keys)
}

mask_ring <- function(mask, call = rlang::caller_env()) {
  if (is.null(mask)) {
    return(NULL)
  }
  if (inherits(mask, c("sf", "sfc", "sfg"))) {
    geom <- as_geometry_list(mask)[[1]]
    ring <- switch(
      geom$type,
      Polygon = geom$coords[[1]],
      MultiPolygon = geom$coords[[1]][[1]],
      cli::cli_abort("{.arg mask} must be a polygon.", call = call)
    )
    return(ring)
  }
  if (is.matrix(mask) && is.numeric(mask) && ncol(mask) >= 2) {
    return(mask[, 1:2, drop = FALSE])
  }
  cli::cli_abort(
    "{.arg mask} must be an {.cls sf} polygon or a two-column matrix.",
    call = call
  )
}

spacer_keys_from_ring <- function(bounds, rows, cols, ring, grid_type) {
  cells <- grid_cells(rows, cols)
  r <- cells$row
  c <- cells$col
  if (grid_type == "hex") {
    offset <- ifelse(r %% 2 == 1, 0.5, 0)
    tx <- if (cols > 1) {
      (c + offset) / ((cols - 1) + 0.5)
    } else {
      rep(0.5, nrow(cells))
    }
    ty <- if (rows > 1) {
      (r * HEX_Y) / ((rows - 1) * HEX_Y)
    } else {
      rep(0.5, nrow(cells))
    }
  } else {
    tx <- if (cols > 1) c / (cols - 1) else rep(0.5, nrow(cells))
    ty <- if (rows > 1) r / (rows - 1) else rep(0.5, nrow(cells))
  }
  gx <- bounds$min_x + tx * bounds$width
  gy <- bounds$min_y + ty * bounds$height
  inside <- vapply(
    seq_along(gx),
    \(i) point_in_polygon(gx[i], gy[i], ring),
    logical(1)
  )
  paste0(r, "_", c)[!inside]
}
