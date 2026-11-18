#' Convert a gridmap to grid cell polygons
#'
#' Builds one square or hexagonal tile per feature, laid out within the
#' bounding box of the input coordinates, so the cartogram can be drawn or
#' exported alongside the original map.
#'
#' @param x A `gridmap` object from [allocate_grid()].
#' @param bounds Optional bounding box to lay the grid out in: a list or
#'   named numeric vector with `min_x`, `max_x`, `min_y` and `max_y`.
#'   Defaults to the bounds of the input coordinates.
#' @param crs Coordinate reference system of the result. Defaults to that of
#'   the input data if it is an `sf` object.
#'
#' @returns An [sf][sf::st_sf] object with one polygon per feature and the
#'   columns of [as.data.frame.gridmap()].
#' @export
#' @examples
#' gm <- allocate_grid(us_states, x = "lon", y = "lat", grid_type = "hex")
#' tiles <- grid_to_sf(gm)
#' plot(tiles["region"])
grid_to_sf <- function(x, bounds = NULL, crs = NULL) {
  if (!inherits(x, "gridmap")) {
    cli::cli_abort("{.arg x} must be a {.cls gridmap} object.")
  }
  layout <- tile_layout(x, bounds)
  if (is.null(crs)) {
    crs <- if (inherits(x$data, "sf")) sf::st_crs(x$data) else sf::NA_crs_
  }

  centers <- tile_centers(x, bounds)
  polys <- lapply(seq_len(nrow(centers)), function(i) {
    ring <- if (x$grid_type == "hex") {
      hexagon_ring(centers[i, 1], centers[i, 2], layout$radius)
    } else {
      cell_ring(x$assignments$row[i] - 1, x$assignments$col[i] - 1, layout)
    }
    sf::st_polygon(list(ring))
  })

  sf::st_sf(
    as.data.frame(x),
    geometry = sf::st_sfc(polys, crs = crs)
  )
}

tile_layout <- function(x, bounds = NULL) {
  b <- bounds %||% x$bounds
  b <- as.list(b)
  width <- b$max_x - b$min_x
  height <- b$max_y - b$min_y
  if (x$grid_type == "hex") {
    scale_x <- width / (x$cols - 1 + 0.5)
    scale_y <- height / ((x$rows - 1) * HEX_Y)
    scale <- min(scale_x, scale_y)
    total_w <- (x$cols - 1 + 0.5) * scale
    total_h <- (x$rows - 1) * HEX_Y * scale
    list(
      scale = scale,
      radius = scale / 2,
      offset_x = b$min_x + (width - total_w) / 2,
      offset_y = b$min_y + (height - total_h) / 2
    )
  } else {
    list(
      min_x = b$min_x,
      min_y = b$min_y,
      col_w = width / x$cols,
      row_h = height / x$rows
    )
  }
}

tile_centers <- function(x, bounds = NULL) {
  layout <- tile_layout(x, bounds)
  r <- x$assignments$row - 1
  c <- x$assignments$col - 1
  if (x$grid_type == "hex") {
    ctr <- cell_center(r, c, "hex")
    cbind(
      x = layout$offset_x + ctr$x * layout$scale,
      y = layout$offset_y + ctr$y * layout$scale
    )
  } else {
    cbind(
      x = layout$min_x + (c + 0.5) * layout$col_w,
      y = layout$min_y + (r + 0.5) * layout$row_h
    )
  }
}

# Rectangular tile computed from its edges, so neighbours share them exactly.
cell_ring <- function(r, c, layout) {
  x0 <- layout$min_x + c * layout$col_w
  x1 <- layout$min_x + (c + 1) * layout$col_w
  y0 <- layout$min_y + r * layout$row_h
  y1 <- layout$min_y + (r + 1) * layout$row_h
  cbind(c(x0, x1, x1, x0, x0), c(y0, y0, y1, y1, y0))
}

# Shapes ------------------------------------------------------------------------

# Pointy-top hexagon (vertices at -30, 30, ..., 270 degrees), closed.
hexagon_ring <- function(cx, cy, radius) {
  a <- (pi / 3) * 0:5 - pi / 6
  m <- cbind(cx + radius * cos(a), cy + radius * sin(a))
  rbind(m, m[1, ])
}

circle_ring <- function(cx, cy, radius, segments = 32) {
  a <- 2 * pi * (0:segments) / segments
  cbind(cx + radius * cos(a), cy + radius * sin(a))
}

rectangle_ring <- function(cx, cy, width, height) {
  hw <- width / 2
  hh <- height / 2
  cbind(
    c(cx - hw, cx + hw, cx + hw, cx - hw, cx - hw),
    c(cy - hh, cy - hh, cy + hh, cy + hh, cy - hh)
  )
}
