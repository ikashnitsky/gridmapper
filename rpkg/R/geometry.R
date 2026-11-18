# Internal geometry representation ---------------------------------------
#
# Features are handled as GeoJSON-like lists, `list(type = , coords = )`:
# * "Point": numeric(2)
# * "Polygon": list of ring matrices (outer ring first, then holes)
# * "MultiPolygon": list of polygons
# All computations are planar on raw coordinates, exactly as in the
# JavaScript library (no spherical geometry, CRS is ignored).

as_geometry_list <- function(x, call = rlang::caller_env()) {
  if (inherits(x, "sf")) {
    x <- sf::st_geometry(x)
  } else if (inherits(x, "sfg")) {
    x <- sf::st_sfc(x)
  }
  if (!inherits(x, "sfc")) {
    cli::cli_abort(
      "{.arg x} must be an {.cls sf} or {.cls sfc} object, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  lapply(x, sfg_to_geom)
}

sfg_to_geom <- function(g) {
  xy <- function(m) unname(m[, 1:2, drop = FALSE])
  if (inherits(g, "POINT")) {
    list(type = "Point", coords = unname(as.numeric(g)[1:2]))
  } else if (inherits(g, "POLYGON")) {
    list(type = "Polygon", coords = lapply(unclass(g), xy))
  } else if (inherits(g, "MULTIPOLYGON")) {
    list(
      type = "MultiPolygon",
      coords = lapply(unclass(g), \(p) lapply(p, xy))
    )
  } else {
    list(type = class(g)[2], coords = NULL)
  }
}

# All vertices of a geometry as a 2-column matrix.
geom_vertices <- function(geom) {
  switch(
    geom$type,
    Point = matrix(geom$coords, ncol = 2),
    Polygon = do.call(rbind, geom$coords),
    MultiPolygon = do.call(rbind, unlist(geom$coords, recursive = FALSE)),
    matrix(numeric(), ncol = 2)
  )
}

geom_bbox <- function(geom) {
  v <- geom_vertices(geom)
  if (nrow(v) == 0) {
    min_x <- Inf
    max_x <- -Inf
    min_y <- Inf
    max_y <- -Inf
  } else {
    min_x <- min(v[, 1])
    max_x <- max(v[, 1])
    min_y <- min(v[, 2])
    max_y <- max(v[, 2])
  }
  c(
    min_x = min_x,
    max_x = max_x,
    min_y = min_y,
    max_y = max_y,
    width = max_x - min_x,
    height = max_y - min_y,
    area = (max_x - min_x) * (max_y - min_y)
  )
}

geoms_bbox <- function(geoms) {
  do.call(rbind, lapply(geoms, geom_bbox))
}

# Shoelace area including the (i, i+1 mod n) closing segment, unsigned.
ring_area <- function(ring) {
  n <- nrow(ring)
  if (n < 3) {
    return(0)
  }
  nxt <- c(2:n, 1)
  abs(sum(ring[, 1] * ring[nxt, 2] - ring[nxt, 1] * ring[, 2])) * 0.5
}

ring_signed_area <- function(ring) {
  n <- nrow(ring)
  prv <- c(n, 1:(n - 1))
  sum(ring[prv, 1] * ring[, 2] - ring[, 1] * ring[prv, 2]) / 2
}

ring_length <- function(ring) {
  n <- nrow(ring)
  if (n < 2) {
    return(0)
  }
  nxt <- c(2:n, 1)
  sum(sqrt((ring[nxt, 1] - ring[, 1])^2 + (ring[nxt, 2] - ring[, 2])^2))
}

polygon_stats <- function(rings) {
  if (length(rings) == 0) {
    return(c(area = 0, perimeter = 0, circularity = 0, hole_area = 0))
  }
  areas <- vapply(rings, ring_area, numeric(1))
  perimeter <- sum(vapply(rings, ring_length, numeric(1)))
  hole_area <- sum(areas[-1])
  area <- max(0, areas[1] - hole_area)
  circularity <- if (perimeter > 0) 4 * pi * area / perimeter^2 else 0
  c(
    area = area,
    perimeter = perimeter,
    circularity = max(0, min(1, circularity)),
    hole_area = max(0, hole_area)
  )
}

geom_stats <- function(geom) {
  if (geom$type == "Polygon") {
    return(polygon_stats(geom$coords))
  }
  if (geom$type == "MultiPolygon") {
    s <- vapply(geom$coords, polygon_stats, numeric(4))
    area <- sum(s["area", ])
    perimeter <- sum(s["perimeter", ])
    circularity <- if (perimeter > 0) 4 * pi * area / perimeter^2 else 0
    return(c(
      area = area,
      perimeter = perimeter,
      circularity = circularity,
      hole_area = sum(s["hole_area", ])
    ))
  }
  c(area = 0, perimeter = 0, circularity = 0, hole_area = 0)
}

# Pole of inaccessibility --------------------------------------------------

#' Interior label points of polygons
#'
#' Computes an interior point for each feature using the pole of
#' inaccessibility (the "polylabel" algorithm). For multipolygons the largest
#' part (by outer ring area) is used; points are returned unchanged. These
#' are the coordinates [allocate_grid()] and [estimate_parameters()] use for
#' polygon inputs.
#'
#' Computation is planar on the raw coordinates, so for longitude/latitude
#' data the result is in degrees.
#'
#' @param x An [sf][sf::st_sf] or [sfc][sf::st_sfc] object with `POINT`,
#'   `POLYGON` or `MULTIPOLYGON` geometries.
#' @param precision Precision of the search, in coordinate units. The
#'   default of `1` matches the JavaScript library; for longitude/latitude
#'   data that is coarse, and smaller values give a more central point.
#'
#' @returns A numeric matrix with columns `x` and `y` and one row per feature.
#'   Unsupported geometry types give `NA`.
#' @export
#' @examples
#' london <- sf::st_read(
#'   system.file("extdata", "london.geojson", package = "gridmapper"),
#'   quiet = TRUE
#' )
#' head(polygon_centroid(london))
polygon_centroid <- function(x, precision = 1) {
  check_number(precision, min = 0)
  geoms <- as_geometry_list(x)
  out <- t(vapply(geoms, geom_centroid, numeric(2), precision = precision))
  colnames(out) <- c("x", "y")
  out
}

geom_centroid <- function(geom, precision = 1) {
  switch(
    geom$type,
    Point = geom$coords,
    Polygon = polylabel(geom$coords, precision),
    MultiPolygon = {
      areas <- vapply(
        geom$coords,
        \(p) abs(ring_signed_area(p[[1]])),
        numeric(1)
      )
      polylabel(geom$coords[[which.max(areas)]], precision)
    },
    c(NA_real_, NA_real_)
  )
}

# Port of mapbox/polylabel 1.x (ISC licence).
polylabel <- function(polygon, precision = 1) {
  if (precision == 0) {
    precision <- 1
  }
  outer <- polygon[[1]]
  min_x <- min(outer[, 1])
  max_x <- max(outer[, 1])
  min_y <- min(outer[, 2])
  max_y <- max(outer[, 2])
  width <- max_x - min_x
  height <- max_y - min_y
  cell_size <- min(width, height)
  h <- cell_size / 2

  if (cell_size == 0) {
    return(c(min_x, min_y))
  }

  segs <- polygon_segments(polygon)
  new_cell <- function(x, y, h) {
    d <- point_polygon_dist(x, y, segs)
    c(x = x, y = y, h = h, d = d, max = d + h * sqrt(2))
  }

  queue <- list()
  x <- min_x
  while (x < max_x) {
    y <- min_y
    while (y < max_y) {
      queue[[length(queue) + 1]] <- new_cell(x + h, y + h, h)
      y <- y + cell_size
    }
    x <- x + cell_size
  }
  queue <- do.call(rbind, queue)

  best <- centroid_cell(outer, new_cell)
  bbox_cell <- new_cell(min_x + width / 2, min_y + height / 2, 0)
  if (bbox_cell[["d"]] > best[["d"]]) {
    best <- bbox_cell
  }

  while (NROW(queue) > 0) {
    i <- which.max(queue[, "max"])
    cell <- queue[i, ]
    queue <- queue[-i, , drop = FALSE]

    if (cell[["d"]] > best[["d"]]) {
      best <- cell
    }
    if (cell[["max"]] - best[["d"]] <= precision) {
      next
    }
    h <- cell[["h"]] / 2
    queue <- rbind(
      queue,
      new_cell(cell[["x"]] - h, cell[["y"]] - h, h),
      new_cell(cell[["x"]] + h, cell[["y"]] - h, h),
      new_cell(cell[["x"]] - h, cell[["y"]] + h, h),
      new_cell(cell[["x"]] + h, cell[["y"]] + h, h)
    )
  }

  unname(best[c("x", "y")])
}

# Segments (a = ring[i], b = ring[i - 1], wrapping) of all rings.
polygon_segments <- function(polygon) {
  do.call(
    rbind,
    lapply(polygon, function(ring) {
      n <- nrow(ring)
      prv <- c(n, seq_len(n - 1))
      cbind(
        ax = ring[, 1],
        ay = ring[, 2],
        bx = ring[prv, 1],
        by = ring[prv, 2]
      )
    })
  )
}

point_polygon_dist <- function(x, y, segs) {
  ax <- segs[, "ax"]
  ay <- segs[, "ay"]
  bx <- segs[, "bx"]
  by <- segs[, "by"]

  crosses <- ((ay > y) != (by > y)) &
    (x < (bx - ax) * (y - ay) / (by - ay) + ax)
  inside <- sum(crosses, na.rm = TRUE) %% 2 == 1

  dx <- bx - ax
  dy <- by - ay
  len2 <- dx * dx + dy * dy
  t <- ifelse(len2 > 0, ((x - ax) * dx + (y - ay) * dy) / len2, 0)
  px <- ifelse(t > 1, bx, ifelse(t > 0, ax + dx * t, ax))
  py <- ifelse(t > 1, by, ifelse(t > 0, ay + dy * t, ay))
  min_dist2 <- min((x - px)^2 + (y - py)^2)

  if (min_dist2 == 0) 0 else (if (inside) 1 else -1) * sqrt(min_dist2)
}

centroid_cell <- function(ring, new_cell) {
  n <- nrow(ring)
  prv <- c(n, seq_len(n - 1))
  f <- ring[, 1] * ring[prv, 2] - ring[prv, 1] * ring[, 2]
  area <- sum(f * 3)
  if (area == 0) {
    return(new_cell(ring[1, 1], ring[1, 2], 0))
  }
  new_cell(
    sum((ring[, 1] + ring[prv, 1]) * f) / area,
    sum((ring[, 2] + ring[prv, 2]) * f) / area,
    0
  )
}

# Convex hull (Andrew's monotone chain). Returns hull vertices as a matrix.
convex_hull <- function(pts) {
  if (nrow(pts) <= 1) {
    return(pts)
  }
  pts <- pts[order(pts[, 1], pts[, 2]), , drop = FALSE]
  cross <- function(o, a, b) {
    (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
  }
  half <- function(idx) {
    hull <- integer()
    for (i in idx) {
      while (
        length(hull) >= 2 &&
          cross(
            pts[hull[length(hull) - 1], ],
            pts[hull[length(hull)], ],
            pts[i, ]
          ) <=
            0
      ) {
        hull <- hull[-length(hull)]
      }
      hull <- c(hull, i)
    }
    hull[-length(hull)]
  }
  n <- nrow(pts)
  pts[c(half(seq_len(n)), half(rev(seq_len(n)))), , drop = FALSE]
}

convex_hull_area <- function(pts) {
  if (nrow(pts) < 3) {
    return(0)
  }
  hull <- convex_hull(pts)
  if (nrow(hull) < 3) {
    return(0)
  }
  ring_area(hull)
}

# Ray-casting point-in-polygon test for a single ring matrix.
point_in_polygon <- function(x, y, ring) {
  if (NROW(ring) == 0) {
    return(FALSE)
  }
  n <- nrow(ring)
  prv <- c(n, seq_len(n - 1))
  xi <- ring[, 1]
  yi <- ring[, 2]
  xj <- ring[prv, 1]
  yj <- ring[prv, 2]
  hits <- ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
  sum(hits, na.rm = TRUE) %% 2 == 1
}
