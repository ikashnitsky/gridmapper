test_that("grid_dimensions() covers n cells with the requested shape", {
  expect_equal(grid_dimensions(50, 0.5), c(rows = 5L, cols = 10L))
  expect_equal(
    grid_dimensions(33, 1, target_cells = 63, min_rows = 7),
    c(rows = 8L, cols = 8L)
  )
  dims <- grid_dimensions(7, aspect_ratio = 100)
  expect_gte(prod(dims), 7)
})

test_that("grid_dimensions() rounds halves up like JavaScript", {
  # sqrt(cells / aspect) = 2.5 exactly
  expect_equal(grid_dimensions(25, aspect_ratio = 4)[["cols"]], 3L)
})

test_that("pca_angle() and pca_rotate() align points with the x axis", {
  x <- c(0, 1, 2, 3)
  y <- c(0, 1, 2, 3)
  angle <- pca_angle(x, y)
  expect_equal(angle, pi / 4)
  rot <- pca_rotate(x, y, -angle)
  expect_equal(unname(rot[, 2]), rep(1.5, 4))
  expect_equal(pca_angle(1, 1), 0)
})

test_that("compute_spacers() returns cells outside the mask", {
  x <- c(0, 1, 0, 1, 0.5)
  y <- c(0, 0, 1, 1, 0.5)
  # Triangle x + y <= 1.2, away from all cell positions
  tri <- cbind(c(-0.1, 1.3, -0.1, -0.1), c(-0.1, -0.1, 1.3, -0.1))
  sp <- compute_spacers(x, y, rows = 3, cols = 3, mask = tri)
  expect_equal(sp, cbind(row = c(2L, 3L, 3L), col = c(3L, 2L, 3L)))
  hex <- compute_spacers(
    x,
    y,
    rows = 3,
    cols = 3,
    mask = tri,
    grid_type = "hex"
  )
  expect_equal(hex, cbind(row = c(2L, 3L, 3L), col = c(3L, 2L, 3L)))
})
