test_that("polygon_centroid() matches JavaScript polylabel", {
  ref <- js_reference()$london
  expect_equal(
    unname(polygon_centroid(london())),
    as_num_matrix(ref$centroids),
    tolerance = 1e-9
  )
})

test_that("polygon_centroid() finds an interior point", {
  # U-shaped polygon whose area centroid lies outside it
  u <- sf::st_polygon(list(cbind(
    c(0, 3, 3, 2, 2, 1, 1, 0, 0),
    c(0, 0, 3, 3, 1, 1, 3, 3, 0)
  )))
  pt <- polygon_centroid(sf::st_sfc(u), precision = 0.01)
  inside <- sf::st_intersects(
    sf::st_sfc(sf::st_point(pt[1, ])),
    sf::st_sfc(u),
    sparse = FALSE
  )
  expect_identical(inside[1, 1], TRUE)
  expect_equal(
    unname(polygon_centroid(sf::st_sfc(sf::st_point(c(2, 3))))),
    matrix(c(2, 3), 1)
  )
})

test_that("polygon_centroid() uses the largest part of multipolygons", {
  small <- list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))
  big <- list(cbind(c(10, 20, 20, 10, 10), c(0, 0, 10, 10, 0)))
  mp <- sf::st_sfc(sf::st_multipolygon(list(small, big)))
  expect_equal(unname(polygon_centroid(mp)[1, ]), c(15, 5))
})

test_that("polygon_centroid() rejects non-spatial input", {
  expect_snapshot(error = TRUE, polygon_centroid(data.frame(x = 1)))
})
