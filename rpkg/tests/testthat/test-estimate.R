test_that("estimate_parameters() matches the JavaScript library", {
  ref <- js_reference()$london$params
  p <- estimate_parameters(london())
  expect_equal(p$rows, ref$rows)
  expect_equal(p$cols, ref$cols)
  expect_equal(p$compactness, ref$compactness)
  expect_identical(p$rotate_pca, ref$rotateByPCA)
})

test_that("estimate_parameters() matches JavaScript on a collinear strip", {
  # A diagonal strip of square "regions"
  sq <- function(i) {
    sf::st_polygon(list(cbind(
      c(i, i + 1, i + 1, i, i),
      c(i, i, i + 1, i + 1, i) * 0.3
    )))
  }
  strip <- sf::st_sf(id = 1:20, geometry = sf::st_sfc(lapply(1:20, sq)))
  p <- estimate_parameters(strip)
  # Value from the JavaScript library. Perfectly collinear centroids are a
  # degenerate case of its heuristics, hence the tall grid.
  expect_equal(
    p,
    list(rows = 20L, cols = 15L, compactness = 0.438, rotate_pca = TRUE)
  )
})

test_that("estimate_parameters() can return embedding targets", {
  set.seed(1)
  p <- estimate_parameters(
    london(),
    embedding = TRUE,
    embedding_options = list(iterations = 50)
  )
  expect_s3_class(p$adjacency_graph, "gridmapper_adjacency")
  expect_equal(dim(p$embedding_targets), c(33, 2))
  gm <- allocate_grid(
    london(),
    rows = p$rows,
    cols = p$cols,
    embedding_targets = p$embedding_targets,
    embedding_weight = 0.2
  )
  expect_identical(gm$feasible, TRUE)
})

test_that("estimate_parameters_from_data() matches the JavaScript library", {
  ref <- js_reference()
  p <- estimate_parameters_from_data(us_states, "lon", "lat")
  expect_equal(p$rows, ref$us_states$from_data$rows)
  expect_equal(p$cols, ref$us_states$from_data$cols)
  expect_equal(p$compactness, ref$us_states$from_data$compactness)
  expect_identical(p$rotate_pca, ref$us_states$from_data$rotateByPCA)

  xy <- polygon_centroid(london())
  p <- estimate_parameters_from_data(xy)
  expect_equal(p$compactness, ref$london$from_data$compactness)
})
