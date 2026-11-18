test_that("adjacency_graph() matches the JavaScript library", {
  ref <- as_num_matrix(js_reference()$london$adj_edges)
  g <- adjacency_graph(london())
  expect_s3_class(g, "gridmapper_adjacency")
  expect_equal(g$edges$from, ref[, 1] + 1)
  expect_equal(g$edges$to, ref[, 2] + 1)
  expect_equal(g$edges$weight, ref[, 3])
  expect_equal(g$nodes$id[1], "Kingston upon Thames")
})

test_that("adjacency_graph() links touching boxes only", {
  sq <- function(x0) {
    sf::st_polygon(list(cbind(c(x0, x0 + 1, x0 + 1, x0, x0), c(0, 0, 1, 1, 0))))
  }
  g <- adjacency_graph(sf::st_sfc(sq(0), sq(1), sq(5)))
  expect_equal(g$edges[c("from", "to")], data.frame(from = 1L, to = 2L))
  expect_equal(g$edges$weight, 1)
})

test_that("embed_graph() places linked nodes closer than unlinked ones", {
  sq <- function(x0) {
    sf::st_polygon(list(cbind(c(x0, x0 + 1, x0 + 1, x0, x0), c(0, 0, 1, 1, 0))))
  }
  g <- adjacency_graph(sf::st_sfc(lapply(c(0, 1, 2, 3, 4, 5), sq)))
  set.seed(3)
  pos <- embed_graph(g, iterations = 300)
  d <- as.matrix(dist(pos))
  expect_lt(d[1, 2], d[1, 6])
  expect_all_true(as.vector(pos >= 0 & pos <= 1))
})
