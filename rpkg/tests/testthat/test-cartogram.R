weighted_london <- function() {
  lon <- london()
  lon$w <- 1 + ((seq_len(nrow(lon)) - 1) * 7) %% 13
  lon
}

test_that("dorling_cartogram() matches the JavaScript library with a seed", {
  ref <- as_num_matrix(js_reference()$london$dorling)
  d <- dorling_cartogram(
    weighted_london(),
    "w",
    seed = 42,
    iterations = 300,
    refinement_iterations = 200
  )
  expect_s3_class(d, "sf")
  expect_equal(
    cbind(d$cartogram_x, d$cartogram_y, d$cartogram_radius),
    ref,
    tolerance = 1e-6
  )
  expect_identical(attr(d, "cartogram")$has_overlaps, FALSE)
})

test_that("demers_cartogram() matches the JavaScript library with a seed", {
  ref <- as_num_matrix(js_reference()$london$demers)
  d <- demers_cartogram(
    weighted_london(),
    "w",
    seed = 7,
    iterations = 300,
    refinement_iterations = 200
  )
  expect_equal(
    cbind(d$cartogram_x, d$cartogram_y, d$cartogram_radius),
    ref,
    tolerance = 1e-6
  )
  expect_identical(attr(d, "cartogram")$shape, "square")
})

test_that("dorling_cartogram() is reproducible with set.seed()", {
  lon <- weighted_london()
  a <- withr::with_seed(
    1,
    dorling_cartogram(lon, "w", shape = "hexagon", iterations = 100)
  )
  b <- withr::with_seed(
    1,
    dorling_cartogram(lon, "w", shape = "hexagon", iterations = 100)
  )
  expect_equal(a$cartogram_x, b$cartogram_x)
  expect_equal(nrow(sf::st_coordinates(a[1, ])), 7)
})

test_that("dorling_cartogram() drops invalid weights", {
  lon <- weighted_london()
  lon$w[1:3] <- c(NA, 0, -1)
  d <- dorling_cartogram(lon, "w", seed = 1, iterations = 50, min_weight = 2)
  expect_equal(nrow(d), sum(lon$w >= 2, na.rm = TRUE))
  expect_snapshot(error = TRUE, {
    dorling_cartogram(lon, "nope")
    dorling_cartogram(lon, "id")
    dorling_cartogram(lon, "w", min_weight = 100)
  })
})
