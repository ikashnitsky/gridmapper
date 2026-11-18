test_that("grid_to_sf() builds one tile per feature", {
  lon <- london()
  gm <- allocate_grid(lon, rows = 7, cols = 9)
  tiles <- grid_to_sf(gm)
  expect_s3_class(tiles, "sf")
  expect_equal(nrow(tiles), 33)
  expect_equal(sf::st_crs(tiles), sf::st_crs(lon))
  expect_contains(names(tiles), c("id", "row", "col"))
  expect_equal(nrow(sf::st_coordinates(tiles[1, ])), 5)

  # Tiles tessellate the bounding box of the input
  bb <- sf::st_bbox(tiles)
  expect_gte(bb[["xmin"]], gm$bounds$min_x - 1e-12)
  expect_lte(bb[["xmax"]], gm$bounds$max_x + 1e-12)
  planar <- sf::st_set_crs(tiles, NA)
  expect_false(any(sf::st_overlaps(planar, sparse = FALSE)))
})

test_that("grid_to_sf() draws hexagons", {
  gm <- allocate_grid(us_states, "lon", "lat", grid_type = "hex")
  tiles <- grid_to_sf(gm)
  expect_equal(nrow(sf::st_coordinates(tiles[1, ])), 7)
  expect_true(is.na(sf::st_crs(tiles)))
  expect_snapshot(error = TRUE, grid_to_sf(us_states))
})
