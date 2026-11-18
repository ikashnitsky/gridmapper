test_that("allocate_grid() matches the JavaScript library on London", {
  ref <- js_reference()$london
  xy <- polygon_centroid(london())
  p <- ref$params

  cases <- list(
    est = list(rows = p$rows, cols = p$cols, compactness = p$compactness),
    hex = list(
      rows = p$rows,
      cols = p$cols,
      compactness = 0.6,
      grid_type = "hex"
    ),
    auto = list(compactness = 0.5),
    adj = list(
      rows = p$rows,
      cols = p$cols,
      compactness = 0.5,
      adjacency_weight = 0.1
    )
  )
  for (case in names(cases)) {
    gm <- do.call(allocate_grid, c(list(xy), cases[[case]]))
    js <- ref$res[[case]]
    expect_equal(gm$score, js$score, tolerance = 1e-8, info = case)
    expect_equal(c(gm$rows, gm$cols), c(js$rows, js$cols), info = case)
    expect_equal(gm$assignments$col, unlist(js$gx) + 1, info = case)
    expect_equal(gm$assignments$row, unlist(js$gy) + 1, info = case)
  }
})

test_that("allocate_grid() matches the JavaScript objective with ties", {
  ref <- js_reference()$london
  gm <- allocate_grid(
    polygon_centroid(london()),
    rows = 7,
    cols = 9,
    distance = "manhattan"
  )
  expect_equal(gm$score, ref$res$manh$score, tolerance = 1e-8)
})

test_that("allocate_grid() matches the JavaScript library on US states", {
  ref <- js_reference()$us_states
  auto <- allocate_grid(us_states, "lon", "lat")
  pca <- allocate_grid(
    us_states,
    "lon",
    "lat",
    rows = 6,
    cols = 11,
    compactness = 0.6,
    rotate_pca = TRUE
  )
  expect_equal(auto$score, ref$auto, tolerance = 1e-8)
  expect_equal(pca$score, ref$pca, tolerance = 1e-8)
})

test_that("allocate_grid() assigns every feature to a distinct cell", {
  gm <- allocate_grid(us_states, "lon", "lat", grid_type = "hex")
  cells <- paste(gm$assignments$row, gm$assignments$col)
  expect_length(unique(cells), nrow(us_states))
  expect_all_true(gm$assignments$row %in% seq_len(gm$rows))
  expect_all_true(gm$assignments$col %in% seq_len(gm$cols))
  expect_identical(gm$feasible, TRUE)
})

test_that("allocate_grid() accepts sf, matrix and data frame input", {
  lon <- london()
  xy <- polygon_centroid(lon)
  df <- data.frame(name = lon$id, x = xy[, 1], y = xy[, 2])
  a <- allocate_grid(lon, rows = 7, cols = 9)
  b <- allocate_grid(xy, rows = 7, cols = 9)
  c <- allocate_grid(df, rows = 7, cols = 9)
  expect_equal(a$assignments, b$assignments)
  expect_equal(b$assignments, c$assignments)
})

test_that("hard spacers stay empty and soft spacers are penalised", {
  spacers <- cbind(row = 1:3, col = 1:3)
  hard <- allocate_grid(
    us_states,
    "lon",
    "lat",
    rows = 6,
    cols = 10,
    spacers = spacers
  )
  used <- paste(hard$assignments$row, hard$assignments$col)
  expect_disjoint(used, paste(spacers[, 1], spacers[, 2]))
  expect_equal(hard$spacers, spacers)

  # 50 states on 50 cells forces one feature into the soft spacer
  soft <- allocate_grid(
    us_states,
    "lon",
    "lat",
    rows = 5,
    cols = 10,
    spacers = cbind(1, 1),
    spacer_mode = "soft",
    mask_penalty = 7
  )
  plain <- allocate_grid(us_states, "lon", "lat", rows = 5, cols = 10)
  expect_equal(soft$score, plain$score + 7, tolerance = 1e-6)
})

test_that("allocate_grid() derives spacers from a mask", {
  # Excludes the top three rows of the grid (north of 50 degrees)
  ring <- cbind(c(-160, -65, -65, -160, -160), c(18, 18, 50, 50, 18))
  gm <- allocate_grid(us_states, "lon", "lat", rows = 8, cols = 12, mask = ring)
  expect_equal(nrow(gm$spacers), 36)
  expect_all_true(gm$spacers[, "row"] >= 6)
  used <- paste(gm$assignments$row, gm$assignments$col)
  expect_disjoint(used, paste(gm$spacers[, 1], gm$spacers[, 2]))
})

test_that("allocate_grid() uses embedding targets", {
  xy <- polygon_centroid(london())
  # Targets equal to the input positions don't change the allocation...
  same <- allocate_grid(
    xy,
    rows = 7,
    cols = 9,
    embedding_targets = xy,
    embedding_weight = 1
  )
  base <- allocate_grid(xy, rows = 7, cols = 9)
  expect_equal(same$assignments, base$assignments)
  # ...while mirrored targets pull the layout apart.
  mirrored <- cbind(-xy[, 1], xy[, 2])
  pulled <- allocate_grid(
    xy,
    rows = 7,
    cols = 9,
    embedding_targets = mirrored,
    embedding_weight = 5
  )
  expect_false(identical(pulled$assignments$col, base$assignments$col))
})

test_that("pairwise adjacency adds auxiliary variables", {
  lon <- london()
  graph <- adjacency_graph(lon)
  gm <- allocate_grid(
    lon,
    rows = 7,
    cols = 9,
    adjacency_graph = graph,
    pairwise_adjacency = pairwise_adjacency_options(
      edge_limit = 50,
      nearest_cells = 6,
      weight = 0.1
    )
  )
  expect_gt(gm$pairwise_vars, 0)
  expect_identical(gm$feasible, TRUE)
})

test_that("allocate_grid() validates its input", {
  expect_snapshot(error = TRUE, {
    allocate_grid(us_states, "lon", "lat", rows = 5, cols = 5)
    allocate_grid(us_states, "long", "lat")
    allocate_grid(us_states, "lon", "lat", compactness = 2)
    allocate_grid(us_states, "lon", "lat", grid_type = "tri")
    allocate_grid(
      us_states,
      "lon",
      "lat",
      pairwise_adjacency = pairwise_adjacency_options()
    )
    allocate_grid("a")
  })
})

test_that("allocate_grid_simple() returns cells in point order", {
  pts <- as.matrix(us_states[c("lon", "lat")])
  res <- allocate_grid_simple(
    pts,
    6,
    10,
    compactness = 0.6,
    spacers = cbind(1, 1)
  )
  expect_equal(dim(res$cells), c(50, 2))
  expect_length(unique(paste(res$cells[, 1], res$cells[, 2])), 50)
  expect_disjoint(paste(res$cells[, 1], res$cells[, 2]), "1 1")

  # Equivalent to allocate_grid() on a rectangular grid
  gm <- allocate_grid(
    pts,
    rows = 6,
    cols = 10,
    compactness = 0.6,
    spacers = cbind(1, 1)
  )
  expect_equal(
    unname(res$cells),
    unname(as.matrix(gm$assignments[c("row", "col")]))
  )
  expect_snapshot(error = TRUE, allocate_grid_simple(pts, 5, 5))
})

test_that("gridmap methods work", {
  gm <- allocate_grid(us_states, "lon", "lat", post_process = TRUE)
  expect_snapshot(print(gm))
  df <- as.data.frame(gm)
  expect_named(df, c(names(us_states), "row", "col"))
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot(gm, labels = "abbr"), gm)
})
