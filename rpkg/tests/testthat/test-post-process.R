test_that("post-processing does not worsen an optimal allocation", {
  ref <- js_reference()$us_states
  gm <- allocate_grid(
    us_states,
    "lon",
    "lat",
    rows = 6,
    cols = 11,
    compactness = 0.6,
    post_process = TRUE
  )
  expect_equal(gm$score, ref$post, tolerance = 1e-8)
  expect_equal(gm$post_process$swaps, 0L)
})

test_that("greedy swaps repair a scrambled allocation", {
  norm <- list(x = c(0, 1, 2), y = c(0, 0, 0))
  res <- post_process_assignments(
    row = c(0L, 0L, 0L),
    col = c(2L, 1L, 0L),
    norm = norm,
    rows = 1,
    cols = 3,
    opts = post_process_options()
  )
  expect_equal(res$col, c(0L, 1L, 2L))
  # First-improvement search: (2, 3), then (1, 2), then (2, 3) again
  expect_equal(res$meta$swaps, 3L)
})

test_that("simulated annealing keeps cells distinct", {
  set.seed(42)
  gm <- allocate_grid(
    us_states,
    "lon",
    "lat",
    rows = 6,
    cols = 11,
    post_process = post_process_options(anneal = TRUE, anneal_iter = 300)
  )
  expect_length(unique(paste(gm$assignments$row, gm$assignments$col)), 50)
})

test_that("the adjacency fixer reduces cut edges", {
  edges <- data.frame(from = c(1, 2), to = c(2, 3))
  res <- reduce_cut_edges(
    row = c(0, 0, 0),
    col = c(0, 1, 4),
    edges = edges,
    rows = 1,
    cols = 5
  )
  expect_equal(res$meta$initial_cuts, 1)
  expect_equal(res$meta$final_cuts, 0)
  expect_equal(res$col, c(0, 1, 2))

  lon <- london()
  gm <- allocate_grid(
    lon,
    rows = 7,
    cols = 9,
    compactness = 0.6,
    adjacency_graph = adjacency_graph(lon),
    post_process = post_process_options(fix_adjacency = TRUE)
  )
  fix <- gm$post_process$adjacency_fix
  expect_lte(fix$final_cuts, fix$initial_cuts)
  expect_length(unique(paste(gm$assignments$row, gm$assignments$col)), 33)
})

test_that("fix_adjacency requires an adjacency graph", {
  expect_snapshot(error = TRUE, {
    allocate_grid(
      us_states,
      "lon",
      "lat",
      post_process = post_process_options(fix_adjacency = TRUE)
    )
    allocate_grid(us_states, "lon", "lat", post_process = "yes")
  })
})
