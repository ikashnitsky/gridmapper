# fix_adjacency requires an adjacency graph

    Code
      allocate_grid(us_states, "lon", "lat", post_process = post_process_options(
        fix_adjacency = TRUE))
    Condition
      Error in `allocate_grid()`:
      ! `fix_adjacency` requires `adjacency_graph`, created by `adjacency_graph()`.
    Code
      allocate_grid(us_states, "lon", "lat", post_process = "yes")
    Condition
      Error in `allocate_grid()`:
      ! `post_process` must be `TRUE`, `FALSE`, `NULL` or created by `post_process_options()`.

