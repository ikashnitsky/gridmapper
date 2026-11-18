# allocate_grid() validates its input

    Code
      allocate_grid(us_states, "lon", "lat", rows = 5, cols = 5)
    Condition
      Error in `allocate_grid()`:
      ! Can't allocate 50 features to 25 available cells.
      i The 5 x 5 grid has 0 spacers.
      i Increase `rows` or `cols`, or reduce the spacers.
    Code
      allocate_grid(us_states, "long", "lat")
    Condition
      Error in `allocate_grid()`:
      ! Column long not found in `data`.
    Code
      allocate_grid(us_states, "lon", "lat", compactness = 2)
    Condition
      Error in `allocate_grid()`:
      ! `compactness` must be a single number between 0 and 1, not 2.
    Code
      allocate_grid(us_states, "lon", "lat", grid_type = "tri")
    Condition
      Error in `allocate_grid()`:
      ! `grid_type` must be one of "rect" or "hex", not "tri".
    Code
      allocate_grid(us_states, "lon", "lat", pairwise_adjacency = pairwise_adjacency_options())
    Condition
      Error in `allocate_grid()`:
      ! `pairwise_adjacency` requires `adjacency_graph`, created by `adjacency_graph()`.
    Code
      allocate_grid("a")
    Condition
      Error in `allocate_grid()`:
      ! `data` must be a data frame, matrix or <sf> object, not a string.

# allocate_grid_simple() returns cells in point order

    Code
      allocate_grid_simple(pts, 5, 5)
    Condition
      Error in `allocate_grid_simple()`:
      ! Can't allocate 50 points to a grid with only 25 cells.

# gridmap methods work

    Code
      print(gm)
    Message
      <gridmap> 50 features on a 5 x 10 rect grid
      Score: 378.633 (Optimal)
      Post-processing: 0 swaps in 1 iteration

