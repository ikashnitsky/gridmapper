# polygon_centroid() rejects non-spatial input

    Code
      polygon_centroid(data.frame(x = 1))
    Condition
      Error in `polygon_centroid()`:
      ! `x` must be an <sf> or <sfc> object, not a data frame.

