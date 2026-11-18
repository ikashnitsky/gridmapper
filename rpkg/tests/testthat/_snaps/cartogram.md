# dorling_cartogram() drops invalid weights

    Code
      dorling_cartogram(lon, "nope")
    Condition
      Error in `dorling_cartogram()`:
      ! `weight` must be the name of a column of `x`.
    Code
      dorling_cartogram(lon, "id")
    Condition
      Error in `dorling_cartogram()`:
      ! Column id must be numeric.
    Code
      dorling_cartogram(lon, "w", min_weight = 100)
    Condition
      Error in `dorling_cartogram()`:
      ! No features with valid positive weights found in column w.

