# `points_normalised()` centres y on the column midpoint, shifting points upward on non-square grids

## Description

In `points_normalised()`, the y coordinates are rescaled around `cc` (the
column midpoint) instead of `cr` (the row midpoint). `cr` is computed but
never used:

```r
cr <- (n_rows + 1) / 2
cc <- (n_cols + 1) / 2
...
y = rescale(y, to = c(cc - r_height / 2, cc + r_height / 2), from = c(min_y, max_y)),
x = rescale(x, to = c(cc - r_width / 2, cc + r_width / 2), from = c(min_x, max_x))
```

This is harmless on square grids, where `cr == cc`. When `n_col > n_row`,
the y targets are moved up, past the top row. When `n_col < n_row`, they
are moved down. The solver then pulls points towards the top (or bottom)
edge of the grid, and the layout no longer reflects the geographic
arrangement.

## Reproducible example

```r
library(gridmappr)
library(tibble)

grd <- grid_locations(n_row = 5, n_col = 10)
pts <- tibble(x = c(0, 1), y = c(0, 1))

points_normalised(pts, grd, compactness = 0.5)
#> # A tibble: 2 × 2
#>         x     y
#>     <dbl> <dbl>
#> 1  0.0200  2.51
#> 2 11.0     8.49
```

At `compactness = 0.5`, the targets should span the grid. The x range
(0.02 to 10.98) does span columns 1 to 10. The y range (2.51 to 8.49) is
centred on 5.5, but the grid only has rows 1 to 5. Most of the y range
lies above the top row.

## Expected behaviour

Rescaled around `cr = 3`, the y targets would run from about 0.01 to
5.99, which covers rows 1 to 5 the same way x covers the columns.

## Suggested fix

```diff
-      y = rescale(y, to = c(cc - r_height / 2, cc + r_height / 2), from = c(min_y, max_y)),
+      y = rescale(y, to = c(cr - r_height / 2, cr + r_height / 2), from = c(min_y, max_y)),
```

This changes allocations on non-square grids, so saved layouts and
vignette figures that use rectangular grids may need regenerating.
