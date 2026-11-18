#' Options for post-processing an allocation
#'
#' Heuristic refinement applied by [allocate_grid()] after the MIP is
#' solved. Greedy pairwise swaps always run; simulated annealing and the
#' adjacency fixer are optional. Simulated annealing uses R's random number
#' generator, so call [set.seed()] for reproducible results.
#'
#' Since the MIP solution is optimal for its objective, greedy swaps and
#' annealing only change the result when that objective differs from the one
#' used here (for example when embedding or pairwise terms were used in the
#' MIP), or when the solver stopped early.
#'
#' @param max_swap_iter Maximum number of greedy swap iterations.
#' @param anneal If `TRUE`, run simulated annealing with random swaps and
#'   moves to empty cells.
#' @param anneal_iter Number of annealing steps.
#' @param anneal_temp Initial annealing temperature.
#' @param anneal_cooling Multiplicative cooling factor per step.
#' @param fix_adjacency If `TRUE`, run a local search that reduces the number
#'   of "cut" edges: adjacent features that are not in neighbouring cells.
#'   Requires the `adjacency_graph` argument of [allocate_grid()].
#' @param fix_iter Maximum number of adjacency fixer iterations.
#'
#' @returns A list of class `gridmapper_post_process_options`.
#' @export
#' @examples
#' post_process_options(anneal = TRUE, anneal_iter = 5000)
post_process_options <- function(
  max_swap_iter = 1000,
  anneal = FALSE,
  anneal_iter = 2000,
  anneal_temp = 1,
  anneal_cooling = 0.995,
  fix_adjacency = FALSE,
  fix_iter = 500
) {
  check_count(max_swap_iter)
  check_flag(anneal)
  check_count(anneal_iter)
  check_number(anneal_temp, min = 0)
  check_number(anneal_cooling, min = 0, max = 1)
  check_flag(fix_adjacency)
  check_count(fix_iter)
  structure(
    list(
      max_swap_iter = max_swap_iter,
      anneal = anneal,
      anneal_iter = anneal_iter,
      anneal_temp = anneal_temp,
      anneal_cooling = anneal_cooling,
      fix_adjacency = fix_adjacency,
      fix_iter = fix_iter
    ),
    class = "gridmapper_post_process_options"
  )
}

# Positions (row, col) are 0-based throughout.
post_process_assignments <- function(
  row,
  col,
  norm,
  rows,
  cols,
  opts,
  grid_type = "rect",
  distance = "euclidean",
  compactness_weight = 1,
  adjacency_weight = 0,
  adjacency_diagonal = FALSE,
  spacer_mode = "hard",
  mask_penalty = 1000,
  spacer_keys = character(),
  adjacency_graph = NULL
) {
  n <- length(row)
  soft <- spacer_mode == "soft"

  # Cost of placing points `i` at cells (r, c).
  dist_cost <- function(i, r, c) {
    ctr <- cell_center(r, c, grid_type)
    dx <- norm$x[i] - ctr$x
    dy <- norm$y[i] - ctr$y
    d <- if (distance == "manhattan") abs(dx) + abs(dy) else dx^2 + dy^2
    pen <- if (soft) {
      ifelse(paste0(r, "_", c) %in% spacer_keys, mask_penalty, 0)
    } else {
      0
    }
    d * compactness_weight + pen
  }

  count_adjacent <- function(row, col) {
    occ <- matrix(FALSE, rows + 2, cols + 2)
    occ[cbind(row + 2, col + 2)] <- TRUE
    inner <- occ[2:(rows + 1), 2:(cols + 1)]
    total <- sum(inner & occ[2:(rows + 1), 3:(cols + 2)]) +
      sum(inner & occ[3:(rows + 2), 2:(cols + 1)])
    if (adjacency_diagonal) {
      total <- total +
        sum(inner & occ[3:(rows + 2), 3:(cols + 2)]) +
        sum(inner & occ[3:(rows + 2), 1:cols])
    }
    total
  }

  # Greedy pairwise swaps: first improving swap, in index order. Swapping two
  # occupied cells leaves the occupancy (and so the adjacency term) intact.
  iter <- 0L
  swaps <- 0L
  improved <- TRUE
  idx <- seq_len(n)
  while (improved && iter < opts$max_swap_iter) {
    improved <- FALSE
    iter <- iter + 1L
    cost <- outer(idx, idx, \(i, j) dist_cost(i, row[j], col[j]))
    d0 <- diag(cost)
    delta <- cost + t(cost) - outer(d0, d0, "+")
    cand <- which(upper.tri(delta) & delta < -1e-9, arr.ind = TRUE)
    if (nrow(cand) > 0) {
      best <- cand[order(cand[, 1], cand[, 2])[1], ]
      i <- best[[1]]
      j <- best[[2]]
      tmp <- c(row[i], col[i])
      row[i] <- row[j]
      col[i] <- col[j]
      row[j] <- tmp[1]
      col[j] <- tmp[2]
      swaps <- swaps + 1L
      improved <- TRUE
    }
  }

  if (opts$anneal) {
    all_cells <- grid_cells(rows, cols)
    keys <- paste0(all_cells$row, "_", all_cells$col)
    keep <- !(keys %in% paste0(row, "_", col))
    if (spacer_mode == "hard") {
      keep <- keep & !(keys %in% spacer_keys)
    }
    free <- as.matrix(all_cells[keep, , drop = FALSE])

    total <- function() {
      adj <- if (adjacency_weight != 0) count_adjacent(row, col) else 0
      sum(dist_cost(idx, row, col)) - adj * adjacency_weight
    }
    temp <- opts$anneal_temp
    best_cost <- total()
    best_row <- row
    best_col <- col
    for (s in seq_len(opts$anneal_iter)) {
      if (nrow(free) == 0) {
        i <- sample.int(n, 1)
        j <- sample.int(n, 1)
        if (i == j) {
          next
        }
        prev <- total()
        old <- c(row[i], col[i], row[j], col[j])
        row[c(i, j)] <- old[c(3, 1)]
        col[c(i, j)] <- old[c(4, 2)]
        post <- total()
        delta <- post - prev
        if (delta <= 0 || exp(-delta / temp) > stats::runif(1)) {
          if (post < best_cost) {
            best_cost <- post
            best_row <- row
            best_col <- col
          }
        } else {
          row[c(i, j)] <- old[c(1, 3)]
          col[c(i, j)] <- old[c(2, 4)]
        }
      } else {
        i <- sample.int(n, 1)
        fi <- sample.int(nrow(free), 1)
        prev <- total()
        old <- c(row[i], col[i])
        row[i] <- free[fi, 1]
        col[i] <- free[fi, 2]
        post <- total()
        delta <- post - prev
        if (delta <= 0 || exp(-delta / temp) > stats::runif(1)) {
          free <- rbind(free[-fi, , drop = FALSE], old)
          if (post < best_cost) {
            best_cost <- post
            best_row <- row
            best_col <- col
          }
        } else {
          row[i] <- old[1]
          col[i] <- old[2]
        }
      }
      temp <- temp * opts$anneal_cooling
    }
    row <- best_row
    col <- best_col
  }

  meta <- list(iterations = iter, swaps = swaps)
  if (opts$fix_adjacency) {
    fix <- reduce_cut_edges(
      row,
      col,
      adjacency_graph$edges,
      rows,
      cols,
      max_iter = opts$fix_iter,
      diagonal = adjacency_diagonal,
      spacer_keys = spacer_keys
    )
    row <- fix$row
    col <- fix$col
    meta$adjacency_fix <- fix$meta
  }

  list(row = as.integer(row), col = as.integer(col), meta = meta)
}

# Local search that moves endpoints of cut edges into empty neighbouring
# cells, or swaps them with a neighbour, whenever that reduces the number of
# cut edges.
reduce_cut_edges <- function(
  row,
  col,
  edges,
  rows,
  cols,
  max_iter = 500,
  diagonal = FALSE,
  spacer_keys = character()
) {
  from <- edges$from
  to <- edges$to
  n <- length(row)
  ok <- from >= 1 & from <= n & to >= 1 & to <= n
  from <- from[ok]
  to <- to[ok]

  adjacent <- function(r1, c1, r2, c2) {
    dr <- abs(r1 - r2)
    dc <- abs(c1 - c2)
    if (diagonal) pmax(dr, dc) == 1 else dr + dc == 1
  }
  is_cut <- function(e) {
    !adjacent(row[from[e]], col[from[e]], row[to[e]], col[to[e]])
  }
  count_cuts <- function() sum(is_cut(seq_along(from)))

  occ <- matrix(0L, rows, cols)
  occ[cbind(row + 1, col + 1)] <- seq_len(n)
  is_spacer <- matrix(FALSE, rows, cols)
  if (length(spacer_keys)) {
    sp <- keys_to_spacer_matrix(spacer_keys)
    sp <- sp[sp[, 1] <= rows & sp[, 2] <= cols, , drop = FALSE]
    is_spacer[sp] <- TRUE
  }
  offsets <- neighbour_offsets(diagonal)

  # Change in cut count if the points in `who` moved to (new_r, new_c).
  cut_delta <- function(who, new_r, new_c) {
    e <- which(from %in% who | to %in% who)
    before <- sum(is_cut(e))
    old_r <- row
    old_c <- col
    row[who] <<- new_r
    col[who] <<- new_c
    after <- sum(is_cut(e))
    row <<- old_r
    col <<- old_c
    after - before
  }

  try_move <- function(mover, anchor) {
    for (o in seq_len(nrow(offsets))) {
      nr <- row[anchor] + offsets[o, 1]
      nc <- col[anchor] + offsets[o, 2]
      if (nr < 0 || nc < 0 || nr >= rows || nc >= cols) {
        next
      }
      if (is_spacer[nr + 1, nc + 1] || occ[nr + 1, nc + 1] != 0) {
        next
      }
      if (cut_delta(mover, nr, nc) < 0) {
        occ[row[mover] + 1, col[mover] + 1] <<- 0L
        occ[nr + 1, nc + 1] <<- mover
        row[mover] <<- nr
        col[mover] <<- nc
        return(TRUE)
      }
    }
    FALSE
  }

  try_swap <- function(a, b) {
    for (o in seq_len(nrow(offsets))) {
      nr <- row[b] + offsets[o, 1]
      nc <- col[b] + offsets[o, 2]
      if (nr < 0 || nc < 0 || nr >= rows || nc >= cols) {
        next
      }
      other <- occ[nr + 1, nc + 1]
      if (other == 0 || other == a) {
        next
      }
      pair <- c(a, other)
      if (cut_delta(pair, row[rev(pair)], col[rev(pair)]) < 0) {
        new_r <- row[rev(pair)]
        new_c <- col[rev(pair)]
        row[pair] <<- new_r
        col[pair] <<- new_c
        occ[cbind(new_r + 1, new_c + 1)] <<- pair
        return(TRUE)
      }
    }
    FALSE
  }

  initial <- count_cuts()
  iter <- 0L
  improved <- TRUE
  while (improved && iter < max_iter) {
    improved <- FALSE
    iter <- iter + 1L
    cut <- which(is_cut(seq_along(from)))
    if (length(cut) == 0) {
      break
    }
    for (e in cut) {
      a <- from[e]
      b <- to[e]
      if (try_move(a, b) || try_move(b, a) || try_swap(a, b)) {
        improved <- TRUE
        break
      }
    }
  }

  list(
    row = row,
    col = col,
    meta = list(
      initial_cuts = initial,
      final_cuts = count_cuts(),
      iterations = iter
    )
  )
}
