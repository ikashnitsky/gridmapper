# Builds and solves the assignment MIP with HiGHS.
#
# Variables, in order:
# 1. x[p, k]: point p assigned to cell k (column-major: (k - 1) * n + p)
# 2. occ[k]: cell k occupied (only when adjacency_weight != 0)
# 3. pij: both cells of an adjacent cell pair occupied
# 4. y: feature pair (a, b) of an adjacency edge placed in adjacent cells
#
# All costs are summed per variable. (The JavaScript library passes repeated
# objective terms to glpk.js, which keeps only the last one, so its embedding
# and soft-spacer penalties replace the distance cost instead of adding to
# it.)
solve_allocation <- function(
  norm_x,
  norm_y,
  cells,
  soft_spacer,
  grid_type = "rect",
  distance = "euclidean",
  compactness_weight = 1,
  mask_penalty = 1000,
  adjacency_weight = 0,
  adjacency_diagonal = FALSE,
  embedding = NULL,
  embedding_weight = 0,
  pairwise = NULL,
  solver_options = list(),
  call = rlang::caller_env()
) {
  n <- length(norm_x)
  k <- nrow(cells)
  centers <- cell_center(cells$row, cells$col, grid_type)

  dx <- outer(norm_x, centers$x, "-")
  dy <- outer(norm_y, centers$y, "-")
  cost <- if (distance == "manhattan") abs(dx) + abs(dy) else dx^2 + dy^2
  cost <- cost * compactness_weight

  if (!is.null(embedding) && embedding_weight != 0) {
    ex <- outer(embedding$x, centers$x, "-")
    ey <- outer(embedding$y, centers$y, "-")
    emb <- sqrt(ex^2 + ey^2) * embedding_weight
    emb[is.na(emb)] <- 0
    cost <- cost + emb
  }
  if (any(soft_spacer)) {
    cost[, soft_spacer] <- cost[, soft_spacer] + mask_penalty
  }

  x_index <- matrix(seq_len(n * k), nrow = n, ncol = k)
  obj <- as.vector(cost)
  n_vars <- n * k

  # Constraint triplets
  ci <- list()
  cj <- list()
  cv <- list()
  lhs <- list()
  rhs <- list()
  n_rows <- 0L
  add_rows <- function(i, j, v, lo, hi) {
    ci[[length(ci) + 1]] <<- i + n_rows
    cj[[length(cj) + 1]] <<- j
    cv[[length(cv) + 1]] <<- v
    lhs[[length(lhs) + 1]] <<- lo
    rhs[[length(rhs) + 1]] <<- hi
    n_rows <<- n_rows + length(lo)
  }

  # Each point in exactly one cell
  add_rows(
    rep(seq_len(n), times = k),
    as.vector(x_index),
    1,
    rep(1, n),
    rep(1, n)
  )
  # Each cell holds at most one point
  add_rows(
    rep(seq_len(k), each = n),
    as.vector(x_index),
    1,
    rep(0, k),
    rep(1, k)
  )

  cell_key <- paste0(cells$row, "_", cells$col)
  offsets <- neighbour_offsets(adjacency_diagonal)
  neighbour_index <- function(kk) {
    nk <- paste0(
      cells$row[kk] + offsets[, 1],
      "_",
      cells$col[kk] + offsets[, 2]
    )
    idx <- match(nk, cell_key)
    idx[!is.na(idx)]
  }

  if (adjacency_weight != 0) {
    occ <- n_vars + seq_len(k)
    obj <- c(obj, rep(0, k))
    n_vars <- n_vars + k
    add_rows(
      c(rep(seq_len(k), each = n), seq_len(k)),
      c(as.vector(x_index), occ),
      c(rep(1, n * k), rep(-1, k)),
      rep(0, k),
      rep(0, k)
    )

    pairs <- do.call(
      rbind,
      lapply(seq_len(k), \(kk) {
        nb <- neighbour_index(kk)
        if (length(nb)) cbind(kk, nb) else NULL
      })
    )
    pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
    n_pairs <- nrow(pairs)
    if (n_pairs > 0) {
      pij <- n_vars + seq_len(n_pairs)
      obj <- c(obj, rep(-adjacency_weight, n_pairs))
      n_vars <- n_vars + n_pairs
      add_and_constraints(add_rows, pij, occ[pairs[, 1]], occ[pairs[, 2]])
    }
  }

  pairwise_vars <- 0L
  if (!is.null(pairwise)) {
    edges <- pairwise$edges
    edges <- edges[order(-edges$weight), , drop = FALSE]
    edge_limit <- pairwise$edge_limit %||% nrow(edges)
    edges <- edges[seq_len(min(edge_limit, nrow(edges))), , drop = FALSE]
    nearest_k <- min(pairwise$nearest_cells, k)

    cdist <- sqrt(
      outer(norm_x, centers$x, "-")^2 + outer(norm_y, centers$y, "-")^2
    )
    nearest <- lapply(seq_len(n), \(p) order(cdist[p, ])[seq_len(nearest_k)])
    neighbours <- lapply(seq_len(k), neighbour_index)

    ya <- integer()
    yb <- integer()
    yw <- numeric()
    for (e in seq_len(nrow(edges))) {
      a <- edges$from[e]
      b <- edges$to[e]
      if (a < 1 || a > n || b < 1 || b > n) {
        next
      }
      cells_b <- nearest[[b]]
      for (ca in nearest[[a]]) {
        for (cb in neighbours[[ca]]) {
          if (cb %in% cells_b) {
            ya <- c(ya, x_index[a, ca])
            yb <- c(yb, x_index[b, cb])
            yw <- c(yw, edges$weight[e] %||% 1)
          }
        }
      }
    }
    pairwise_vars <- length(ya)
    if (pairwise_vars > 0) {
      yv <- n_vars + seq_len(pairwise_vars)
      obj <- c(obj, -pairwise$weight * yw)
      n_vars <- n_vars + pairwise_vars
      add_and_constraints(add_rows, yv, ya, yb)
    } else {
      cli::cli_warn(
        c(
          "Pairwise adjacency produced no auxiliary variables.",
          i = "Consider increasing {.arg nearest_cells} or {.arg edge_limit}."
        ),
        call = call
      )
    }
  }

  A <- Matrix::sparseMatrix(
    i = unlist(ci),
    j = unlist(cj),
    x = rep_len_each(cv, ci),
    dims = c(n_rows, n_vars)
  )
  control <- do.call(highs::highs_control, solver_options)
  res <- highs::highs_solve(
    L = obj,
    lower = rep(0, n_vars),
    upper = rep(1, n_vars),
    A = A,
    lhs = unlist(lhs),
    rhs = unlist(rhs),
    types = rep("I", n_vars),
    control = control
  )

  if (length(res$primal_solution) != n_vars || anyNA(res$primal_solution)) {
    cli::cli_abort(
      c(
        "The MIP solver did not return a solution.",
        i = "HiGHS status: {res$status_message}."
      ),
      call = call
    )
  }

  sol <- matrix(res$primal_solution[seq_len(n * k)] > 0.5, nrow = n)
  assigned <- apply(sol, 1, \(r) if (any(r)) which(r)[1] else NA_integer_)

  list(
    cell = assigned,
    score = res$objective_value,
    status = res$status_message,
    feasible = identical(res$status_message, "Optimal"),
    pairwise_vars = pairwise_vars,
    pairwise_constraints = 3L * pairwise_vars
  )
}

# Linearises z = a AND b:  z <= a,  z <= b,  a + b - z <= 1.
add_and_constraints <- function(add_rows, z, a, b) {
  m <- length(z)
  idx <- seq_len(m)
  add_rows(
    c(idx, idx),
    c(z, a),
    c(rep(1, m), rep(-1, m)),
    rep(-Inf, m),
    rep(0, m)
  )
  add_rows(
    c(idx, idx),
    c(z, b),
    c(rep(1, m), rep(-1, m)),
    rep(-Inf, m),
    rep(0, m)
  )
  add_rows(
    c(idx, idx, idx),
    c(a, b, z),
    c(rep(1, m), rep(1, m), rep(-1, m)),
    rep(-Inf, m),
    rep(1, m)
  )
}

rep_len_each <- function(values, along) {
  unlist(Map(\(v, i) rep_len(v, length(i)), values, along))
}
