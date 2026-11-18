# JavaScript's Math.round() rounds halves towards +Inf; R's round() does not.
js_round <- function(x) {
  floor(x + 0.5)
}

# Equivalent of JS `x || 1` for a numeric width/height.
or_one <- function(x) {
  if (is.na(x) || x == 0) 1 else x
}

calc_bounds <- function(x, y) {
  min_x <- min(x)
  max_x <- max(x)
  min_y <- min(y)
  max_y <- max(y)
  list(
    min_x = min_x,
    max_x = max_x,
    min_y = min_y,
    max_y = max_y,
    width = max_x - min_x,
    height = max_y - min_y
  )
}

check_number <- function(
  x,
  min = -Inf,
  max = Inf,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x > max) {
    range <- if (is.finite(min) && is.finite(max)) {
      paste0(" between ", min, " and ", max)
    } else if (is.finite(min)) {
      paste0(" >= ", min)
    } else {
      ""
    }
    if (is.numeric(x) && length(x) == 1) {
      cli::cli_abort(
        "{.arg {arg}} must be a single number{range}, not {x}.",
        call = call
      )
    }
    cli::cli_abort(
      "{.arg {arg}} must be a single number{range}, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  invisible(x)
}

check_count <- function(
  x,
  allow_null = FALSE,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (allow_null && is.null(x)) {
    return(invisible(x))
  }
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != round(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a single positive whole number, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  invisible(x)
}

check_flag <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be `TRUE` or `FALSE`, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  invisible(x)
}

# Coerce spacers given as a 2-column matrix/data frame of 1-based
# (row, col) positions into 0-based integer keys "r_c".
spacer_keys <- function(spacers, call = rlang::caller_env()) {
  if (is.null(spacers) || NROW(spacers) == 0) {
    return(character())
  }
  spacers <- as_spacer_matrix(spacers, call = call)
  paste0(spacers[, 1] - 1L, "_", spacers[, 2] - 1L)
}

as_spacer_matrix <- function(spacers, call = rlang::caller_env()) {
  if (is.null(spacers) || NROW(spacers) == 0) {
    return(matrix(integer(), ncol = 2, dimnames = list(NULL, c("row", "col"))))
  }
  if (is.data.frame(spacers)) {
    cols <- if (all(c("row", "col") %in% names(spacers))) {
      spacers[c("row", "col")]
    } else {
      spacers[1:2]
    }
    spacers <- as.matrix(cols)
  }
  if (!is.matrix(spacers) || ncol(spacers) < 2 || !is.numeric(spacers)) {
    cli::cli_abort(
      "{.arg spacers} must be a two-column matrix or data frame of {.field row} and {.field col} positions.",
      call = call
    )
  }
  out <- matrix(as.integer(spacers[, 1:2]), ncol = 2)
  colnames(out) <- c("row", "col")
  out
}

keys_to_spacer_matrix <- function(keys) {
  if (length(keys) == 0) {
    return(as_spacer_matrix(NULL))
  }
  parts <- do.call(rbind, strsplit(keys, "_", fixed = TRUE))
  out <- matrix(as.integer(parts) + 1L, ncol = 2)
  colnames(out) <- c("row", "col")
  out
}
