#' Apply rolling functions
#'
#' @param fun The function to apply.
#' @param data The data values
#' @param time The time values
#' @param width Window width.
#' @param ... Arguments passed to `fun`.
#' @noRd
roll <- function(fun, data, time, window, align, ...) {
  n <- length(data)
  out <- rep(NA, n)
  left <- (window - 1L) %/% 2
  right <- window - 1L - left
  start <- 1L
  end <- n
  if (align == "center") {
    start <- 1 + left
    end <- n - right
    w <- seq(start - left, start + right)
  } else if (align == "right") {
    start <- window
    w <- seq(start - window + 1L, start)
  } else {
    end <- n - window + 1
    w <- seq(start, start + window - 1L)
  }
  for (i in seq(start, end)) {
    tmp <- ifelse_(
      missing(time),
      try_(fun(data[w], ...)),
      try_(fun(data = data[w], time = time[w], ...))
    )
    if (!inherits(tmp, "try-error") && is.finite(tmp)) {
      out[i] <- tmp
    }
    w <- w + 1L
  }
  out
}

# Residual Mean Square Differences
rmsqd <- function(x) {
  sqrt(mean(diff(x)^2, na.rm = TRUE))
}

rescale <- function(x, scale) {
  x_range <- range(x, na.rm = TRUE, finite = TRUE)
  (x - x_range[1L]) / (x_range[2L] - x_range[1L]) *
    (scale[2] - scale[1]) + scale[1L]
}

# Functions borrowed from the `dynamite` and `tna` packages -------------------

#' Shorthand for `try(., silent = TRUE)`
#'
#' @param expr An \R expression to try.
#' @noRd
try_ <- function(expr) {
  try(expr, silent = TRUE)
}

# Define the null coalescing operator for older R versions
if (base::getRversion() < "4.4.0") {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

#' Number of unique elements in a vector
#'
#' @param x A `vector`.
#' @noRd
n_unique <- function(x) {
  length(unique(x))
}

#' Shorthand for `if (test) yes else no`
#'
#' @param test A `logical` value of the condition to evaluate.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @param no An \R object to return when `test` evaluates to `FALSE`.
#' @noRd
ifelse_ <- function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}

#' Return `yes` if `test` is `TRUE`, otherwise return `NULL`
#'
#' @param test \[`logical(1)`] Condition to evaluate.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @noRd
onlyif <- function(test, yes) {
  if (test) {
    yes
  } else {
    NULL
  }
}

#' Generate a Warning Message
#'
#' @param message See [cli::cli_warn()].
#' @param ... See [cli::cli_warn()].
#' @noRd
warning_ <- function(message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}

#' Stop Function Execution Without Displaying the Call
#'
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}

#' Stop function execution unless a condition is true
#'
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See [cli::cli_abort()].
#' @noRd
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}

#' Generate an Informative Message
#'
#' @param message See [cli::cli_inform()]
#' @param ... See [cli::cli_inform()]
#' @noRd
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}

#' Create a Comma-separated Character String
#'
#' @param x A `character` vector.
#' @noRd
cs <- function(...) {
  paste0(c(...), collapse = ", ")
}
