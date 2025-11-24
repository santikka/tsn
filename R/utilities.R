#' Apply rolling functions
#'
#' @param fun The function to apply.
#' @param values The data values
#' @param time The time values
#' @param window Window width.
#' @param align Window alignment.
#' @param ... Arguments passed to `fun`.
#' @noRd
roll <- function(fun, values, time, window, align, ...) {
  n <- length(values)
  out <- rep(NA, n)
  left <- (window - 1L) %/% 2
  right <- window - 1L - left
  start <- 1L
  end <- n
  if (align == "center") {
    start <- 1 + left
    end <- n - right
    w <- (start - left):(start + right)
  } else if (align == "right") {
    start <- window
    w <- (start - window + 1L):start
  } else {
    end <- n - window + 1
    w <- start:(start + window - 1L)
  }
  if (missing(time)) {
    for (i in start:end) {
      out[i] <- fun(values[w], ...)
      w <- w + 1L
    }
  } else {
    for (i in start:end) {
      out[i] <- fun(values = values[w], time = time[w], ...)
      w <- w + 1L
    }
  }
  out
}

#' Compute rolling mean
#' @param values The data values
#' @param window Window width.
#' @param align Window alignment.
#' @param ... Arguments passed to `sum`.
#' @noRd
rollmean <- function(values, window, align, ...) {
  n <- length(values)
  out <- rep(NA, n)
  left <- (window - 1L) %/% 2
  right <- window - 1L - left
  start <- 1L
  end <- n
  if (align == "center") {
    start <- 1 + left
    end <- n - right
    w <- (start - left):(start + right)
  } else if (align == "right") {
    start <- window
    w <- (start - window + 1L):start
  } else {
    end <- n - window + 1
    w <- start:(start + window - 1L)
  }
  idx <- start:end
  out[idx] <- values[window:n] - values[c(1, seq_len(n - window))]
  out[start] <- sum(values[w], ...)
  out[idx] <- cumsum(out[idx]) / window
  out
}

#' Statistical mode
#'
#' @param x A `vector`.
#' @noRd
stat_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  if (length(ux) == 0) {
    NA
  } else {
    ux[which.max(tabulate(match(x, ux)))]
  }
}

#' Compute entropy by binning
#'
#' @param x A `numeric` vector.
#' @param bins An `integer` for the number of bins.
#' @noRd
entropy <- function(x, bins) {
  x <- x[!is.na(x)]
  r <- range(x, na.rm = TRUE)
  breaks <- seq(r[1L] - 1e-9, r[2L] + 1e-9 , length.out = bins + 1)
  counts <- graphics::hist(x, breaks = breaks, plot = FALSE)$counts
  prob <- counts / sum(counts)
  prob <- prob[prob > 0]
  -sum(prob * log2(prob))
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

#' Default value operator for a missing argument
#'
#' @param x An \R object
#' @param y An \R object to assign if `x` is missing
#' @noRd
`%m%` <- function(x, y) {
  if (missing(x)) y else x
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
