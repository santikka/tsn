# Some validation functions from the `tna` package.

#' Check that `x` is of specific class
#'
#' @param x An \R object.
#' @inheritParams class
#' @noRd
check_class <- function(x, what) {
  arg <- deparse(substitute(x))
  stopifnot_(
    inherits(x, what),
    "Argument {.arg {arg}} must be a {.cls {what}} object."
  )
}

#' Check that columns can be found in the data
#'
#' @param cols A `character` vector of column names.
#' @param data_names Column names of a data frame.
#' @noRd
check_cols <- function(cols, data_names) {
  cols_obs <- cols %in% data_names
  cols_mis <- cols[!cols_obs]
  stopifnot_(
    all(cols_obs),
    c(
      "The columns {.val {cols}} must exist in the data.",
      `x` = "The following columns were
             not found in the data: {.val {cols_mis}}."
    )
  )
}

#' Check That `x` is a Logical Value
#'
#' @param x An \R object expected to be a `logical` value.
#' @noRd
check_flag <- function(x) {
  arg <- deparse(substitute(x))
  stopifnot_(
    checkmate::test_flag(x = x),
    "Argument {.arg {arg}} must be a single {.cls logical} value."
  )
}

#' Check if argument matches given choices ignoring case
#'
#' @param x A `character` string.
#' @inheritParams match.arg
#' @noRd
check_match <- function(x, choices, several.ok = FALSE) {
  arg <- deparse(substitute(x))
  x <- onlyif(is.character(x), tolower(x))
  x <- try_(match.arg(arg = x, choices = choices, several.ok = several.ok))
  n_choices <- length(choices)
  prefix <- ifelse_(
    several.ok,
    "Elements of",
    "Argument"
  )
  stopifnot_(
    !inherits(x, "try-error"),
    "{prefix} {.arg {arg}} must be either
    {cli::qty(n_choices)} {.or {.val {choices}}}."
  )
  x
}

#' Check if argument is missing
#'
#' @param x An \R object.
#' @noRd
check_missing <- function(x) {
  arg <- deparse(substitute(x))
  stopifnot_(
    !missing(x),
    "Argument {.arg {arg}} is missing."
  )
}

#' Check that `x` is between a minimum and a maximum value
#'
#' @param x An \R object expected to be a single  `numeric` or `integer` value.
#' @noRd
check_range <- function(x, type = "numeric", scalar = TRUE,
                        min = 0.0, max = 1.0) {
  arg <- deparse(substitute(x))
  prefix <- ifelse_(scalar, "be a single", "only contain")
  suffix <- ifelse_(
    scalar,
    ifelse_(type == "integer", "", " value"),
    " values"
  )
  test_fun <- ifelse_(
    type == "numeric",
    ifelse_(scalar, checkmate::test_number, checkmate::test_numeric),
    ifelse_(scalar, checkmate::test_int, checkmate::test_integer)
  )
  stopifnot_(
    test_fun(x = x, lower = min, upper = max),
    "Argument {.arg {arg}} must {prefix}
    {.cls {type}} {suffix} between {min} and {max}."
  )
}

#' Check if argument is a character string
#'
#' @param x An \R object.
#' @noRd
check_string <- function(x) {
  if (missing(x)) {
    return()
  }
  arg <- deparse(substitute(x))
  stopifnot_(
    is.character(x) && length(x) == 1L,
    "Argument {.arg {arg}} must be a {.cls character} vector of length 1."
  )
}



#' Check that `x` is a non-negative
#'
#' @param x An \R object expected to be a `numeric` or `integer`
#' value or a vector.
#' @param type A `character` string corresponding to
#' the type that `x` should be.
#' @param strict A `logical` value. If `FALSE` (the default), expects
#' non-negative values and positive otherwise.
#' @param scalar A `logical` value indicating if `x` should be expected
#' to be a single value.
#' @noRd
check_values <- function(x, type = "integer", strict = FALSE,
                         scalar = TRUE) {
  arg <- deparse(substitute(x))
  suffix <- ifelse_(
    scalar,
    ifelse_(type == "integer", "", " value"),
    " vector"
  )
  test_fun <- ifelse_(
    type == "numeric",
    ifelse_(scalar, checkmate::test_number, checkmate::test_numeric),
    ifelse_(scalar, checkmate::test_int, checkmate::test_integer)
  )
  strictness <- ifelse_(strict, "positive", "non-negative")
  stopifnot_(
    test_fun(x = x, lower = as.integer(strict)),
    "Argument {.arg {arg}} must be a {strictness} {.cls {type}}{suffix}."
  )
}

check_network_dots <- function(...) {
  dots <- list(...)
  if (!is.null(dots$k)) {
    check_values(k, strict = TRUE)
  }
  if (!is.null(dots$percentile)) {
    check_range(percentile)
  }
  if (!is.null(dots$threshold)) {
    check_values(threshold, type = "numeric")
  }
  if (!is.null(dots$sigma)) {
    check_values(sigma, type = "numeric", strict = TRUE)
  }
}
