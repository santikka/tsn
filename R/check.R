# Some validation functions from the `tna` package.

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
