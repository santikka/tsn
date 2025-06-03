#' @export
as.tsn <- function(x) {
  UseMethod("as.tsn")
}

#' @export
as.tsn.data.frame <- function(x) {
  cls <- class(x)
  if ("tsn" %in% cls) {
    return(x)
  }
  stopifnot_(
    ncol(x) == 1L,
    "Time-series in data frame format must have a single column."
  )
  as_tsn(as.numeric(x[[1L]]), seq_len(nrow(x)))
}

#' @export
as.tsn.ts <- function(x) {
  as_tsn(as.numeric(x), stats::time(x))
}

#' @export
as.tsn.numeric <- function(x) {
  as_tsn(x, seq_along(x))
}

#' @export
as.tsn.default <- function(x) {
  cls <- class(x)
  if ("tsn" %in% cls) {
    return(x)
  }
  stop_(
    "Unable to coerce an object of class {.cls cls} into a {.cls tsn} object."
  )
}

as_tsn <- function(x, time) {
  structure(
    data.frame(
      series = 1L,
      time = time,
      value = x
    ),
    id_col = "series",
    value_col = "value",
    time_col = "time",
    class = c("tsn", "data.frame")
  )
}

get_values <- function(x) {
  x[[attr(x, "value_col")]]
}

get_time <- function(x) {
  x[[attr(x, "time_col")]]
}
