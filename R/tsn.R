as.tsn <- function(x) {
  UseMethod("as.tsn")
}

as.tsn.data.frame <- function(x) {
  as_tsn(as.numeric(x[[1L]]), seq_len(nrow(x)))
}

as.tsn.ts <- function(x) {
  as_tsn(as.numeric(x), stats::time(x))
}

as.tsn.numeric <- function(x) {
  as_tsn(x, seq_along(x))
}

as.tsn.default <- function(x) {
  cls <- class(x)
  stop_(
    "Unable to coerce an object of class {.cls cls} into a {.cls tsn} object."
  )
}

as.tsn.tsn <- identity

as_tsn <- function(x, time) {
  structure(
    list(
      timeseries = data.frame(
        series = 1L,
        time = time,
        value = x
      ),
      network = NULL
    ),
    id_col = "series",
    value_col = "value",
    time_col = "time",
    class = "tsn"
  )
}

get_ts <- function(x) {
  x$timeseries[[attr(x, "value_col")]]
}
