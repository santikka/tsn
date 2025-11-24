tsn <- function(x, value_col, id_col, time_col) {
  cols_req <- c(
    value_col,
    onlyif(!missing(id_col), id_col),
    onlyif(!missing(time_col), time_col)
  )
  check_cols(cols_req, names(x))
  x <- x[, cols_req, drop = FALSE]
  if (missing(id_col)) {
    id_col <- "id"
    x$id <- 1L
  }
  x <- x |>
    dplyr::group_by(!!rlang::sym(id_col))
  if (missing(time_col)) {
    time_col <- "time"
    x <- x |>
      dplyr::mutate(time = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col))
  } else {
    x <- x |>
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col)) |>
      dplyr::ungroup()
  }
  x <- x |>
    dplyr::rename(
      value = !!rlang::sym(value_col),
      id = !!rlang::sym(id_col),
      time = !!rlang::sym(time_col)
    )
  structure(
    x,
    class = c("tsn", "data.frame")
  )
}

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
      id = 1L,
      value = x,
      time = time
    ),
    class = c("tsn", "data.frame")
  )
}

