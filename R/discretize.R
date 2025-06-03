#' Convert Time-Series Data into Wide Format Sequence Data
#'
#' Converts time-series data into sequence data via discretization.
#' Various methods for discretization are available including gaussian
#' mixtures, K-means clustering and kernel density based binning.
#'
#' @export
#' @param data \[`data.frame`, `ts`, `tsn`, `numeric()`]\cr Either time-series
#'   data in long format (`data.frame`, `tsn`), a time-series object (`ts`), or
#'   a vector of values.
#' @param id_col \[`character(1)`]\cr The name of the column that contains
#'   the unique identifiers.
#' @param value_col \[`character(1)`]\cr The name of the column that contains
#'   the data values.
#' @param order_col \[`character(1)`]\cr The name of the column that contains
#'   the time values (not required if the data is already in order),
#' @param n_states \[`integer(1)`]\cr The number of states to discretize the
#'   data into.
#' @param method \[`character(1)`]\cr The name of the discretization method to
#'   use.
#'
#'   * `kmeans`: for K-means clustering (the default).
#'   * `width`: for equal width binning.
#'   * `quantile`: for quantile-based binning.
#'   * `kde`: for binning based on kernel density estimation.
#'   * `gaussian`: for a Gaussian mixture model.
#'
#' @param labels \[`character()`]\cr A vector of names for the states. The
#'   length must be `n_states` The defaults is consecutive numbering,
#'   i.e. `1:n_states`.
#' @param unused_fn \[`function`]\cr How to handle extra columns when pivoting
#'   to wide format. See [tidyr::pivot_wider()]. The default is to keep all
#'   columns and to use the first value.
#' @param ... Additional arguments passed to the discretization method
#'   ([stats::kmeans()] for `kmeans`, [stats::density()] and
#'   [pracma::findpeaks()] for `kde`, and
#'   [mclust::Mclust()] for `gaussian`).
#' @return TODO
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @examples
#' # Long format data
#' ts_data <- data.frame(
#'   id = gl(10, 100),
#'   series = c(
#'     replicate(
#'       10,
#'       stats::arima.sim(list(order = c(2, 1, 0), ar = c(0.5, 0.2)), n = 99)
#'     )
#'   )
#' )
#' data <- discretize(ts_data, "id", "series", n_states = 3)
#'
#' # Time-series data
#' data <- discretize(EuStockMarkets, n_states = 3)
#'
discretize <- function(data, n_states, labels = 1:n_states,
                       method = "kmeans", unused_fn = dplyr::first, ...) {
  UseMethod("discretize")
}

#' @export
#' @rdname discretize
discretize.ts <- function(data, n_states, labels = 1:n_states,
                          method = "kmeans", unused_fn = dplyr::first, ...) {
  df <- data.frame(
    series = 1L,
    value = as.numeric(data),
    time = stats::time(data)
  )
  discretize(
    df,
    id_col = "series",
    value_col = "value",
    order_col = "time",
    n_states = n_states,
    labels = labels,
    method = method,
    unused_fn = unused_fn,
    ...
  )
}

#' @export
#' @rdname discretize
discretize.default <- function(data, id_col, value_col, order_col, n_states,
                               labels = 1:n_states, method = "kmeans",
                               unused_fn = dplyr::first, ...) {
  check_missing(data)
  check_missing(value_col)
  check_class(data, "data.frame")
  check_string(id_col)
  check_string(value_col)
  labels <- try_(as.character(labels))
  stopifnot_(
    checkmate::test_character(
      x = labels,
      any.missing = FALSE,
      min.len = n_states,
      max.len = n_states
    ),
    "Argument {.arg labels} must be coercible to {.cls character} and
     provide a label for each state."
  )
  stopifnot_(
    checkmate::test_int(x = n_states, lower = 2L),
    "Argument {.arg n_states} must be an integer greater than 1."
  )
  stopifnot_(
    length(labels) == n_states,
    "Argument {.arg labels} must have length {n_states}
    (same as {.arg n_states})."
  )
  method <- check_match(method, names(discretization_funs))
  cols_req <- c(
    value_col,
    onlyif(!missing(id_col), id_col),
    onlyif(!missing(order_col), order_col)
  )
  check_cols(cols_req, names(data))
  if (missing(id_col)) {
    id_col <- ".id"
    data$.id <- 1L
  }
  complete <- stats::complete.cases(data[, c(id_col, value_col)])
  values <- data[[value_col]][complete]
  # TODO warn if number of states is less than n_states
  discretized <- discretization_funs[[method]](values, n_states, ...)
  states <- discretized$states
  output <- discretized$output
  state_col <- paste0(value_col, "_state")
  data[[state_col]] <- NA
  data[[state_col]][complete] <- states
  data[[state_col]] <- factor(
    data[[state_col]],
    levels = seq_len(n_states),
    labels = labels
  )
  # stats <- compute_state_statistics(data, id_col, value_col, state_col)
  data <- data |>
    dplyr::group_by(!!rlang::sym(id_col))
  if (missing(order_col)) {
     data <- data |>
      dplyr::mutate(.time = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::arrange(!!rlang::sym(id_col), .time)
  } else {
    data <- data |>
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(order_col)) |>
      dplyr::mutate(.time = dplyr::row_number()) |>
      dplyr::ungroup()
  }
  # wide_data <- tidyr::pivot_wider(
  #   timed_data,
  #   id_cols = !!rlang::sym(id_col),
  #   names_from = !!rlang::sym(".time"),
  #   names_prefix = "T",
  #   values_from = !!rlang::sym(state_col),
  #   unused_fn = unused_fn
  # )
  # data$.id <- NULL
  # wide_data$.id <- NULL
  # time_cols <- grepl("^T[0-9]+$", names(wide_data), perl = TRUE)
  # sequence_data <- wide_data[, time_cols]
  # meta_data <- wide_data[, !time_cols]
  structure(
    data,
    id_col = id_col,
    value_col = value_col,
    state_col = state_col,
    time_col = ifelse_(missing(order_col), ".time", order_col),
    output = output,
    class = c("tsn", "data.frame")
  )
}

# #' Calculate comprehensive statistics for states
# #'
# #' @inheritParams prepare_ts
# #' @param state_col A `character` string naming the column that contains the
# #' state information.
# #' @return A `list` of global and local statistics.
# #' @noRd
# compute_state_statistics <- function(data, id_col, value_col, state_col) {
#   # For R CMD Check
#   group_size_ <- NULL
#   global <- data |>
#     dplyr::group_by(!!rlang::sym(state_col)) |>
#     dplyr::summarize(
#       freq = dplyr::n(),
#       prop = dplyr::n() / nrow(data),
#       mean = mean(!!rlang::sym(value_col)),
#       median = stats::median(!!rlang::sym(value_col)),
#       sd = stats::sd(!!rlang::sym(value_col)),
#       min = min(!!rlang::sym(value_col)),
#       max = max(!!rlang::sym(value_col)),
#       q25 = unname(stats::quantile(!!rlang::sym(value_col), 0.25)),
#       q75 = unname(stats::quantile(!!rlang::sym(value_col), 0.75)),
#     )
#   if (id_col == ".id") {
#     local <- global
#   } else {
#     local <- data |>
#       dplyr::group_by(!!rlang::sym(id_col)) |>
#       dplyr::mutate(
#         group_size_ = dplyr::n()
#       ) |>
#       dplyr::group_by(!!rlang::sym(id_col), !!rlang::sym(state_col)) |>
#       dplyr::summarize(
#         freq = dplyr::n(),
#         prop = dplyr::n() / dplyr::first(group_size_),
#         mean = mean(!!rlang::sym(value_col)),
#         median = stats::median(!!rlang::sym(value_col)),
#         sd = stats::sd(!!rlang::sym(value_col)),
#         min = min(!!rlang::sym(value_col)),
#         max = max(!!rlang::sym(value_col)),
#         q25 = unname(stats::quantile(!!rlang::sym(value_col), 0.25)),
#         q75 = unname(stats::quantile(!!rlang::sym(value_col), 0.75)),
#       )
#   }
#   list(global = global, local = local)
# }

# Discretization function wrappers --------------------------------------------

discretization_funs <- list()

discretization_funs$width <- function(x, n_states, ...) {
  r <- range(x)
  k <- n_states + 1L
  breaks <- seq(r[1], r[2], length.out = k)
  breaks[1L] <- breaks[1L] - 1.0
  breaks[k] <- breaks[k] + 1.0
  list(
    states = cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE),
    output = breaks
  )
}

discretization_funs$quantile <- function(x, n_states, ...) {
  k <- n_states + 1L
  probs <- seq(0, 1, length.out = k)
  breaks <- stats::quantile(x, probs = probs)
  breaks[1L] <- breaks[1L] - 1.0
  breaks[k] <- breaks[k] + 1.0
  list(
    states = cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE),
    output = breaks
  )
}

discretization_funs$kde <- function(x, n_states, ...) {
  stopifnot_(
    requireNamespace("pracma", quietly = TRUE),
    "Please install the {.pkg pracma} package
     to use discretization based on kernel density estimation."
  )
  dots <- list(...)
  is_arg <- names(dots) %in% methods::formalArgs(pracma::findpeaks)
  density_args <- dots[!is_arg]
  density_args$x <- x
  dens <- do.call(stats::density, args = density_args)
  findpeaks_args <- dots[is_arg]
  findpeaks_args$npeaks <- n_states - 1L
  findpeaks_args$x <- -dens$y
  valleys <- do.call(pracma::findpeaks, args = findpeaks_args)
  breaks <- sort(unique(dens$x[valleys[, 2L]]))
  breaks <- c(min(x) - 1.0, breaks, max(x) + 1.0)
  list(
    states = cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE),
    output = breaks
  )
}

discretization_funs$gaussian <- function(x, n_states, ...) {
  stopifnot_(
    requireNamespace("mclust", quietly = TRUE),
    "Please install the {.pkg mclust} package
     to use gaussian mixture-based discretization."
  )
  model <- mclust::Mclust(
    data = x,
    G = n_states,
    modelNames = "V",
    verbose = FALSE
  )
  ord <- order(model$parameters$mean)
  list(
    states = ord[model$classification],
    output = model
  )
}

discretization_funs$kmeans <- function(x, n_states, ...) {
  km <- stats::kmeans(x, centers = n_states, ...)
  ord <- order(km$centers)
  list(
    states = ord[km$cluster],
    output = km
  )
}
