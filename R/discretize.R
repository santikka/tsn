#' Convert Time-Series Data into Wide Format Sequence Data
#'
#' Converts time-series data into sequence data via discretization.
#' Various methods for discretization are available including gaussian
#' mixtures, K-means clustering and kernel density based binning.
#'
#' @export
#' @param x \[`data.frame`, `ts`, `tsn`]\cr Either a time-series data object in
#'   long format (`data.frame`, `tsn`) or a time-series object (`ts`).
#' @param id_col \[`character(1)`]\cr The name of the column that contains
#'   the unique identifiers.
#' @param time_col \[`character(1)`]\cr The name of the column that contains
#'   the time values (not required if the data is already in order).
#' @param value_col \[`character(1)`]\cr The name of the column that contains
#'   the data values.
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
#' @return A `tsn` object which is a `data.frame` containing the original
#'   time series data and the discretized sequence of states.
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
#'
#' discr <- discretize(
#'   ts_data, id_col = "id", value_col = "series", n_states = 3
#' )
#'
#' # Time-series data
#' discr2 <- discretize(EuStockMarkets, n_states = 3)
#'
discretize <- function(x, ...) {
  UseMethod("discretize")
}

#' @export
#' @rdname discretize
discretize.default <- function(x, ...) {
  df <- data.frame(value = as.numeric(x), id = 1L, time = seq_along(x))
  discretize(x = tsn(df, "value", "id", "time"), ...)
}

#' @export
#' @rdname discretize
discretize.ts <- function(x, ...) {
  df <- data.frame(value = as.numeric(x), id = 1L, time = stats::time(x))
  discretize(x = tsn(df, "value", "id", "time"), ...)
}

#' @export
#' @rdname discretize
discretize.tsn <- function(x, ...) {
  discretize_(x, ...)
}

#' @export
#' @rdname discretize
discretize.data.frame <- function(x, value_col, id_col, time_col, n_states,
                                  labels = 1:n_states, method = "kmeans",
                                  unused_fn = dplyr::first, ...) {
  check_missing(x)
  check_missing(value_col)
  check_class(x, "data.frame")
  check_string(value_col)
  check_string(id_col)
  check_string(time_col)
  discretize_(
    x = tsn(x, value_col, id_col, time_col),
    n_states = n_states,
    labels = labels,
    method = method,
    unused_fn = unused_fn,
    ...
  )
}

discretize_ <- function(x, n_states, labels = 1:n_states, method = "kmeans",
                        unused_fn = dplyr::first, ...) {
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
  complete <- stats::complete.cases(x[, c("id", "value")])
  values <- x$value[complete]
  # TODO warn if number of states is less than n_states
  discretized <- discretization_funs[[method]](values, n_states, ...)
  states <- discretized$states
  output <- discretized$output
  x$state <- NA
  x$state[complete] <- states
  x$state <- factor(
    x$state,
    levels = seq_len(n_states),
    labels = labels
  )
  attr(x, "output") <- output
  x
}

# Discretization function wrappers --------------------------------------------

discretization_funs <- list()

discretization_funs$width <- function(x, n_states, ...) {
  r <- range(x)
  k <- n_states + 1L
  breaks <- seq(r[1L], r[2L], length.out = k)
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
  breaks <- c(min(x) - 1, breaks, max(x) + 1)
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
