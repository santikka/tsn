#' Time Series Network Construction
#'
#' Build networks from distance matrices using various construction methods
#'
#' @export
#' @param data \[`tsn`, `ts`, `data.frame`, `numeric()`]\cr Time-series data.
#' @param measure \[`character(1)`]\cr The distance measure to use.
#' @param window \[`integer(1)`]\cr Width of the sliding window (default 1).
#' @param step \[`integer(1)`]\cr Step size for windows (default 1)
#' @param pairwise \[`logical(1)`]\cr If `TRUE` (default), calculates all
#'   pairwise distances between windows. If `FALSE`, uses consecutive windows
#'   only.
#' @param symmetric \[`logical(1)`]\cr If `TRUE` (default), ensures that the
#'   distance matrix is symmetric. If `FALSE`, keeps edge directions
#'   when `pairwise = FALSE`.
#' @param normalize \[`logical(1)`]\cr If `TRUE`, the distance matrix is
#'   max-normalized before network construction. If `FALSE`, uses the
#'   distances as is (default).
#' @param method \[`character(1)`]\cr Network construction method
#'   The available options are:
#'
#'   * `"full"`: Uses the raw similarity matrix (default).
#'   * `"knn"`: K-nearest neighbors network (KNN).
#'   * `"percentile"`: Percentile threshold network.
#'   * `"threshold"`: Distance threshold network.
#'   * `"gaussian"`: Gaussian kernel network.
#'
#' @param ... Additional arguments passed to the network construction method.
#'   These include:
#'
#'   * `k` Number of nearest neighbors (for KNN method)
#'   * `percentile` Percentile threshold (for percentile method)
#'   * `threshold` Distance threshold (for threshold method)
#'   * `sigma` Gaussian kernel parameter (for Gaussian method)
#'
#' @return TODO
#' @examples
#' set.seed(123)
#' x <- cumsum(rnorm(100))
#' net <- build_network(x, window = 10, step = 10)
#'
build_network <- function(data, measure = "euclidean", window = 1L, step = 1L,
                          pairwise = TRUE, symmetric = TRUE, normalize = FALSE,
                          method = "full", ...) {
  check_missing(data)
  measure <- check_match(measure, get_distance_measures())
  method <- check_match(method, get_network_methods())
  window <- default(window, 1L)
  step <- default(step, 1L)
  check_flag(pairwise)
  check_flag(symmetric)
  check_flag(normalize)
  check_network_dots(...)
  data <- as.tsn(data)
  values <- get_values(data)
  dist_mat <- distance_matrix(
    values, measure, window, step, pairwise, symmetric, ...
  )
  if (normalize) {
    max_dist <- max(dist_mat[is.finite(dist_mat)])
    dist_mat <- dist_mat / max_dist
  }
  adj_mat <- network_funs[[method]](dist_mat, ...)
  adj_mat
}

# Methods for building the network based on distances ---------------------

network_funs <- list()

network_funs$full <- function(x, ...) {
  out <- x / (1 + x)
  diag(out) <- 0
  out
}

network_funs$threshold <- function(x, threshold, ...) {
  n <- nrow(x)
  out <- matrix(0L, n, n)
  out[x <= threshold & x > 0] <- 1L
  out
}

network_funs$knn <- function(x, k, ...) {
  n <- nrow(x)
  out <- matrix(0L, n, n)
  for (i in seq_len(n)) {
    neighbors <- x[i,-i]
    if (all(is.finite(neighbors))) {
      out[i, order(neighbors)[1:k]] <- 1L
    }
  }
  out
}

network_funs$percentile <- function(x, percentile, ...) {
  q <- quantile(x[x > 0], percentile, na.rm = TRUE)
  x[x <= q & x > 0] <- 1
}

network_funs$gaussian <- function(x, sigma, ...) {
  if (missing(sigma)) {
    sigma <- stats::median(x[x > 0], na.rm = TRUE)
  }
  out <- exp(-x^2 / (2 * sigma^2))
  diag(out) <- 0L
  threshold <- exp(-4) # TODO argument?
  out[out < threshold] <- 0
  out
}

get_network_methods <- function() {
  names(network_funs)
}
