#' Time Series Distance Matrix Calculation
#'
#' Calculate distances between time series windows using various measures
#'
#' @param values Time-series data values
#' @inheritParams build_network
#' @param ... Additional arguments passed to distance calculation
#' @return Distance matrix
#' @noRd
distance_matrix <- function(values, measure = "euclidean", window,
                            step, pairwise = TRUE, symmetric = TRUE,
                            ...) {
  n <- length(values)
  step <- ifelse_(missing(step), 1L, step)
  windows <- generate_windows(values, window, step)
  n_windows <- length(windows)
  dist_mat <- matrix(Inf, n_windows, n_windows)
  diag(dist_mat) <- 0.0
  if (pairwise) {
    for (i in 1L:(n_windows - 1L)) {
      for (j in (i + 1L):n_windows) {
        dist_mat[i, j] <- calculate_distance(
          windows[[i]],
          windows[[j]],
          measure,
          ...
        )
      }
    }
  } else {
    for (i in 1:(n_windows - 1)) {
      j <- i + 1L
      dist_mat[i, j] <- calculate_distance(
        windows[[i]],
        windows[[j]],
        method,
        ...
      )
    }
  }
  if (symmetric) {
    dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
  }
  window_names <- paste0("W", 1:n_windows)
  rownames(dist_mat) <- window_names
  colnames(dist_mat) <- window_names
  dist_mat
}

get_distance_measures <- function(all = FALSE) {
  c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  )
}

#' Calcuate distances between two vectors
#'
#' @param x The first vector.
#' @param y The fecond vector.
#' @param measure Distance measure.
#' @param ... Additional arguments passed to the distance method.
#' @return Distance value.
#' @noRd
calculate_distance <- function(x, y, measure, ...) {
  dist_measures <- get_distance_measures()
  if (measure %in% dist_measures) {
    dist <- as.numeric(stats::dist(rbind(x, y), method = measure, ...))
  } # TODO other methods
  dist
}

#' Generate sliding windows from time-series data
#'
#' @param values Time-series data as a vector
#' @param window Size of each window
#' @param step Step size between windows
#' @noRd
generate_windows <- function(values, window, step) {
  n <- length(values)
  k <- n - window
  m <- k %/% step + 1L
  partial <- (k %% step) > 0
  windows <- vector(mode = "list", length = m + as.integer(partial))
  idx <- 1L
  w <- 1:window
  for (i in seq_len(m)) {
    windows[[idx]] <- values[w]
    idx <- idx + 1L
    w <- w + step
  }
  if (partial) {
    windows[[m + 1L]] <- values[(n - step):n]
  }
  windows
}
