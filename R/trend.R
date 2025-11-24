#' Compute trend classification based on various metrics
#'
#' Calculates rolling metrics for a time series and classifies trends
#' as ascending, descending, flat, or turbulent.
#'
#' Computes rolling metrics. Trend classifications ("ascending", "descending",
#' "flat") are first determined using `epsilon`. Then, a "turbulent"
#' classification can override these if the rolling volatility of the metric
#' exceeds a dynamically adjusted turbulence threshold. For segments initially
#' classified as "flat", this threshold is
#' `turbulence_threshold * flat_to_turbulent_factor`,
#' making them more stable against reclassification as "turbulent" due to minor
#' noise.
#'
#' @export
#' @param x \[`tsn`, `ts`, `data.frame`, `numeric()`]\cr Time-series data.
#' @param window \[`integer(1)`]\cr The window size for metric calculation.
#'   If missing, uses the following adaptive sizing:
#'   `max(3, min(length(data), round(length(data) / 10)))`.
#' @param method \[`character(1)`]\cr The method for trend calculation.
#'   The available options ares: `"slope"` (default) and `"growth_factor"`.
#' @param slope \[`character(1)`]\cr The method for slope calculation for
#'   `method = "slope"`. The available options are:
#'   `"ols"` (ordinary least squares), `"robust"` (Theil-Sen estimator),
#'   `"spearman"` (Spearman rank correlation based),
#'   and `"kendall"` (Kendall's tau based). Default: `"robust"`.
#' @param epsilon \[`numeric(1)`]\cr A threshold value for defining flat trends
#'   based on the metric value. For `method = "slope"`, values between
#'   `(-epsilon, +epsilon)` are considered flat. For `"growth_factor"`, values
#'   between `(1 - epsilon)` and `(1 + epsilon)` are considered flat.
#'   Default: `0.05`.
#' @param turbulence_threshold \[`numeric(1)`]\cr The baseline threshold value
#'   for classifying a segment as "turbulent". Based on a custom combined
#'   volatility metric (CV + 0.5 * range factor of rolling metric values).
#'   Default: 5.
#' @param flat_to_turbulent_factor \[`numeric(1)`]\cr A multiplier for
#'   `turbulence_threshold` when assessing if an already "flat" segment should
#'   be reclassified as "turbulent". A value > 1 makes "flat" trends more
#'   resistant to becoming "turbulent". Default: 1.5.
#' @param align \[`character(1)`]\cr Alignment of the window. The available
#'   options are: `"center"` (default), `"right"`, and `"left"`. The calculated
#'   metric is assigned to the center, rightmost, or leftmost point of the
#'   window, respectively.
#' @return A `tsn` object whose `series` column is a factor with
#'   the following classes: `"Ascending"`, `"Descending"`, `"Flat"`,
#'   `"Turbulent"`, `"Missing Data"`, or `"Initial"`.
#' @examples
#' set.seed(123)
#' x <- cumsum(rnorm(200)) # Longer series to see more varied trends
#' # Using a slightly larger epsilon to catch more "flat" regions
#' tr <- trend(
#'   x,
#'   window = 15,
#'   method = "slope",
#'   slope = "ols",
#'   epsilon = 0.1,
#'   turbulence_threshold = 5,
#'   flat_to_turbulent_factor = 2
#' )
#'
trend <- function(x, ...) {
  UseMethod("trend")
}

#' @export
#' @rdname trend
trend.ts <- function(x, ...) {
  df <- data.frame(value = as.numeric(x), id = 1L, time = stats::time(x))
  trend(
    x = tsn(df, "value", "id", "time"),
    ...
  )
}

#' @rdname trend
#' @export
trend.default <- function(x, ...) {
  df <- data.frame(value = as.numeric(x), id = 1L, time = seq_along(x))
  trend(
    x = tsn(df, "value", "id", "time"),
    ...
  )
}

#' @rdname trend
#' @export
trend.data.frame <- function(x, value_col, id_col, time_col, ...) {
  check_missing(x)
  check_missing(value_col)
  check_class(x, "data.frame")
  check_string(value_col)
  check_string(id_col)
  check_string(time_col)
  trend(
    tsn(x, value_col, id_col, time_col),
    ...
  )
}

#' @rdname trend
#' @export
trend.tsn <- function(x, window, method = "slope", slope = "robust",
                      epsilon = 0.05, turbulence_threshold = 5,
                      flat_to_turbulent_factor = 1.5, align = "center") {
  check_missing(x)
  check_class(x, "tsn")
  method <- check_match(method, c("slope", "growth_factor"))
  slope <- check_match(slope, c("ols", "robust", "spearman", "kendall"))
  align <- check_match(align, c("center", "right", "left"))
  check_values(epsilon, type = "numeric")
  check_values(turbulence_threshold, type = "numeric")
  check_values(flat_to_turbulent_factor, type = "numeric")
  ids <- unique(x$id)
  n <- min(table(x$id))
  window <- ifelse_(missing(window), max(3, min(n, round(n / 10))), window)
  stopifnot_(
    window > 2 && window < n,
    "Argument {.arg window} must be between 2 and {n}
    (the number of observations)."
  )
  window_method <- ifelse_(method == "slope", slope, method)
  trend_codes <- character(nrow(x))
  for (i in seq_along(ids)) {
    idx <- which(x$id == ids[i])
    values <- x$value[idx]
    time <- x$time[idx]
    n <- length(values)
    state <- rep("Initial", n)
    metric_values <- roll(
      fun = metric_funs[[window_method]],
      values = values,
      time = time,
      window = window,
      align = align
    )
    values_na <- is.na(values)
    state[values_na] <- "Missing Data"
    neutral_val <- ifelse_(method == "growth_factor", 1, 0)
    lower <- neutral_val - epsilon
    upper <- neutral_val + epsilon
    valid <- !is.na(metric_values) & !values_na
    valid_metrics <- metric_values[valid]
    state[valid] <- ifelse(
      valid_metrics > upper,
      "Ascending",
      ifelse(
        valid_metrics < lower,
        "Descending",
        "Flat"
      )
    )
    volatility_window <- min(max(3, window %/% 2), sum(valid))
    valid_idx <- which(valid)
    n_valid <- length(valid_metrics)
    if (n_valid >= volatility_window) {
      for (j in seq(volatility_window, n_valid)) {
        window_idx <- seq(j - volatility_window + 1L, j)
        window_metric <- valid_metrics[window_idx]
        if (sum(!is.na(window_metric)) < 2L) {
          next
        }
        metric_sd <- stats::sd(window_metric, na.rm = TRUE)
        metric_am <- abs(mean(window_metric, na.rm = TRUE))
        metric_range <- diff(range(window_metric, na.rm = TRUE))
        if (is.na(metric_sd) ||
            is.na(metric_am) ||
            is.na(metric_range) ||
            metric_sd == 0 || metric_am == 0) {
          next
        }
        volatility_cv <- metric_sd / metric_am
        volatility_range_factor <- metric_range / metric_am
        combined_vol <- volatility_cv + 0.5 * volatility_range_factor
        k <- valid_idx[j]
        base_trend <- state[k]
        effective_threshold <- ifelse_(
          base_trend == "Flat",
          turbulence_threshold * flat_to_turbulent_factor,
          turbulence_threshold
        )
        if (base_trend != "Missing Data" &&
            combined_vol > effective_threshold) {
          state[k] <- "Turbulent"
        }
      }
    }
    trend_codes[idx] <- state
  }
  x$state <- factor(
    trend_codes,
    levels = c(
      "Ascending",
      "Descending",
      "Flat",
      "Turbulent",
      "Missing Data",
      "Initial"
    )
  )
  x
}


# Metrics -----------------------------------------------------------------

metric_funs <- list()

metric_funs$growth_factor <- function(values, time) {
  values <- values[!is.na(y)]
  n <- length(values)
  if (n < 2) {
    return(NA_real_)
  }
  values[n] / values[1L]
}

metric_funs$ols <- function(values, time) {
  stats::cov(time, values) / stats::var(time)
}

metric_funs$robust <- function(values, time) {
  n <- length(time)
  s <- numeric(n * (n - 1) %/% 2L)
  idx <- 0L
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      idx <- idx + 1L
      s[idx] <- (values[j] - values[i]) / (time[j] - time[i])
    }
  }
  stats::median(s, na.rm = TRUE)
}

metric_funs$spearman <- function(values, time) {
  corr <- stats::cor(time, values, method = "spearman", use = "complete.obs")
  corr * stats::sd(values, na.rm = TRUE) / stats::sd(time, na.rm = TRUE)
}

metric_funs$kendall <- function(values, time) {
  corr <- stats::cor(time, values, method = "kendall", use = "complete.obs")
  corr * stats::sd(values, na.rm = TRUE) / stats::sd(time, na.rm = TRUE)
}
