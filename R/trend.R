#' Compute trend classification based on various metrics
#'
#' Calculates rolling metrics for a time series and classifies trends
#' as ascending, descending, flat, or turbulent.
#'
#' @export
#' @param data \[`ts`, `data.frame`, `numeric()`]\cr Time-series data.
#' @param window_size \[`integer(1)`]\cr The window size for metric calculation.
#'   If missing, uses the following adaptive sizing:
#'   `max(3, min(length(data), round(length(data)/10)))`.
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
#'   between `(1-epsilon)` and `(1+epsilon)` are considered flat.
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
#' @return A `tsn` object with an added `trend` column for the `timeseries`
#'   data with the following classes: `"ascending"`, `"descending"`, `"flat"`,
#'   `"turbulent"`, `"Missing_Data"`, or `"Initial"`.
#'
#' @details Computes rolling metrics. Trend classifications
#' ("ascending", "descending", "flat") are first determined using `epsilon`.
#' Then, a "turbulent" classification can override these if the rolling
#' volatility of the metric exceeds a dynamically adjusted turbulence threshold.
#' For segments initially classified as "flat", this threshold is
#' `turbulence_threshold * flat_to_turbulent_factor`,
#' making them more stable against reclassification as "turbulent" due to minor
#' noise.
#'
#' @examples
#' set.seed(123)
#' x <- cumsum(rnorm(200)) # Longer series to see more varied trends
#' # Using a slightly larger epsilon to catch more "flat" regions
#' tr <- trend(
#'   x, window = 15, method = "slope",
#'   slope = "ols", epsilon = 0.1,
#'   turbulence_threshold = 5, flat_to_turbulent_factor = 2
#' )
#'
trend.default <- function(data, window, method = "slope", slope = "robust",
                          epsilon = 0.05, turbulence_threshold = 5,
                          flat_to_turbulent_factor = 1.5, align = "center") {
  data <- get_ts(as.tsn(data))
  method <- check_match(method, c("slope", "growth_factor"))
  slope <- check_match(slope, c("ols", "robust", "spearman", "kendall"))
  align <- check_match(align, c("center", "right", "left"))
  check_values(epsilon, type = "numeric")
  check_values(turbulence_threshold, type = "numeric")
  check_values(flat_to_turbulent_factor, type = "numeric")
  n <- length(data)
  window <- ifelse_(missing(window), max(3, min(n, round(n / 10))), window)
  stopifnot_(
    window > 2 && window < n,
    "Argument {.arg window} must be between 2 and {n}
    (the number of observations)."
  )

  time_idx <- ifelse_(stats::is.ts(data), stats::time(data), seq_along(data))
  # TODO try
  data <- as.numeric(data)

  # TODO test
  window_method <- ifelse_(method == "slope", slope, method)
  metric_values <- roll(
    fun = metric_funs[[window_method]],
    data = data,
    time = time,
    window = window,
    align = align
  )


  # if (align == "center" && sum(!is.na(metric_values)) > 0L) {
  #   first_valid <- min(which(!is.na(metric_values)))
  #   if (first_valid > 1L) {
  #     metric_values[seq(1L, first_valid - 1L)] <- metric_values[first_valid]
  #   }
  #   last_valid <- max(which(!is.na(metric_values)))
  #   if (last_valid < n) {
  #     metric_values[seq(last_valid + 1L, n)] <- metric_values[last_valid]
  #   }
  # }

  trend_codes <- rep("Initial", n)
  data_na <- is.na(data)
  trend_codes[data_na] <- "Missing_Data"
  neutral_val <- ifelse_(method == "growth_factor", 1, 0)
  lower <- neutral_val - epsilon
  upper <- neutral_val + epsilon
  valid <- !is.na(metric_values) & !data_na
  valid_metrics <- metric_values[valid]
  trend_codes[valid] <- ifelse(
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
    for (i in seq(volatility_window, n_valid)) {
      window_idx <- seq(i - volatility_window + 1L, i)
      window_metric <- valid_metrics[window_idx]
      if (sum(!is.na(window_metric)) < 2L) next
      metric_sd <- stats::sd(window_metric, na.rm = TRUE)
      metric_am <- abs(mean(window_metric, na.rm = TRUE))
      metric_range <- diff(range(window_metric, na.rm = TRUE))
      if (is.na(metric_sd) || is.na(metric_am) || is.na(metric_range)) next
      volatility_cv <- metric_sd / metric_am
      volatility_range_factor <- metric_range / metric_am
      combined_vol <- volatility_cv + 0.5 * volatility_range_factor
      j <- valid_idx[i]
      # Adjust turbulence threshold if current trend is "flat"
      base_trend <- trend_codes[j]
      effective_threshold <- ifelse_(
        base_trend == "Flat",
        turbulence_threshold * flat_to_turbulent_factor,
        turbulence_threshold
      )
      if (base_trend != "Missing_Data" && combined_vol > effective_threshold) {
        trend_codes[j] <- "Turbulent"
      }
    }
  }

  structure(
    list(
      timeseries = data.frame(
        id = 1L,
        time = time_idx,
        value = data,
        state = factor(
          trend_codes,
          levels = c(
            "Ascending",
            "Descending",
            "Flat",
            "Turbulent",
            "Missing_Data",
            "Initial"
          )
        )
      ),
      network = NULL
    ),
    id_col = "id",
    value_col = "value",
    state_col = "state",
    time_col = "time",
    class = "tsn"
  )
}


# Metrics -----------------------------------------------------------------

metric_funs <- list()

# metric_funs$ar1 <- function(x, y) {
#   y <- stats::na.omit(y)
#   model <- stats::ar.ols(
#     y, aic = FALSE, order.max = 1, demean = FALSE, intercept = TRUE
#   )
#   model$ar[1L]
# }

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
