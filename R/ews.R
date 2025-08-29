#' Detect Early Warning Signals in a Time Series
#'
#' @param data A `data.frame` that contains the time series data and time
#' points.
#' @param ts_col A `character` string naming the column for the actual
#' time series data values.
#' @param time_col A `character` string naming the column for the time points.
#' @param method A `character` string for the analysis method.
#' Either `"rolling"` or `"expanding"` for rolling window and expanding window,
#' respectively.
#' @param metrics A `character` vector for the EWS metrics to compute.
#' The default is `"all"` computing all metrics.
#' @param window A `numeric` value (percentage) for the window size of
#' the rolling window method.
#' @param burnin Integer. Burn-in period for the expanding method.
#' @param demean A `logical` value. If `TRUE` (the default), the `"ar1"`
#' metric will be based on an AR1 model where the mean of the observations
#' is first subtracted. See [stats::ar.ols()] for details.
#' @param detrend TODO
#' @param threshold A `numeric` value giving the z-score threshold for the
#' expanding window method..
#' @return An object of class `detect_warning` containing the EWS results,
#'   parameters, original data, and a summary.
#' @export
detect_warnings <- function(data, ts_col, time_col, method = "rolling",
                            metrics = "all", window = 50, burnin = 30,
                            detrend = "none", demean = TRUE,
                            threshold = 2, bandwidth,
                            span, degree, ...) {
  data <- as_tsn(data[[ts_col]], data[[time_col]])
  values <- get_values(data)
  time <- get_time(data)
  method <- check_match(method, c("rolling", "expanding"))
  available_metrics <- c("ar1", "sd", "skew", "kurt", "cv", "rr")
  metrics <- metrics %m% available_metrics
  metrics <- check_match(
    metrics,
    c(available_metrics, "all"),
    several.ok = TRUE
  )
  metrics <- ifelse_(
    "all" %in% metrics,
    available_metrics,
    metrics
  )
  check_range(window, min = 0.0, max = 100.0)
  check_range(burnin, min = 0.0, max = 100.0)
  check_flag(demean)
  detrend <- check_match(
    detrend,
    c("none", "gaussian", "loess", "linear", "first-diff")
  )
  window <- floor(0.01 * window * length(values))
  bandwidth <- bandwidth %m% round(window / 2)
  span <- span %m% 0.25
  degree <- degree %m% 2
  values <- detrend_ts(values, time, detrend, window, ...)
  ifelse_(
    method == "rolling",
    rolling_ews(values, time, metrics, window, demean),
    expanding_ews(values, time, metrics, burnin, demean, threshold)
  )
}

rolling_ews <- function(x, time, metrics, window, demean) {
  w <- window
  n <- length(x)
  m <- n - w + 1L
  idx <- seq_len(m)
  rolling_ar1 <- numeric(m)
  rolling_ar1_demean <- numeric(m)
  rolling_mean <- numeric(m)
  rolling_var <- numeric(m)
  rolling_skew <- numeric(m)
  rolling_kurt <- numeric(m)
  y <- x[1:w]
  s1 <- sum(y)
  s2 <- sum(y^2)
  s3 <- sum(y^3)
  s4 <- sum(y^4)
  s_cur <- s1 - y[1]
  s_lag <- s1 - x[w]
  s_lag2 <- s2 - x[w]^2
  s_prod <- sum(y[-1] * y[-w])
  mu <- s1 / w
  m2 <- (s2 - s1^2 / w) / w
  m3 <- (s3 - 3 * mu * s2 + 2 * w * mu^3) / w
  m4 <- (s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 3 * w * mu^4) / w
  mu_sq <- (w - 1) * mu^2
  rolling_ar1[1] <- ifelse_(
    demean,
    (s_prod - mu * (s_cur + s_lag) + mu_sq) / (s_lag2 - 2 * mu * s_lag + mu_sq),
    s_prod / s_lag2
  )
  rolling_mean[1] <- mu
  rolling_var[1] <- m2 * w / (w - 1)
  rolling_skew[1] <- m3 / (m2^(3/2))
  rolling_kurt[1] <- m4 / (m2^2)
  for (i in seq(2, m)) {
    x_new <- x[i + w - 1L]
    x_old <- x[i - 1L]
    s1 <- s1 + x_new - x_old
    s2 <- s2 + x_new^2 - x_old^2
    s3 <- s3 + x_new^3 - x_old^3
    s4 <- s4 + x_new^4 - x_old^4
    s_lag2 <- s_lag2 + x[i + w - 2L]^2 - x_old^2
    s_prod <- s_prod - x_old * x[i] + x_new * x[i + w - 2L]
    mu <- s1 / w
    m2 <- (s2 - s1^2 / w) / w
    m3 <- (s3 - 3 * mu * s2 + 2 * w * mu^3) / w
    m4 <- (s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 3 * w * mu^4) / w
    if (demean) {
      s_lag <- s_lag + x[i + w - 2L] - x_old
      s_cur <- s_cur + x_new - x[i]
      mu_sq <- (w - 1) * mu^2
      rolling_ar1[i] <- (s_prod - mu * (s_cur + s_lag) + mu_sq) /
        (s_lag2 - 2 * mu * s_lag + mu_sq)
    } else {
      rolling_ar1[i] <- s_prod / s_lag2
    }
    rolling_mean[i] <- mu
    rolling_var[i] <- m2 * w / (w - 1)
    rolling_skew[i] <- m3 / (m2^(3/2))
    rolling_kurt[i] <- m4 / (m2^2)
  }
  rolling_sd <- sqrt(rolling_var)
  rolling_metrics <- data.frame(
    time = time[w:n],
    ar1 = rolling_ar1,
    mean = rolling_mean,
    sd = rolling_sd,
    skew = rolling_skew,
    kurt = rolling_kurt,
    cv = rolling_sd / rolling_mean,
    rr = 1.0 - rolling_ar1
  )
  rolling_metrics <- rolling_metrics[, c("time", metrics)]
  kendall_tau <- apply(rolling_metrics[, -1, drop = FALSE], 2, function(z) {
    stats::cor.test(
      x = idx,
      y = z,
      alternative = "two.sided",
      conf.level = 0.95,
      method = "kendall"
    )$estimate
  })
  long <- tidyr::pivot_longer(
    rolling_metrics,
    cols = !(!!rlang::sym("time")),
    names_to = "metric",
    values_to = "score"
  ) |>
    dplyr::group_by(!!rlang::sym("metric")) |>
    dplyr::mutate(std = as.numeric(scale(!!rlang::sym("score")))) |>
    dplyr::ungroup()
  structure(
    long,
    orig_values = x,
    orig_time = time,
    cor = kendall_tau,
    method = "rolling",
    class = c("tsn_ews", "tbl_df", "tbl", "data.frame")
  )
}

expanding_ews <- function(x, time, metrics, burnin, demean, threshold) {
  w <- burnin + 1
  n <- length(x)
  m <- n - w + 1L
  idx <- seq_len(m)
  expanding_ar1 <- numeric(m)
  expanding_mean <- numeric(m)
  expanding_var <- numeric(m)
  expanding_skew <- numeric(m)
  expanding_kurt <- numeric(m)
  y <- x[1:w]
  s1 <- sum(y)
  s2 <- sum(y^2)
  s3 <- sum(y^3)
  s4 <- sum(y^4)
  s_cur <- s1 - y[1]
  s_lag <- s1 - x[w]
  s_lag2 <- s2 - x[w]^2
  s_prod <- sum(y[-1] * y[-w])
  mu <- s1 / w
  m2  <- (s2 - s1^2 / w) / w
  m3  <- (s3 - 3 * mu * s2 + 2 * w * mu^3) / w
  m4  <- (s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 3 * w * mu^4) / w
  mu_sq <- (w - 1) * mu^2
  expanding_ar1[1] <- ifelse_(
    demean,
    (s_prod - mu * (s_cur + s_lag) + mu_sq) / (s_lag2 - 2 * mu * s_lag + mu_sq),
    s_prod / s_lag2
  )
  expanding_mean[1] <- mu
  expanding_var[1] <- m2 * w / (w - 1)
  expanding_skew[1] <- m3 / (m2^(3/2))
  expanding_kurt[1] <- m4 / (m2^2)
  for (i in seq(2, m)) {
    w <- w + 1L
    x_new <- x[w]
    s1 <- s1 + x_new
    s2 <- s2 + x_new^2
    s3 <- s3 + x_new^3
    s4 <- s4 + x_new^4
    s_lag2 <- s_lag2 + x[w - 1L]^2
    s_prod <- s_prod + x_new * x[w - 1L]
    mu <- s1 / w
    m2  <- (s2 - s1^2 / w) / w
    m3  <- (s3 - 3 * mu * s2 + 2 * w * mu^3) / w
    m4  <- (s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 3 * w * mu^4) / w
    if (demean) {
      s_lag <- s_lag + x[w - 1L]
      s_cur <- s_cur + x_new
      mu_sq <- (w - 1) * mu^2
      expanding_ar1[i] <- (s_prod - mu * (s_cur + s_lag) + mu_sq) /
        (s_lag2 - 2 * mu * s_lag + mu_sq)
    } else {
      expanding_ar1[i] <- s_prod / s_lag2
    }
    expanding_mean[i] <- mu
    expanding_var[i] <- m2 * w / (w - 1)
    expanding_skew[i] <- m3 / (m2^(3/2))
    expanding_kurt[i] <- m4 / (m2^2)
  }
  expanding_sd <- sqrt(expanding_var)
  expanding_metrics <- data.frame(
    time = time[(burnin + 1):n],
    ar1 = expanding_ar1,
    mean = expanding_mean,
    sd = expanding_sd,
    skew = expanding_skew,
    kurt = expanding_kurt,
    cv = expanding_sd / expanding_mean,
    rr = 1.0 - expanding_ar1
  )
  expanding_metrics <- expanding_metrics[, c("time", metrics)]
  signs <- rep(1.0, ncol(expanding_metrics) - 1)
  names(signs) <- names(expanding_metrics[-1])
  if ("rr" %in% names(signs)) {
    signs["rr"] <- -1.0
  }
  long <- tidyr::pivot_longer(
    expanding_metrics,
    cols = !(!!rlang::sym("time")),
    names_to = "metric",
    values_to = "score"
  ) |>
    dplyr::group_by(!!rlang::sym("metric")) |>
    dplyr::mutate(
      z_score = expanding_z(!!rlang::sym("score")),
      crossed = !!rlang::sym("z_score") *
        signs[!!rlang::sym("metric")] > threshold,
      detected = factor(
        as.integer(crossed & dplyr::lag(crossed)),
        levels = c(1, 0)
      )
    ) |>
    dplyr::ungroup()
  structure(
    long,
    orig_values = x,
    orig_time = time,
    threshold = threshold,
    method = "expanding",
    class = c("tsn_ews", "tbl_df", "tbl", "data.frame")
  )
}

expanding_z <- function(x) {
  n <- length(x)
  var <- numeric(n)
  cent <- numeric(n)
  y <- x[1:2]
  s1 <- sum(y)
  s2 <- sum(y^2)
  mu <- s1 / 2
  m2 <- (s2 - s1^2 / 2) / 2
  cent[2] <- (x[2] - mu)
  var[2] <- m2 * 2
  for (i in 3:n) {
    x_new <- x[i]
    s1 <- s1 + x_new
    s2 <- s2 + x_new^2
    mu <- s1 / i
    m2  <- (s2 - s1^2 / i) / i
    cent[i] <- (x[i] - mu)
    var[i] <- m2 * i / (i - 1)
  }
  c(0, cent[-1] / sqrt(var[-1]))
}

detrend_ts <- function(values, time, method, window, bandwith, span, degree) {
  switch(method,
    `gaussian` = {
      smoothed <- stats::ksmooth(
        x = time,
        y = value,
        kernel = "normal",
        bandwidth = bandwidth,
        x.points = time
      )$y
      values - smoothed
    },
    `loess` = {
      fit <- stats::loess(
        values ~ time,
        span = span,
        degree = degree,
        normalize = FALSE
      )
      stats::residuals(fit)
    },
    `linear` = {
      fit <- stats::lm(values ~ time)
      stats::residuals(fit)
    },
    `first-diff` = {
      c(0, diff(values))
    },
    `none` = {
      values
    }
  )
}

#' Analyzes the expanding window EWS results to classify the system
#' state at each time point based on the number of metrics showing warnings.
#' This provides a qualitative assessment of system stability.
#'
#' @param x A `tsn_ews` object from the expanding window method.
#' @return A `data.frame` with time, warning count, and state classifications.
#' @noRd
classify_ews <- function(x) {
  n_metrics <- length(unique(x$metric))
  d <- x |>
    dplyr::filter(!is.na(!!rlang::sym("z_score"))) |>
    dplyr::group_by(!!rlang::sym("time")) |>
    dplyr::summarize(
      time = !!rlang::sym("time"),
      count = sum(as.integer(!!rlang::sym("detected")))
    )
  if (n_metrics == 1) {
   d$state <- cut(
      d$count,
      breaks = c(-Inf, 0, 1, Inf),
      labels = c("Stable", "Warning", "Critical"),
      right = TRUE
    )
  } else if (n_metrics <= 3) {
    d$state <- cut(
      d$count,
      breaks = c(-Inf, 0, 1, n_metrics, Inf),
      labels = c("Stable", "Vulnerable", "Warning", "Critical"),
      right = TRUE
    )
  } else {
    d$state <- cut(
      d$count,
      breaks = c(-Inf, 0, 1, 2, 0.5 * n_metrics, Inf),
      labels = c(
        "Stable", "Vulnerable", "Weak Warning", "Strong Warning", "Failing"
      ),
      right = TRUE
    )
  }
  d
}
