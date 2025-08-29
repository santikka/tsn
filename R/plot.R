#' Plot Time-Series Data with State Frequencies
#'
#' @export
#' @param data \[`tsn`, `ts`, `data.frame`, `numeric()`]\cr
#'   Time-series data to be plotted.
#' @param selected \[`character()`]\cr A vector of indices or names of
#'   individual time-series to plot. If not provided (default), all
#'   time-series are plotted up to `max_series` number of plots.
#' @param overlay \[`logical(1)`]\cr An option for plotting the overlay
#'   that indicates the state assigned to each time point. Can be either `"h"`
#'   for a horizontal overlay or `"v"` for a vertical overlay (default). If
#'   `NULL`, no overlay is plotted.
#' @param points \[`logical(1)`]\cr Should a point be added for each
#'   observation? Defaults to `FALSE`. The points are colored according to the
#'   assigned state.
#' @param ncol An `integer` giving the number of columns to use for the facets.
#' @param max_series \[`integer(1)`]\cr The maximum number  of time-series to
#'   plot. The default is 10.
#' @param trend \[`logical(1)`]\cr Should are trend line be added to the plot?
#'   The default is `FALSE` for no trend line.
#' @param scales \[`character(1)`]\cr Any of `"fixed"`, `"free_x"`, `"free_y"`,
#'   or `"free"` (default).
#' @return A `ggplot` object.
#' @examples
#' ts_data <- data.frame(
#'   id = gl(10, 100),
#'   time = rep(1:100,10),
#'   series = c(
#'     replicate(
#'       10,
#'       stats::arima.sim(list(order = c(2, 1, 0), ar = c(0.5, 0.2)), n = 99)
#'     )
#'   )
#' )
#'
#' ts_data_disc <- discretize(ts_data, "id", "series", "time", n_states = 5)
#' plot_series(ts_data_disc)
#'
plot_series <- function(data, selected, overlay = "v", points = FALSE,
                        ncol = NULL, max_series = 10, trend = FALSE,
                        scales = c("free", "free_x", "free_y", "fixed")) {
  data <- as.tsn(data)
  if (!is.null(overlay)) {
    overlay <- check_match(overlay, c("h", "v"))
  }
  check_flag(points)
  id_col <- attr(data, "id_col")
  value_col <- attr(data, "value_col")
  state_col <- attr(data, "state_col")
  time_col <- attr(data, "time_col")
  ids <- unique(data[[id_col]])
  selected <- ifelse_(
    missing(selected),
    ids[seq_len(min(length(ids), max_series))],
    selected[seq_len(min(length(selected), max_series))]
  )
  data <- data[data[[id_col]] %in% selected, ]
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = !!rlang::sym(time_col), y = !!rlang::sym(value_col))
  )
  # Create segments where the state remains constant
  if (!is.null(overlay) && !is.null(state_col)) {
    if (overlay == "v") {
      segment_col <- time_col
      xmin <- rlang::sym(".min")
      xmax <- rlang::sym(".max")
      ymin <- rlang::sym(".neginf")
      ymax <- rlang::sym(".posinf")
    } else if (overlay == "h") {
      segment_col <- value_col
      xmin <- rlang::sym(".neginf")
      xmax <- rlang::sym(".posinf")
      ymin <- rlang::sym(".min")
      ymax <- rlang::sym(".max")
    }
    rects <- data |>
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(segment_col)) |>
      dplyr::group_by(!!rlang::sym(id_col)) |>
      dplyr::mutate(
        .grouping_var = cumsum(
          !!rlang::sym(state_col) != dplyr::lag(
            !!rlang::sym(state_col),
            default = dplyr::first(!!rlang::sym(state_col))
          )
        ),
        .lag = dplyr::lag(
          !!rlang::sym(segment_col),
          default = dplyr::first(!!rlang::sym(segment_col))
        ),
        .lead = dplyr::lead(
          !!rlang::sym(segment_col),
          default = dplyr::last(!!rlang::sym(segment_col))
        )
      ) |>
      dplyr::group_by(
        !!rlang::sym(id_col),
        .grouping_var,
        !!rlang::sym(state_col)
      ) |>
      dplyr::summarise(
        .neginf = -Inf,
        .posinf = Inf,
        .min = 0.5 * (min(!!rlang::sym(segment_col)) + min(.lag)),
        .max = 0.5 * (max(!!rlang::sym(segment_col)) + max(.lead)),
        .groups = "drop"
      )
    p <- p + ggplot2::geom_rect(
      data = rects,
      ggplot2::aes(
        xmin = !!xmin,
        xmax = !!xmax,
        ymin = !!ymin,
        ymax = !!ymax,
        fill = !!rlang::sym(state_col)
      ),
      alpha = 0.5,
      inherit.aes = FALSE
    )
  }
  # TODO fill
  p <- p + ggplot2::geom_line(linewidth = .5)
  if (points) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = !!rlang::sym(state_col)),
      show.legend = FALSE,
      pch = 21
    )
  }
  if (trend) {
    p <- p + ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "darkgray",
      lty = 2
    )
  }
  if (length(selected) > 1) {
    p <- p + ggplot2::facet_wrap(
      id_col,
      ncol = ncol,
      scales = scales
    )
  }
  p +
    ggplot2::scale_fill_brewer(
      palette = ifelse(
        n_unique(data[[state_col]]) <= 8,
        "Accent",
        "Set3"
      ),
      limits = levels(
        factor(base::sort(dplyr::pull(data[, state_col, drop = FALSE], 1L)))
      ),
      name = "State"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "") +
    ggplot2::theme(legend.position = "bottom")
}

#' @export
plot.tsn_ews <- function(x, ...) {
  d <- data.frame(
    value = attr(x, "orig_values"),
    time = attr(x, "orig_time")
  )
  p_ts <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  if (attr(x, "method") == "rolling") {
    p_ts <- p_ts +
      ggplot2::labs(title = "Time Series", y = "Value", x = NULL)
    p_metrics <- plot_rolling_ews(x, ...)
    #patchwork::wrap_plots(p_ts, p_metrics, ncol = 1L) +
    #  patchwork::plot_layout(guides = "collect")
  } else {
    warn <- data.frame(time = x$time[x$detected == 1])
    p_ts <- p_ts +
      ggplot2::geom_rug(
        data = warn,
        ggplot2::aes(
          x = !!rlang::sym("time"), color = "EWS"
        ),
        show.legend = TRUE,
        sides = "b",
        length = ggplot2::unit(0.05, "npc"),
        inherit.aes = FALSE
      ) +
      # Add point to control shape
      ggplot2::geom_point(
        data = data.frame(x = -Inf, y = -Inf),
        ggplot2::aes(
          x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = "Void"
        ),
        show.legend = TRUE,
        na.rm = TRUE
      ) +
      ggplot2::labs(
        title = "Time Series with Detected Warnings",
        y = "Value",
        x = NULL
      ) +
      ggplot2::scale_color_manual(
        name = "EWS",
        values = c(EWS = "red", Void = "white"),
        labels = c("Detected", "")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(
            linetype = c(0, 0),
            shape = c("|", "."),
            size = 4
          )
        )
      )

    p_metrics <- plot_expanding_ews(x,  ...)
    #p_cls <- plot_classification(x)
  }
  #patchwork::wrap_plots(p_ts, p_metrics, p_cls, ncol = 1L)
  patchwork::wrap_plots(p_ts, p_metrics, ncol = 1L)
}

plot_rolling_ews <- function(x, ...) {
  cor <- attr(x, "cor")
  x$metric <- factor(
    x$metric,
    levels = names(cor),
    labels = lapply(
      paste0(
        names(cor), "~(tau == ", round(cor, 2), ")"
      ),
      str2lang
    )
  )
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("std"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(
      as.formula("~metric"),
      scales = "free_y",
      ncol = 3,
      drop = TRUE,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      subtitle = "Rolling EWS Indicators",
      y = "Metric Value",
      x = "Time"
    ) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(
        face = "bold", size = 10, margin = ggplot2::margin(b = 5)
      ),
      panel.spacing = ggplot2::unit(1, "lines"),
      axis.text = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(angle = 90, size = 10),
      plot.subtitle = ggplot2::element_text(
        size = 11, margin = ggplot2::margin(b = 10)
      ),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
    )
}

plot_expanding_ews <- function(x, ...) {
  warn <- x[x$detected == 1, ]
  cls <- classify_ews(x)
  state_colors <- c(
    "Stable" = "#440154FF",
    "Vulnerable" = "#3B528BFF",
    "Weak Warning" = "#21908CFF",
    "Strong Warning" = "#5DC863FF",
    "Failing" = "#FDE725FF",
    "Warning" = "orange",
    "Critical" = "red"
  )
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("z_score"),
      color = !!rlang::sym("metric"),
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      data = warn,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("z_score"),
        alpha = "EWS"
      ),
      size = 2.5
    ) +
    # Add line for alpha
    ggplot2::geom_line(
     data = warn,
     ggplot2::aes(
       x = !!rlang::sym("time"),
       y = !!rlang::sym("z_score"),
       color = !!rlang::sym("metric"),
       alpha = "Void"
     ),
     inherit.aes = FALSE,
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(-4, -2, 0, 2, 4)
    ) +
    ggplot2::labs(
      title = "Strength of Early Warning Signals",
      #subtitle = "Points indicate when a metric's strength crosses the threshold",
      y = "Scaled Metric Value",
      x = NULL
    ) +
    ggplot2::scale_color_viridis_d(
      name = "EWS Indicator",
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          linewidth = 1,
          linetype = "solid",
          shape = NA
        )
      )
    ) +
    ggplot2::scale_alpha_manual(
      name = "EWS",
      values = c(EWS = 1, Void = 1),
      labels = c("Detected", "Not Detected"),
      guide = ggplot2::guide_legend(
        order = 1,
        alpha = c(1, 1),
        linetype = c("blank", "solid")
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::geom_hline(
      yintercept = c(1, -1) * attr(x, "threshold"),
      linetype = "dashed",
      color = "grey50"
    )
}

plot_classification <- function(x) {
  state_colors <- c(
    "Stable" = "#440154FF",
    "Vulnerable" = "#3B528BFF",
    "Weak Warning" = "#21908CFF",
    "Strong Warning" = "#5DC863FF",
    "Failing" = "#FDE725FF",
    "Warning" = "orange",
    "Critical" = "red"
  )
  x$state <- factor(x$state, levels = names(state_colors))
  ggplot2::ggplot(x, ggplot2::aes(x = !!rlang::sym("time"), y = 1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = !!rlang::sym("state"), height = 1)) +
    ggplot2::scale_fill_manual(
      values = state_colors,
      name = "System State",
      drop = FALSE
    ) +
    ggplot2::labs(x = "Time Point", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
