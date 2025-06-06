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

