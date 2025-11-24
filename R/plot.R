#' Plot Time-Series Data with State Frequencies
#'
#' @export
#' @param x \[`tsn`]\cr Time-series data to be plotted.
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
#' @param max_series \[`integer(1)`]\cr The maximum number of time-series to
#'   plot. The default is 10.
#' @param trend \[`logical(1)`]\cr Should trend lines be added to the plot?
#'   The default is `FALSE` for no trend lines.
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
#' plot(ts_data_disc)
#'
plot.tsn <- function(x, selected, overlay = "v", points = FALSE,
                     ncol = NULL, max_series = 10, trend = FALSE,
                     scales = c("free", "free_x", "free_y", "fixed")) {
  check_missing(x)
  check_class(x, "tsn")
  if (!is.null(overlay)) {
    overlay <- check_match(overlay, c("h", "v"))
  }
  check_flag(points)
  ids <- unique(x$id)
  selected <- ifelse_(
    missing(selected),
    ids[seq_len(min(length(ids), max_series))],
    selected[seq_len(min(length(selected), max_series))]
  )
  states <- factor(base::sort(dplyr::pull(x[, "state", drop = FALSE], 1L)))
  data <- x[x$id %in% selected, ]
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
  )
  # Create segments where the state remains constant
  if (!is.null(overlay) && !is.null(x$state)) {
    if (overlay == "v") {
      segment_col <- "time"
      xmin <- rlang::sym(".min")
      xmax <- rlang::sym(".max")
      ymin <- rlang::sym(".neginf")
      ymax <- rlang::sym(".posinf")
    } else if (overlay == "h") {
      segment_col <- "value"
      xmin <- rlang::sym(".neginf")
      xmax <- rlang::sym(".posinf")
      ymin <- rlang::sym(".min")
      ymax <- rlang::sym(".max")
    }
    rects <- data |>
      dplyr::arrange(!!rlang::sym("id"), !!rlang::sym(segment_col)) |>
      dplyr::group_by(!!rlang::sym("id")) |>
      dplyr::mutate(
        .grouping_var = cumsum(
          !!rlang::sym("state") != dplyr::lag(
            !!rlang::sym("state"),
            default = dplyr::first(!!rlang::sym("state"))
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
        !!rlang::sym("id"), .grouping_var, !!rlang::sym("state")
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
        fill = !!rlang::sym("state")
      ),
      alpha = 0.5,
      show.legend = TRUE,
      inherit.aes = FALSE
    )
  }
  # TODO fill
  p <- p + ggplot2::geom_line(linewidth = .5)
  if (points) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = !!rlang::sym("state")),
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
      "id",
      ncol = ncol,
      scales = scales
    )
  }
  p +
    ggplot2::scale_fill_brewer(
      palette = ifelse(
        n_unique(states) <= 8,
        "Accent",
        "Set3"
      ),
      limits = levels(states),
      name = "State",
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "Value") +
    ggplot2::theme(legend.position = "bottom")
}
