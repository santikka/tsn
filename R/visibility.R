#' Compute a Visibility Graph from Time-Series Data
#'
#' Creates visibility graphs from time series data, where nodes represent time
#' points and edges represent visibility between points according to the
#' selected method. The function supports both Natural Visibility Graphs (NVG)
#' and Horizontal Visibility Graphs (HVG), with options for limited penetrable
#' visibility and temporal decay.
#'
#' @export
#' @param data \[`tsn`, `ts`, `data.frame`, `numeric()`]\cr Time-series data.
#' @param method \[`character(1)`]\cr The visibility graph construction
#'   method. The options are `"nvg"` for natural visibility graphs and
#'   `"hvg"` for horizontal visibility graphs.
#' @param directed \[`logical(1)`]\cr Should the graph be a directed graph?
#'   If `TRUE`, edges have direction showing temporal order.
#'   The default is `FALSE`.
#' @param limit \[`integer(1)`]\cr Maximum temporal distance (in time steps)
#'   for visibility connections. If not provided (default), no limit is applied.
#' @param penetrable \[`integer(1)`]\cr Number of points allowed to penetrate
#'   the visibility line for Limited Penetrable Visibility Graphs (LPVG).
#'   The default is 0 (standard visibility graph).
#' @param decay_factor \[`numeric(1)`]\cr Temporal decay factor of
#'   edge weights. The default is 0 (no decay). Higher values cause faster
#'   decay with increasing temporal distance.
#' @return A `matrix` of edge weights.
#'
#' @details
#' ## Visibility Graph Methods
#' **Natural Visibility Graph (NVG):** Two points can "see" each other if
#' a straight line connecting them doesn't intersect with any intermediate
#' points. This preserves more geometric information from the original time
#' series.
#'
#' **Horizontal Visibility Graph (HVG):**  Two points can "see" each other if
#' all intermediate points are strictly lower than both endpoints. This creates
#' sparser graphs that capture different temporal patterns.
#'
#' ## Advanced Features
#' **Limited Penetrable Visibility (LPVG):** Allows a specified number of
#' points to penetrate the visibility line, creating more connected graphs.
#'
#' **Temporal Decay:** Edge weights decrease exponentially with temporal
#' distance, emphasizing local temporal relationships.
#'
#' @references
#' Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C. (2008).
#' From time series to complex networks: The visibility graph.
#' *PNAS*, **105(13)**, 4972-4975.
#'
#' Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009).
#' Horizontal visibility graphs: Exact results for random time series.
#' *Physical Review E*, **80(4)**, 046103.
#'
#' @examples
#' set.seed(123)
#' ts_data <- rnorm(50, mean = 50, sd = 15)
#'
#' # Create natural visibility graph
#' nvg <- visibility_graph(ts_data, method = "nvg")
#'
#' # Create horizontal visibility graph with temporal limit
#' hvg <- visibility_graph(ts_data, method = "hvg", limit = 10)
#'
visibility_graph <- function(data, method = "nvg", directed = FALSE, limit,
                             penetrable = 0, decay_factor = 0) {
  data <- as.tsn(data)
  method <- check_match(method, c("nvg", "hvg"))
  n <- nrow(data)
  limit <- ifelse_(missing(limit), n, limit)
  check_flag(directed)
  check_values(limit) # TODO values check
  check_values(penetrable)
  check_values(decay_factor)
  limit <- as.integer(limit)
  penetrable <- as.integer(penetrable)
  edges <- visibility_edges(
    values = get_values(data),
    time = get_time(data),
    method = method,
    directed = directed,
    limit = limit,
    penetrable = penetrable,
    decay_factor = decay_factor
  )
  edges
}

#' Visibility Edges
#'
#' @param values A `numeric` vector of the time-series data.
#' @param time A `numeric` vector of the time points.
#' @inheritParams visibility_graph
#' @return A weighted adjacency `matrix` of the network.
#' @noRd
visibility_edges <- function(values, time, method, directed, limit,
                             penetrable, decay_factor) {
  n <- length(values)
  edges <- matrix(0.0, nrow = n, ncol = n)
  visibility_fun <- ifelse_(
    method == "nvg",
    natural_visibility,
    horizontal_visibility
  )
  for (i in 1L:(n - 1L)) {
    for (j in (i + 1L):(min(n, i + limit))) {
      if (visibility_fun(x = values, y = time, i, j, penetrable)) {
        weight <- exp(-decay_factor * abs(j - i))
        edges[i, j] <- weight
      }
    }
  }
  if (!directed) {
    edges <- edges + t(edges)
  }
  edges
}

#' Check natural visibility between two points
#'
#' Natural visibility between two points exists if you can draw a straight line
#' between them that doesn't intersect with any intermediate points
#' (or intersects with at most `penetrable` points for LPVG variants).
#'
#' @param x A `numeric` vector of the time series data.
#' @param y A `numeric` vector of the time values.
#' @param i An `integer` index of the first time point.
#' @param j An `integer` index of the second time point.
#' @inheritParams visibility_graph.
#' @return A `logical` value indicating whether the endpoints of `x` have
#' natural visibility
#' @noRd
natural_visibility <- function(x, y, i, j, penetrable) {
  if (abs(i - j) == 1L) {
    return(TRUE)
  }
  penetrations <- 0L
  idx <- (i + 1L):(j - 1L)
  x_start <- x[i]
  x_end <- x[j]
  y_start <- y[i]
  y_end <- y[j]
  b <- (x_start - x_end) / (y_end - y_start)
  for (k in idx) {
    line_height <- x_end + b * (y_end - y[k])
    if (x[k] >= line_height) {
      penetrations <- penetrations + 1L
      if (penetrations > penetrable) {
        return(FALSE)
      }
    }
  }
  TRUE
}

#' Check horizontal visibility between two points
#'
#' Horizontal visibility between two points exists if all intermediate points
#' are strictly lower than both endpoints (or at most `penetrable` points
#' violate this condition for LPVG variants).
#'
#' @param x A `numeric` vector of the time series data.
#' @param y Not used.
#' @inheritParams visibility_graph.
#' @return A `logical` value indicating whether the endpoints of `x` have
#' horizontal visibility
#' @noRd
horizontal_visibility <- function(x, y, i, j, penetrable) {
  if (abs(i - j) == 1L) {
    return(TRUE)
  }
  x_min <- min(x[c(i, j)], na.rm = TRUE)
  sum(x[(i + 1L):(j - 1L)] > x_min) <= penetrable
}
