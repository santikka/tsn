% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexity.R
\name{complexity}
\alias{complexity}
\title{Calculate Dynamic Complexity Measures for Time-Series Data}
\usage{
complexity(data, measures = "complexity", window = 7L, align = "center")
}
\arguments{
\item{data}{[\code{tsn}, \code{ts}, \code{data.frame}, \code{numeric()}]\cr Time-series data.}

\item{measures}{[\code{character()}]\cr A vector of measures to calculate. The
available options are: \code{"complexity"}, \code{"fluctuation"}, \code{"distribution"},
\code{"autocorrelation"}, \code{"max"}, \code{"min"}, \code{"variance"}, \code{"all"}.
The default is "complexity". See 'Details' for more information on these
measures.}

\item{window}{[\code{integer(1)}]\cr A positive integer specifying the rolling
window size. Must be at least 2. The default is 7.}

\item{align}{[\code{character(1)}]\cr Alignment of the window. The available
options are: \code{"center"} (default), \code{"right"}, and \code{"left"}. The calculated
measure is assigned to the center, rightmost, or leftmost point of the
window, respectively.}
}
\value{
A \code{data.frame} with the time index, the original time-series data,
and the calculated measures.
}
\description{
Computes dynamic complexity and other rolling window measures for univariate
or multivariate time series data. Supports various measures including
complexity, fluctuation, distribution, autocorrelation, and basic statistics.
}
\details{
The following measures can be calculated:
\itemize{
\item \code{complexity}: Product of fluctuation and distribution measures.
\item \code{fluctuation}: Root mean square of successive differences.
\item \code{distribution}: Deviation from uniform distribution.
\item \code{autocorrelation}: Lag-1 autocorrelation coefficient.
\item \code{max}: Rolling maximum.
\item \code{min}: Rolling minimum.
\item \code{variance}: Rolling variance.
}

The option \code{"all"} computes all of the above.
}
\examples{
# Basic complexity calculation
set.seed(123)
ts_data <- rnorm(100)
result <- complexity(ts_data, measures = "complexity")

# Multiple measures
result_multi <- complexity(ts_data, measures = c("complexity", "variance"))

}
