% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-sequences.R
\name{ts_f2}
\alias{ts_f2}
\alias{ts_f3}
\alias{ts_f4}
\alias{ts_f4ratio}
\title{Calculate the f2, f3, f4, and f4-ratio statistics}
\usage{
ts_f2(
  ts,
  A,
  B,
  mode = c("site", "branch", "node"),
  span_normalise = TRUE,
  windows = NULL
)

ts_f3(
  ts,
  A,
  B,
  C,
  mode = c("site", "branch", "node"),
  span_normalise = TRUE,
  windows = NULL
)

ts_f4(
  ts,
  W,
  X,
  Y,
  Z,
  mode = c("site", "branch", "node"),
  span_normalise = TRUE,
  windows = NULL
)

ts_f4ratio(
  ts,
  X,
  A,
  B,
  C,
  O,
  mode = c("site", "branch"),
  span_normalise = TRUE
)
}
\arguments{
\item{ts}{Tree sequence object of the class \code{slendr_ts}}

\item{mode}{The mode for the calculation ("sites" or "branch")}

\item{span_normalise}{Divide the result by the span of the window? Default
TRUE, see the tskit documentation for more detail.}

\item{windows}{Coordinates of breakpoints between windows. The first
coordinate (0) and the last coordinate (equal to \code{ts$sequence_length})
do not have to be specified as they are added automatically.}

\item{W, X, Y, Z, A, B, C, O}{Character vectors of individual names (following the
nomenclature of Patterson et al. 2021)}
}
\value{
Data frame with statistics calculated for the given sets of individuals
}
\description{
Calculate the f2, f3, f4, and f4-ratio statistics
}
\examples{
\dontshow{check_dependencies(python = TRUE) # make sure dependencies are present
}
# load an example model with an already simulated tree sequence
slendr_ts <- system.file("extdata/models/introgression.trees", package = "slendr")
model <- read_model(path = system.file("extdata/models/introgression", package = "slendr"))

# load the tree-sequence object from disk
ts <- ts_load(slendr_ts, model, mutate = TRUE, mutation_rate = 1e-8, random_seed = 42)

# calculate f2 for two individuals in a previously loaded tree sequence
ts_f2(ts, A = "AFR_1", B = "EUR_1")

# calculate f2 for two sets of individuals
ts_f2(ts, A = c("AFR_1", "AFR_2"), B = c("EUR_1", "EUR_3"))

# calculate f3 for two individuals in a previously loaded tree sequence
ts_f3(ts, A = "EUR_1", B = "AFR_1", C = "NEA_1")

# calculate f3 for two sets of individuals
ts_f3(ts, A = c("AFR_1", "AFR_2", "EUR_1", "EUR_2"),
          B = c("NEA_1", "NEA_2"),
          C = "CH_1")

# calculate f4 for single individuals
ts_f4(ts, W = "EUR_1", X = "AFR_1", Y = "NEA_1", Z = "CH_1")

# calculate f4 for sets of individuals
ts_f4(ts, W = c("EUR_1", "EUR_2"),
          X = c("AFR_1", "AFR_2"),
          Y = "NEA_1",
          Z = "CH_1")

# calculate f4-ratio for a given set of target individuals X
ts_f4ratio(ts, X = c("EUR_1", "EUR_2", "EUR_4", "EUR_5"),
               A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CH_1")
}
