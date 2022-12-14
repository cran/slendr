% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{check_env}
\alias{check_env}
\title{Check that the active Python environment is setup for slendr}
\usage{
check_env(quiet = FALSE)
}
\arguments{
\item{quiet}{Should a log message be printed? If \code{FALSE}, only a logical
value is returned (invisibly).}
}
\value{
Either \code{TRUE} (slendr Python environment is present) or \code{FALSE}
(slendr Python environment is not present).
}
\description{
This function inspects the Python environment which has been activated by the
reticulate package and prints the versions of all slendr Python dependencies
to the console.
}
\examples{
\dontshow{check_dependencies(python = TRUE) # make sure dependencies are present
}
check_env()
}
