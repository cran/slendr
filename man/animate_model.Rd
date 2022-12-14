% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{animate_model}
\alias{animate_model}
\title{Animate the simulated population dynamics}
\usage{
animate_model(model, file, steps, gif = NULL, width = 800, height = 560)
}
\arguments{
\item{model}{Compiled \code{slendr_model} model object}

\item{file}{Path to the table of saved individual locations}

\item{steps}{How many frames should the animation have?}

\item{gif}{Path to an output GIF file (animation object returned
by default)}

\item{width, height}{Dimensions of the animation in pixels}
}
\value{
If \code{gif = NULL}, return gganimate animation object. Otherwise a GIF
file is saved and no value is returned.
}
\description{
Animate the simulated population dynamics
}
