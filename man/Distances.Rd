\name{Distances}
\alias{Distances}
\title{
Unified interface for computing distance matrices
}
\description{
Computes a distance matrix of class \code{"proximities"} by automatically detecting whether the input variables are continuous, binary, categorical, or mixed-type, or by dispatching to the specified distance function.
}
\usage{
Distances(x, y = NULL,
  type = c("auto", "continuous", "binary", "categorical", "mixed"),
  coef = NULL, ...)
}
\arguments{
  \item{x}{
A matrix or data frame of observations (rows) by variables (columns).
}
  \item{y}{
Optional second matrix of observations (only supported for \code{"continuous"} and \code{"binary"} types).
}
  \item{type}{
Type of data: \code{"auto"} (default), \code{"continuous"}, \code{"binary"}, \code{"categorical"}, or \code{"mixed"}.
}
  \item{coef}{
Distance coefficient passed to the underlying distance function (\code{\link{DistContinuous}}, \code{\link{DistBinary}}, \code{\link{DistCategorical}}, or \code{\link{DistMixed}}). If \code{NULL}, the default coefficient for the detected data type is used.
}
  \item{...}{
Additional arguments passed to the underlying distance function.
}
}
\value{
An object of class \code{"proximities"} ready to be used as input in \code{\link{PERMANOVA}}.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
data(wine)
D <- Distances(wine[, 4:21])
}
