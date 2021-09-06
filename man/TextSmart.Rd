\name{TextSmart}
\alias{TextSmart}

\title{
Labels of a Scatter
}
\description{
Plots labels of points in a scattergram. labels for points with positive \code{x} are placed on the right of
the points, and labels for points with negative values on the left.
}
\usage{
TextSmart(A, Labels, CexPoints = 1, ColorPoints = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Coordinates of the points for the scattergram.
}
  \item{Labels}{
Labels for the points.
}
  \item{CexPoints}{
Size of the labels.
}
  \item{ColorPoints}{
Colors of the labels.
}
  \item{\dots}{
Additional graphical arguments.
}
}
\details{
The function is used to improve the readability of the labels in a scattergram.
}
\value{
No value returned.
}

\author{
Jose Luis Vicente-Villardon
}
