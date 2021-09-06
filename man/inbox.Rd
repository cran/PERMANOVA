\name{inbox}
\alias{inbox}

\title{
Checks if a point is inside a box
}
\description{
Checks if a point is inside a box. The point is specified bi its x and y coordinates and the box with the minimum and maximum values on both coordinate axis: xmin, xmax, ymin, ymax. The vertices of the box are then (xmin, ymin), (xmax, ymin), (xmax, ymax) and (xmin, ymax).
}
\usage{
inbox(x, y, xmin, xmax, ymin, ymax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x coordinate of the point.
}
  \item{y}{
y coordinate of the point.
}
  \item{xmin}{
Minimum value of X.
}
  \item{xmax}{
Maximum value of X.
}
  \item{ymin}{
Minimum value of Y.
}
  \item{ymax}{
Maximum value of Y.
}
}

\value{
Returns a logical value: TRUE if the point is inside the box and FALSE otherwise.
}

\author{
Jose Luis Vicente Villardon
}

\examples{
inbox(0, 0, -1, 1, -1, 1)
}
