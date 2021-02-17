\name{Circle2}
\alias{Circle2}

\title{
Draws a circle.
}
\description{
Draws a circle for a given radius at the specified center with the given color.
}
\usage{
Circle2(radius = 1, origin = c(0, 0), col = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{radius}{
Radius of the circle.
}
  \item{origin}{
Centre of the circle.
}
  \item{col}{
Color of the circle.
}
  \item{\dots}{
Additional graphical parameters.
}
}
\details{
Draws a circle for a given radius at the specified center with the given color.
}
\value{
No value is returned.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
plot(0,0)
Circle2(1,c(0,0))
}
