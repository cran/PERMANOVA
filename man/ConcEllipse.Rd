\name{ConcEllipse}
\alias{ConcEllipse}

\title{
Non-parametric concentration ellipses
}
\description{
Calculates non-parametric concentration ellipses for a two-dimensional set of points
}
\usage{
ConcEllipse(data, center = NULL, confidence = 1, npoints = 100)
}

\arguments{
  \item{data}{
A two dimensional set of points.
}
  \item{center}{
Center of the data.
}
  \item{confidence}{
Confidence for the concentration ellipse.
}
  \item{npoints}{
Number of points to plot.
}
}
\details{
Calculates non-parametric concentration ellipses for a two-dimensional set of points using distances to the center.
}
\value{
\item{data}{Original data}
\item{confidence}{Percent of points selected for the ellipse}
\item{ellipse}{points to plot the ellipse}
\item{center}{center of the ellipse}
}
\references{
Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories.

Linting, M., Meulman, J. J., Groenen, P. J., & Van der Kooij, A. J. (2007). Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychological Methods, 12(3), 359.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
x=runif(30)
y=x+rnorm(30)*0.3
plot(x, y)
el=ConcEllipse(cbind(x,y), confidence=0.95)
plot(el)
}
