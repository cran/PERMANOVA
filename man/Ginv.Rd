\name{Ginv}
\alias{Ginv}

\title{
G inverse
}
\description{
Calculates the g-inverse of a squared matrix using the eigen decomposition and removing the eigenvalues smaller than a tolerance.
}
\usage{
Ginv(X, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{X}{
Matrix to calculate the g-inverse.
}
  \item{tol}{
Tolerance.
}
}
\details{
The function is useful to avoid singularities.
}
\value{
Returns the g-inverse.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
data(iris)
x=as.matrix(iris[,1:4])
S= t(x) %*% x
Ginv(S)
}
