\name{MANOVA.Estimation}
\alias{MANOVA.Estimation}

\title{
Estimation of the MANOVA parameters.
}
\description{
The function estimates the parameters of the MANOVA for a set of contrasts and a set of linear combinations of the dependent variables.
}
\usage{
MANOVA.Estimation(Y, X, C, M)
}

\arguments{
  \item{Y}{
The matrix of dependent variables.
}
  \item{X}{
The matrix of independent variables.
}
  \item{C}{
A matrix containing a set of contrasts on the columns of X.
}
  \item{M}{
A matrix containing a set of linear combinations of the columns of Y.
}
}
\details{
The function estimates the parameters of the MANOVA for a set of contrasts contained in a matrix C and a set of linear combinations of the dependent variables in a matrix M. The  function is actually used in a more general procedure called MANOVA.
}
\value{
Returns the main results of the MANOVA.
}
\references{
Seber, G. A. (2009). Multivariate observations (Vol. 252). John Wiley & Sons.
}
\author{
Jose Luis Vicente-Villardon, Laura Vicente-Gonzalez
}
