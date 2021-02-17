\name{PERMANOVA.Estimation}
\alias{PERMANOVA.Estimation}

\title{
Estimation of the PERMANOVA parameters
}
\description{
Estimation of the PERMANOVA parameters.
}
\usage{
PERMANOVA.Estimation(D, X, C, Effects = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
A matrix of distances.
}
  \item{X}{
A matrix of independent variables.
}
  \item{C}{
Contrast matrix. By default it is null and the identity is used.
}
  \item{Effects}{
A factor with the effects in the rows of C. By default it is null and each row is considered as an effect.
}
}
\details{
Estimation of the PERMANOVA parameters.
}
\value{
\item{Global}{Significance for the global effect}
\item{Contrastes}{Significance for each contrast}
}
\references{
Anderson, M. J. (2008). A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46.

Anderson, M. J. (2005). Permanova: a fortran computer program for permutational multivariate analysis of variance. Department of Statistics, University of Auckland, New Zealand, 24.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
