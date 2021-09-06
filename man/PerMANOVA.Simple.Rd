\name{PerMANOVA.Simple}
\alias{PerMANOVA.Simple}

\title{
PERMANOVA from a matrix of distancies
}
\description{
This function makes a PERMANOVA from distancies matrix and factor for groups.
}
\usage{
PerMANOVA.Simple(D, grupo, nperm = 999, seed = NULL, C = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
A matrix containing the distances between individuals.
}
  \item{grupo}{
A factor containing the groups to compare.
}
  \item{nperm}{
Number of permutation that want to perform. By default it is 999.
}
  \item{seed}{
Seed to start permutations. By default it is null.
}
  \item{C}{
Contrast matrix. By default it is null.
}
}

\value{
The PERMANOVA.Simple function create a list that return:
\item{call}{Function}
\item{nperm}{Number of permutation.}
\item{Inicial}{Containing:

* TSS -> Total sum of squares.

* BSS -> Between groups sum of squares.

* WSS -> Within groups sum of squares.

* glt -> Total degrees of freedom.

* glb -> Between groups degrees of freedom.

* glw -> Within groups degrees of freedom.

* Fexp -> Experimental F.
}
\item{Fvals}{F values of the permutation.}
\item{pval}{Estimate p-valor for PERMANOVA.}
}
\references{
Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46.

Anderson, M. J. (2005). Permanova: a fortran computer program for permutational multivariate analysis of variance. Department of Statistics, University of Auckland, New Zealand, 24.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
\donttest{data(wine)
X = wine[,4:21]
Dist = DistContinuous (X)
PerMANOVA.Simple(Dist$D, wine$Group)
  }
}
