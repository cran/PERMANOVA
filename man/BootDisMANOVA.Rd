\name{BootDisMANOVA}
\alias{BootDisMANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Analysis of Variance based on Distances and Bootstrap
}
\description{
Multivariate Analysis of variance based on distances and Bootstrap.
}
\usage{
BootDisMANOVA(Distance, groups, C = NULL, Effects = NULL, nB = 1000, seed = NULL,
CoordPrinc = FALSE, dimens = 2, PCoA = "Standard", ProjectInd = TRUE, tol = 1e-04,
DatosIni = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Distance}{
A matrix of distances.
}
  \item{groups}{
A factor containing the groups to compare.
}
  \item{C}{
A matrix of contrasts (if null the identity is used).
}
  \item{Effects}{
A vector of effects.
}
  \item{nB}{
Number of Bootstrap replicates.
}
  \item{seed}{
Seed for the random numbers.
}
  \item{CoordPrinc}{
Should Principal Coordinates be calculated.
}
  \item{dimens}{
Dimension of the solution.
}
  \item{PCoA}{
Type of Principal Coordinates to calulate.
}
  \item{ProjectInd}{
Should the individuals be projected onto the graph.
}
  \item{tol}{
Tolerance for convergence of the algorithms.
}
  \item{DatosIni}{
Should the initial data be included in the results.
}
}
\details{
Multivariate Analysis of Variance based on distances and Bootstrap.
}

\value{
\item{call}{Function}
\item{Title}{Title of the study}
\item{Type}{BootMANOVA}
\item{Distances}{A matrix containing the distances between individuals.}
\item{C}{Contrasts Matrix.}
\item{Initial}{ Containing two matrices:

    * Global -> Global contrast.

    * Contrastes ->Contrasts for groups.
}
\item{DistMuestral}{Sample distribution of F-exp from permutations.}
\item{pvalue}{Estimate p-valor for PERMANOVA.}
\item{ExplainedVariance}{Explained variance by Principal Coordinates selected.}
\item{Inertias}{Eigenvalue, Explained variance, Cumulative explained variance.}
\item{MeanCoordinates}{Mean Coordinates by groups for the dimensions obtained in the Principal Coordinates Analysis.}
\item{Qualities}{Qualities representation by groups for the dimensions of PCoA.}
\item{CummulativeQualities}{Cummulative qualities representation.}
\item{ClusterType}{Cluster type selected.}
\item{Clusters}{Clusters created.}
\item{ClusterNames}{Names of clusters}
\item{ClusterColors}{Colors of clusters, color name and HTML code.}

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
D = DistContinuous (X)
bootwine=BootDisMANOVA(D, wine$Group)
bootwine
}
}
