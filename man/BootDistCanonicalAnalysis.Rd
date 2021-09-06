\name{BootDistCanonicalAnalysis}
\alias{BootDistCanonicalAnalysis}

\title{Canonical Analysis based on Distances
}
\description{
Canonical Analysis based on distances. Confidence Regions for the mean vectors are calculated using bootstrap.
}
\usage{
BootDistCanonicalAnalysis(Distance, groups, dimens = NULL, nB = 100, seed = NULL,
PCoA = "Standard", ProcrustesRot = TRUE, DatosIni = TRUE, tol = 1e-04)
}

\arguments{
  \item{Distance}{
A list of three elements containing the data, the distances between individuals and type of distance used.
}
  \item{groups}{
A factor containing the groups to compare.
}
  \item{dimens}{
Number of dimensions to choose . By default is 2.
}
  \item{nB}{
Number of Bootstrap samples.
}
  \item{seed}{
Seed for the generation of the random samples (Added for reproducibility)
}
  \item{PCoA}{
The type of Principal Coordinates Analysis.

  There are two possibilities:

   * Standard

   * Weighted

  By default is Standard.
}
  \item{ProcrustesRot}{
Should Procrustes rotation be applied to the configurations?
}
  \item{DatosIni}{
The input object contains the initial data?. By default is TRUE.
}
  \item{tol}{
Tolerance
}
}
\details{
Calculates a Canonical Analysis based on distance matrices with confidence regions based on bootstrap resampling.
}
\value{
\item{D}{The distance matrix used for calculations}
\item{Coefficient }{Dissimilarity  coefficient used to calculate the proximities}
\item{nB}{Number of bootstrap samples}
\item{Groups}{Factor containing the groups to compare}
\item{GroupNames}{Names of the groups}
\item{Inertia}{Variance accounted for ech canonical coordinate}
\item{Inertias}{able with the eigenvalue, variance accounted and cumulated variance}
\item{MeanCoordinates}{Coordinates of the groups in the representation}
\item{Qualities}{Qualities of the representation of the group means}
\item{CumulativeQualities}{Cumulative Qualities of the representation}
\item{CoordBoot}{Coordinates of the bootstrap replicates}
\item{RowCoordinates}{Coordinates of the individuals}
}
\references{
Gower, J. C. and Krzanowski, W. J. (1999). Analysis of distance for structured multivariate
data and extensions to multivariate analysis of variance. Journal of the Royal
Statistical Society: Series C (Applied Statistics), 48(4):505–519.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
\donttest{data(wine)
X=wine[,4:21]
Dist=DistContinuous(X)
canon=BootDistCanonicalAnalysis(Dist, groups=wine$Group, nB=10)
}}

