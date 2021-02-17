\name{PERMANOVA}
\alias{PERMANOVA}

\title{
PERMANOVA: MANOVA based on distances
}
\description{
The correct application of MANOVA needs normal and homocedastic data and the number of variables be much smaller than the number of individuals, but for many applications the conditions do not hold. To extend the application to this data Anderson develops PERMANOVA. This non-parametric test based on distances uses permutation to approximate the sampling distribution of the test statistic.

Contrasts and Effects can be added to the calculations.
}
\usage{
PERMANOVA(Distance, group, C = NULL, Effects = NULL, nperm = 1000, seed = NULL,
CoordPrinc = FALSE, dimens = 2, PCoA = "Standard", ProjectInd = TRUE, tol = 1e-04,
DatosIni = TRUE)
}

\arguments{
  \item{Distance}{
A list of three elements containing the data, the distances between individuals and type of distance used.
}
  \item{group}{
A factor containing the groups to compare.
}
  \item{C}{
Contrast matrix. By default it is null and the identity is used.
}
  \item{Effects}{
A factor with the effects in the rows of C. By default it is null and each row is considered as an effect.
}
  \item{nperm}{
Number of permutations to perform. By default is 1000.
}
  \item{seed}{
Seed to start permutations. By default is null.
}
  \item{CoordPrinc}{
Should the principal coordinates be calculated?. By default is FALSE.
}
  \item{dimens}{
Number of dimensions to choose in the PCoA. By default is 2.
}
  \item{PCoA}{
The type of Principal Coordinates Analysis.

  There are two possibilities:

   * Standard

   * Weighted

  By default is Standard.
}
  \item{ProjectInd}{
Show the row coordinates. By default is TRUE.
}
  \item{tol}{
Tolerance.
}
  \item{DatosIni}{
It contains the initial data. By default is TRUE.
}
}
\details{
The function performs a PERMANOVA Analysis.
}
\value{
The PERMANOVA function create a list that return:
\item{call}{Function.}
\item{Title}{Permutation based MANOVA.}
\item{Type}{PERMANOVA.}
\item{Distances}{A matrix containing the distances between individuals.}
\item{C}{Contrasts Matrix.}
\item{Initial}{ Containing two matrices:

    * Global -> Global contrast.

    * Contrastes ->Contrasts for groups.
}
\item{DistMuestral}{Sample distribution of F-exp from permutations.}
\item{pvalue}{Estimate p-valor for PERMANOVA.}
\item{ExplainedVariance}{Explained variance by Principal Coordinates selected.}
\item{Inertias}{Own value, Explained variance, Cumulative explained variance.}
\item{MeanCoordinates}{Mean Coordinates by groups for the dimensions obtained in the Principal Coordinates Analysis.}
\item{Qualities}{Qualities representation by groups for the dimensions of PCoA.}
\item{CummulativeQualities}{Cumulative qualities representation.}
\item{ClusterType}{Cluster type selected.}
\item{Clusters}{Clusters created.}
\item{ClusterNames}{Names of clusters.}
\item{ClusterColors}{Colors of clusters, color name and HTML code.}
}
\references{
Anderson, M. J. (2008). A new method for non-parametric multivariate analysis of variance. Austral ecology, 26(1):32–46.

Anderson, M. J. (2005). Permanova: a fortran computer program for permutational multivariate analysis of variance. Department of Statistics, University of Auckland, New Zealand, 24.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
\dontrun{data(wine)
X = wine[,4:21]
X=IniTransform(X)
D = DistContinuous (X)
perwine=PERMANOVA(D, wine$Group)
perwine


C = matrix(c(1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1), nrow=3, byrow=TRUE)
rownames(C)=c("C1", "C2", "C3")
colnames(C)=levels(wine$Group)

effects=factor(c(1,2,3))
levels(effects)=c("Origin", "Year", "Interaction")
perwine2=PERMANOVA(D, wine$Group, C=C, Effects=effects, CoordPrinc = TRUE)
summary(perwine2)
}
}
