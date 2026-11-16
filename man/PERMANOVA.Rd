\name{PERMANOVA}
\alias{PERMANOVA}
\title{
PERMANOVA: MANOVA based on distances
}
\description{
The correct application of MANOVA needs normal and homoscedastic data and the number of variables to be much smaller than the number of individuals, but for many applications these conditions do not hold. To extend the application to such data, Anderson developed PERMANOVA. This non-parametric test based on distances uses permutations to approximate the sampling distribution of the test statistic.

Contrasts and Effects can be added to the calculations, and permutations can be executed in parallel for large datasets.
}
\usage{
PERMANOVA(Distance, group, C = NULL, Effects = NULL, nperm = 1000,
  seed = NULL, CoordPrinc = FALSE, dimens = 2, PCoA = "Standard",
  ProjectInd = TRUE, tol = 1e-04, DatosIni = TRUE,
  PostHoc = "bonferroni", parallel = FALSE, cores = NULL)
}
\arguments{
  \item{Distance}{
An object of class \code{"proximities"} (a list containing the data, the distance matrix \code{D}, and the coefficient used).
}
  \item{group}{
A factor containing the groups to compare.
}
  \item{C}{
Contrast matrix. By default it is \code{NULL} and the identity is used.
}
  \item{Effects}{
A factor with the effects in the rows of \code{C}. By default it is \code{NULL} and each row is considered as an effect.
}
  \item{nperm}{
Number of permutations to perform. By default it is 1000.
}
  \item{seed}{
Seed to start permutations. By default it is \code{NULL}.
}
  \item{CoordPrinc}{
Logical; should the principal coordinates be calculated? By default it is \code{FALSE}.
}
  \item{dimens}{
Number of dimensions to choose in the PCoA. By default it is 2.
}
  \item{PCoA}{
The type of Principal Coordinates Analysis: \code{"Standard"} or \code{"Weighted"}. By default it is \code{"Standard"}.
}
  \item{ProjectInd}{
Logical; show the row coordinates. By default it is \code{TRUE}.
}
  \item{tol}{
Numerical tolerance for eigenvalues. By default it is \code{1e-04}.
}
  \item{DatosIni}{
Logical; whether to use the initial continuous data matrix for individual projection when available. By default it is \code{TRUE}.
}
  \item{PostHoc}{
Multiple testing correction method passed to \code{\link[stats]{p.adjust}} (\code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, or \code{"none"}). By default it is \code{"bonferroni"}.
}
  \item{parallel}{
Logical; if \code{TRUE}, permutations are run in parallel using the \pkg{parallel} package. By default it is \code{FALSE}.
}
  \item{cores}{
Number of CPU cores to use when \code{parallel = TRUE}. If \code{NULL}, a safe number of cores is selected automatically.
}
}
\details{
The function performs a permutational multivariate analysis of variance (PERMANOVA) with optional contrast testing, multi-factor effect decomposition, parallel permutation execution, and Principal Coordinates Analysis (PCoA) biplot representation.
}
\value{
An object of class \code{"PERMANOVA"} containing:
\item{call}{Matched call.}
\item{Title}{Analysis title.}
\item{Type}{Analysis type (\code{"PERMANOVA"}).}
\item{Distances}{Distance matrix between individuals.}
\item{C}{Contrasts matrix.}
\item{Initial}{List containing \code{Global}, \code{Contrastes}, and (optionally) \code{Effects} tables.}
\item{DistMuestral}{Sample distribution of global F-statistics from permutations.}
\item{pvalue}{Estimated p-value for global PERMANOVA.}
\item{ExplainedVariance}{Explained variance by Principal Coordinates selected (when \code{CoordPrinc = TRUE}).}
\item{Inertias}{Eigenvalue, explained variance, and cumulative explained variance.}
\item{MeanCoordinates}{Mean coordinates by groups in PCoA.}
\item{Qualities}{Representation qualities by groups.}
\item{CummulativeQualities}{Cumulative representation qualities.}
}
\references{
Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. \emph{Austral Ecology}, 26(1), 32--46.

Anderson, M. J. (2005). \emph{PERMANOVA: a FORTRAN computer program for permutational multivariate analysis of variance}. Department of Statistics, University of Auckland, New Zealand.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
\donttest{
data(wine)
X <- wine[, 4:21]
X <- IniTransform(X)
D <- DistContinuous(X)
perwine <- PERMANOVA(D, wine$Group)
perwine

C <- matrix(c(1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1), nrow = 3, byrow = TRUE)
rownames(C) <- c("C1", "C2", "C3")
colnames(C) <- levels(wine$Group)

effects <- factor(c(1, 2, 3))
levels(effects) <- c("Origin", "Year", "Interaction")
perwine2 <- PERMANOVA(D, wine$Group, C = C, Effects = effects, CoordPrinc = TRUE)
summary(perwine2)
}
}
