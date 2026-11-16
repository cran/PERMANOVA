\name{DistContinuous}
\alias{DistContinuous}
\title{
Distances among individuals with continuous data
}
\description{
Calculates distances among individuals for continuous data. It is possible to introduce two matrices (\code{x}, \code{y}) and calculate the distances between the two sets of rows, or introduce only one matrix (\code{x}) and calculate the distances between its rows. Includes classical distances as well as standardized, correlation-based, and robust distances from \pkg{dbrobust}.
}
\usage{
DistContinuous(x, y = NULL, coef = "Pythagorean", r = 1,
  alpha = 0.1, euclidean = FALSE)
}
\arguments{
  \item{x}{
A matrix or data frame containing continuous data.
}
  \item{y}{
An optional matrix containing continuous data different from \code{x} (only supported for classical distances). By default it is \code{NULL}.
}
  \item{coef}{
Coefficient to calculate continuous distances. Options are: \code{"Pythagorean"}, \code{"Taxonomic"}, \code{"City"}, \code{"Minkowski"}, \code{"Divergence"}, \code{"dif_sum"}, \code{"Camberra"}, \code{"Bray_Curtis"}, \code{"Soergel"}, \code{"Ware_Hedges"}, \code{"Gower"}, \code{"Euclidean_Standardized"}, \code{"Maximum"}, \code{"Cosine"}, \code{"Correlation"}, \code{"Mahalanobis"}, \code{"Robust_GV"}, and \code{"Robust_MCD"}. By default \code{"Pythagorean"} is used.
}
  \item{r}{
Power parameter for Minkowski and Taxonomic distances. By default it is 1.
}
  \item{alpha}{
Trimming proportion used when \code{coef = "Robust_GV"} or \code{"Robust_MCD"}. By default it is 0.1.
}
  \item{euclidean}{
Logical; if \code{TRUE}, applies Lingoes Euclidean correction via \code{\link[dbrobust]{make_euclidean}}. By default it is \code{FALSE}.
}
}
\details{
The function calculates distances among individuals for matrices of continuous data. Classical distances follow Gower (2006), whereas standardized, Mahalanobis, and robust trimmed distances (\code{"Robust_GV"} and \code{"Robust_MCD"}) are computed through \pkg{dbrobust}.
}
\value{
A list of class \code{"proximities"} with:
\item{Data}{A matrix with the initial data (\code{x} matrix).}
\item{D}{The matrix of distances.}
\item{Coefficient}{The coefficient used.}
\item{Type}{Character string indicating \code{"continuous"}.}
}
\references{
Boj, E. and Grane, A. (2024). The robustification of distance-based linear models: Some proposals. \emph{Socio-Economic Planning Sciences}, 95, 101992.

Gower, J. C. (2006). Similarity, dissimilarity and Distance, measures of. \emph{Encyclopedia of Statistical Sciences}, 2nd ed., Volume 12. Wiley.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
data(wine)
X <- wine[, 4:21]
D <- DistContinuous(X)
}
