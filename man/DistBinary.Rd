\name{DistBinary}
\alias{DistBinary}
\title{
Distances for binary data
}
\description{
Calculates distances among individuals for binary data. It is possible to introduce two matrices (\code{x}, \code{y}) and calculate the distances between them, or introduce only one matrix (\code{x}) and calculate the distances among individuals in the single matrix.
}
\usage{
DistBinary(x, y = NULL, coefficient = "Simple_Matching",
  transformation = "sqrt(1-S)", euclidean = FALSE)
}
\arguments{
  \item{x}{
A matrix containing binary data (0/1).
}
  \item{y}{
A matrix containing binary data different from \code{x}. By default it is \code{NULL}.
}
  \item{coefficient}{
Similarity coefficients for binary data. By default it is \code{"Simple_Matching"}. See details.
}
  \item{transformation}{
Transformations of similarities into distances. By default it is \code{"sqrt(1-S)"}. See details.
}
  \item{euclidean}{
Logical; if \code{TRUE}, applies Lingoes Euclidean correction via \code{\link[dbrobust]{make_euclidean}}. By default it is \code{FALSE}.
}
}
\details{
The function calculates coefficients of similarity among individuals for binary data and converts them into measures of distance.

Available coefficients (numbers instead of names can also be used):
1. \code{Kulezynski}, 2. \code{Russell_and_Rao}, 3. \code{Jaccard}, 4. \code{Simple_Matching}, 5. \code{Anderberg}, 6. \code{Rogers_and_Tanimoto}, 7. \code{Sorensen_Dice_and_Czekanowski}, 8. \code{Sneath_and_Sokal}, 9. \code{Hamman}, 10. \code{Kulezynski2}, 11. \code{Anderberg2}, 12. \code{Ochiai}, 13. \code{S13}, 14. \code{Pearson_phi}, 15. \code{Yule}, 16. \code{Sorensen}, 17. \code{Dice}.

Nine transformations of similarities into distances are available:
1. \code{Identity}, 2. \code{1-S}, 3. \code{sqrt(1-S)}, 4. \code{-log(s)}, 5. \code{1/S-1}, 6. \code{sqrt(2(1-S))}, 7. \code{1-(S+1)/2}, 8. \code{1-abs(S)}, 9. \code{1/(S+1)}.
}
\value{
A list of class \code{"proximities"} with:
\item{Data}{A matrix with initial data.}
\item{D}{A matrix containing binary distances.}
\item{Coefficient}{Character containing the name of the coefficient used.}
\item{Type}{Character string indicating \code{"binary"}.}
}
\references{
Gower, J. C. (2006). Similarity, dissimilarity and Distance, measures of. \emph{Encyclopedia of Statistical Sciences}, 2nd ed., Volume 12. Wiley.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
data(spidersb)
DD <- DistBinary(spidersb)
}
