\name{DistCategorical}
\alias{DistCategorical}
\title{
  Distances for categorical data
}
\description{
  Calculates distances among individuals for categorical (nominal) variables using either the simple matching coefficient or the regularized Cramer's V association-corrected distance from \pkg{dbrobust}.
}
\usage{
DistCategorical(x, coef = "cramer", weights = NULL, euclidean = FALSE)
}
\arguments{
  \item{x}{
A data frame or matrix containing categorical variables. Columns are coerced to factors internally.
}
  \item{coef}{
Distance coefficient: \code{"cramer"} (default, corrects for redundancy among associated categorical variables) or \code{"matching_coefficient"}.
}
  \item{weights}{
Optional numeric vector of observation weights of length \code{nrow(x)} when \code{coef = "cramer"}.
}
  \item{euclidean}{
Logical; if \code{TRUE}, applies Lingoes Euclidean correction via \code{\link[dbrobust]{make_euclidean}}. By default it is \code{FALSE}.
}
}
\value{
A list of class \code{"proximities"} with:
\item{Data}{Data frame with the categorical variables.}
\item{D}{Matrix of distances.}
\item{Coefficient}{Character string indicating the coefficient used.}
\item{Type}{Character string indicating \code{"categorical"}.}
}
\references{
Grane, A., Salini, S., and Infante, G. (2026). New distances for mixed-type data able to cope with redundant information. \emph{AStA Advances in Statistical Analysis}.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
df_cat <- data.frame(
  A = factor(c("low", "med", "high", "low")),
  B = factor(c("yes", "no", "yes", "no"))
)
D <- DistCategorical(df_cat, coef = "cramer")
}
