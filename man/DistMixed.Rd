\name{DistMixed}
\alias{DistMixed}
\title{
  Distances for mixed-type data
}
\description{
  Calculates distances among individuals for datasets combining continuous, binary, and categorical variables using Gower's distance or robust distance methods (\code{"ggower"} and \code{"relms"}) from \pkg{dbrobust}.
}
\usage{
DistMixed(x, coef = "gower", cont_vars = NULL, bin_vars = NULL,
  cat_vars = NULL, w = NULL, type = list(), cov_method = "gv",
  alpha = 0.1, euclidean = FALSE)
}
\arguments{
  \item{x}{
A data frame containing mixed-type variables.
}
  \item{coef}{
Distance coefficient: \code{"gower"} (default), \code{"ggower"} (Generalized Gower robust distance), or \code{"relms"} (Related Metric Scaling robust distance).
}
  \item{cont_vars}{
Optional character vector of column names for continuous variables. Detected automatically if \code{NULL}.
}
  \item{bin_vars}{
Optional character vector of column names for binary variables. Detected automatically if \code{NULL}.
}
  \item{cat_vars}{
Optional character vector of column names for categorical variables. Detected automatically if \code{NULL}.
}
  \item{w}{
Optional numeric vector of variable weights (for \code{"gower"}) or observation weights (for \code{"ggower"} and \code{"relms"}).
}
  \item{type}{
Optional list specifying variable types passed to \code{\link[dbrobust]{calculate_distances}} when \code{coef = "gower"}.
}
  \item{cov_method}{
Robust covariance estimation method for continuous variables when \code{coef} is \code{"ggower"} or \code{"relms"}: \code{"gv"} (Geometric Variability, default) or \code{"mcd"} (Minimum Covariance Determinant).
}
  \item{alpha}{
Trimming proportion for robust covariance estimation. By default it is 0.1.
}
  \item{euclidean}{
Logical; if \code{TRUE}, applies Lingoes Euclidean correction via \code{\link[dbrobust]{make_euclidean}}. By default it is \code{FALSE}.
}
}
\value{
A list of class \code{"proximities"} with:
\item{Data}{Data frame with the initial data.}
\item{D}{Matrix of distances.}
\item{Coefficient}{Character string indicating the coefficient used.}
\item{Type}{Character string indicating \code{"mixed"}.}
}
\references{
Boj, E. and Grane, A. (2024). The robustification of distance-based linear models: Some proposals. \emph{Socio-Economic Planning Sciences}, 95, 101992.
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
df_mix <- data.frame(
  v1 = c(1.2, 3.4, 2.1, 5.0),
  v2 = c(0, 1, 1, 0),
  v3 = factor(c("A", "B", "A", "C"))
)
D <- DistMixed(df_mix, coef = "gower")
}
