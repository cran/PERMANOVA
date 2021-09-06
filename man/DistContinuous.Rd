\name{DistContinuous}
\alias{DistContinuous}

\title{
Distances among individuals with continuous data
}
\description{
Calculate distances among individuals for continuous data. It is possible introduce two matrices (x, y) and calculate the distances between the two sets of rows or introduce only one matrix (x) and calculate the distances between its rows.
}

\usage{
DistContinuous(x, y = NULL, coef = "Pythagorean", r = 1)
}

\arguments{
  \item{x}{
A matrix containing continuous data.
}
  \item{y}{
A matrix containing continuous data different from x. By default it is null.
}
  \item{coef}{
Coefficient to calculate continuous distances. By default we use Pythagorean distances.
}
  \item{r}{
For Minkowski distances. By default it is 1.
}
}

\details{
The function calculates distances among individuals for matrices of continuous data using different distance measures. If two matrices are provided, distances among individuals, one from the first matrix and another from the second, are calculated. If only one matrix is provided, interdistances among its rows are calculated.

The following coefficients are calculated

1.- Pythagorean =  sqrt(sum((y[i, ] - x[j, ])^2)/p)

2.- Taxonomic =  sqrt(sum(((y[i,]-x[j,])^2)/r^2)/p)

3.- City = sum(abs(y[i,]-x[j,])/r)/p

4.- Minkowski = (sum((abs(y[i,]-x[j,])/r)^t)/p)^(1/t)

5.- Divergence = sqrt(sum((y[i,]-x[j,])^2/(y[i,]+x[j,])^2)/p)

6.- dif_sum = sum(abs(y[i,]-x[j,])/abs(y[i,]+x[j,]))/p

7.- Camberra = sum(abs(y[i,]-x[j,])/(abs(y[i,])+abs(x[j,])))

8.- Bray_Curtis = sum(abs(y[i,]-x[j,]))/sum(y[i,]+x[j,])

9.- Soergel = sum(abs(y[i,]-x[j,]))/sum(apply(rbind(y[i,],x[j,]),2,max))

10.- Ware_hedges = sum(abs(y[i,]-x[j,]))/sum(apply(rbind(y[i,],x[j,]),2,max))
}

\value{A list with:
\item{Data}{A matrix with the initial data (x matrix).}
\item{SupData}{A matrix with the supplementary data (y matrix).}
\item{D}{The matrix of distances.}
\item{Coefficient}{The coefficient used.}
}

\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley
}

\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
data(wine)
X = wine[,4:21]
D=DistContinuous(X)
}
