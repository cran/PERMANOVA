\name{DistBinary}
\alias{DistBinary}

\title{
Distances for binary data
}
\description{
Calculates distancies among individuals for binary data. It is possible introduce two matrices (x, y) and calculate the distancies between them or introduce only one matrix (x) and calculate the distancies among individuals in the single matrix.
}
\usage{
DistBinary(x, y = NULL, coefficient = "Simple_Matching", transformation = "sqrt(1-S)")
}

\arguments{
  \item{x}{
A matrix containing binary data.
}
  \item{y}{
A matrix containing binary data different than x. By default it is null because if we want calculate the distancies between individuals it's not necessary.
}
  \item{coefficient}{
Similarity coefficients for binary data.By default it is Simple_Matching. See details.
}
  \item{transformation}{
Transformations of similarities into distances. By default it is sqrt(1-S). See details.
}
}
\details{
The function calculates coefficients of similarity (among individuals) for binary data and converts them into measures of distance.

We have 17 different coefficients (numbers instead of names can be used as arguments of the function) :


   * 1.- Kulezynski

   * 2.- Russell_and_Rao

   * 3.- Jaccard

   * 4.- Simple_Matching

   * 5.- Anderberg

   * 6.- Rogers_and_Tanimoto

   * 7.- Sorensen_Dice_and_Czekanowski

   * 8.- Sneath_and_Sokal

   * 9.- Hamman

   * 10.- Kulezynski2

   * 11.- Anderberg2

   * 12.- Ochiai

   * 13.- S13

   * 14.- Pearson_phi

   * 15.- Yule

   * 16.- Sorensen

   * 17.- Dice

   There are nine possible transformations of the initial similarities:

  * 1.- Identity

  * 2.- 1-S

  * 3.- sqrt(1-S)

  * 4.- (-log(S))

  * 5.- (1/S-1)

  * 6.- (sqrt(2(1-S)))

  * 7.- (1-(S+1)/2)

  * 8.- 1-abs(S)

  * 9.- 1/(S+1)
}
\value{
The DistBinary function creates a list that return:
\item{Data}{A matrix with initial data.}
\item{D}{A matrix containing binary distancies.}
\item{Coefficient}{Character containing the name of the coefficient used.}
}

\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
data(spidersb)
DD=DistBinary(spidersb)
}
