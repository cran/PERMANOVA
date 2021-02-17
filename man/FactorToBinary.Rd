\name{FactorToBinary}
\alias{FactorToBinary}

\title{
Converts a Factor into its indicator matrix
}
\description{
Converts a factor into a binary matrix with as many columns as categories of the factor.
}
\usage{
FactorToBinary(y, Name = NULL)
}

\arguments{
  \item{y}{
A factor.
}
  \item{Name}{
Name to use in the final matrix.
}
}

\value{
An indicator binary matrix.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
y=factor(c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1))
FactorToBinary(y)
}
