\name{BinaryVectorCheck}
\alias{BinaryVectorCheck}

\title{
Checks if a vector is binary
}
\description{
Checks if all the entries of a vector are 0 or 1.
}
\usage{
BinaryVectorCheck(x)
}

\arguments{
  \item{x}{
The vector to check
}
}

\value{
The logical result.
}

\author{
Jose luis Vicente-Villardon
}

\examples{
x=c(0, 0, 0, 0,  1, 1, 1, 2)
BinaryVectorCheck(x)
}
