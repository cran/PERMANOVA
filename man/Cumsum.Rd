\name{Cumsum}
\alias{Cumsum}

\title{
Cumulative sums
}
\description{
Cumulative sums.
}
\usage{
Cumsum(X, dimens = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data Matrix.
}
  \item{dimens}{
Dimension for summing.
}
}
\details{
Cumulative sums within rows (dimens=1) or columns (dimens=2) of a data matrix.
}
\value{
A matrix of the same size as X with cumulative sums within each row or each column.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
data(wine)
X=wine[,4:21]
Cumsum(X,1)
Cumsum(X,2)
}

