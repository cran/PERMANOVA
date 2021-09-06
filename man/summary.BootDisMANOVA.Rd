\name{summary.BootDisMANOVA}
\alias{summary.BootDisMANOVA}

\title{
Summarizes the results of a Bootstrap Manova based on distances
}
\description{
Summarizes the results of a Bootstrap Manova based on distances.
}
\usage{
\method{summary}{BootDisMANOVA}(object, Latex = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "BootDisMANOVA".
}
  \item{Latex}{
Should Latex tables be provided?
}
  \item{\dots}{
Any other parameter.
}
}
\details{
Summarizes the results of a Bootstrap Manova based on distances including Latex tables.
}
\value{
Prints the results.
}

\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
\donttest{data(wine)
X = wine[,4:21]
D = DistContinuous (X)
bootwine=BootDisMANOVA(D, wine$Group)
summary(bootwine)
}
}
