\name{ConstructContrasts}
\alias{ConstructContrasts}

\title{
Construction of contrasts for several factors
}
\description{
Constructs the contrasts for main effects and interaction for a dataframe of factors.
}
\usage{
ConstructContrasts(Factors, MaxOrderIter = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Factors}{
A data frame of factors (For general data frames the functions extracts all the columns that are factors).
}
  \item{MaxOrderIter}{
Maximum order of the interaction to construct.
}
}
\details{
Constructs the contrasts for main effects and interaction for a dataframe of factors. For a general data frame the functions extracts all the factors first. The function constructs the contrasts for all the main effects and the interactions up to order MaxOrderIter.
}
\value{
\item{Groups}{The groups resulting of combining the levels of all the factors.}
\item{Contrasts}{A matrix of contrasts containing the contrasts that are associated to the main affects and interaction. The rows contain contrasts and the columns the groups to combine.}
\item{Effects}{A factor describing which contrasts correspond to each effect.}
}
\references{
Bapat, R. B. (2012). Linear algebra and linear models. Springer Science & Business Media.
}
\author{
Jose Luis Vicente-Villardon, Laura Vicente-Gonzalez
}

\examples{
ConstructContrasts(wine[,1:2])
}
