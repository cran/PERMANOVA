\name{PostHocComp}
\alias{PostHocComp}

\title{
Post Hoc pairwise comparisons
}
\description{
Creates a contrast C matrix of post hoc comparisons among groups.
}
\usage{
PostHocComp(grupo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grupo}{
A factor containing the groups or treatments.
}
}
\details{
Creates a contrast C matrix of post hoc comparisons among groups. The rows of the contrast matrix are not orthogonal.
}
\value{
The PostHocComp function return a matrix to compare the levels of groups introduced.
}

\author{
Laura Vicente-Gonzalez and Jose Luis Vicente-Villardon
}

\examples{
data(wine)
PostHocComp(wine$Group)
}
