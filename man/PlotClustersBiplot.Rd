\name{PlotClustersBiplot}
\alias{PlotClustersBiplot}

\title{
Plot clusters on a biplot.
}
\description{
Highlights several groups or clusters on a biplot representation.
}
\usage{
PlotClustersBiplot(A, Groups = Ones(c(nrow(A), 1)), TypeClus = "st", ClusterColors =
NULL, ClusterNames = NULL, centers = TRUE, ClustConf = 1, Legend = FALSE, LegendPos =
"topright", CexClustCenters = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Coordinates of the points in the scattergram.
}
  \item{Groups}{
Factor defining the groups to be highlighted.
}
  \item{TypeClus}{
Type of representation of the clusters. For the moment just a convex hull but in the future ellipses and stars will be added.
}
  \item{ClusterColors}{
A vector of colors with as many elements as clusters. If \code{NULL} the function selects the rainbow colors.
}
  \item{ClusterNames}{
A vector of names with as many elements as clusters.
}
  \item{centers}{
Logical variable to control if centres of the clusters are plotted.
}
  \item{ClustConf}{
Percent of points included in the cluster. Only the ClusConf percent of the points nearest to the center will be used to calculate the cluster.
}
  \item{Legend}{
Should a legend be plotted?
}
  \item{LegendPos}{
Position of the legend.
}
  \item{CexClustCenters}{
Size of the centers.
}
  \item{\dots}{
Any other graphical parameters.
}
}
\details{
The clusters to plot should be added to the biplot object using the function \code{AddClusterToBiplot}.
}
\value{
It takes effects on a plot.
}

\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{AddClusterToBiplot}}
}
