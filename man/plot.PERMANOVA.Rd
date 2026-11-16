\encoding{UTF-8}
\name{plot.PERMANOVA}
\alias{plot.PERMANOVA}
\title{
Plots the results of the PERMANOVA function
}
\description{
Plots the principal coordinates representation of group centroids and individuals from a \code{PERMANOVA} object (computed with \code{CoordPrinc = TRUE}).
}
\usage{
\method{plot}{PERMANOVA}(x, A1 = 1, A2 = 2, ShowAxes = FALSE, margin = 0.1, ShowBox = TRUE,
PlotGroups = TRUE, LabelGroups = TRUE, CexGroup = 1.5, PchGroup = 16,
ColorGroup = NULL, PlotInd = FALSE, LabelInd = TRUE, CexInd = 0.8,
PchInd = 3, ColorInd = NULL, WhatInds = NULL, IndLabels = NULL,
SmartLabels = TRUE, AddLegend = TRUE, LegendPos = "topright",
ConvexHulls = FALSE, PlotClus = TRUE, TypeClus = "ch",
ClustConf = 1, CexClustCenters = 1, ClustCenters = FALSE,
UseClusterColors = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{"PERMANOVA"} obtained with \code{CoordPrinc = TRUE}.}
  \item{A1}{Dimension for the first axis. 1 is the default.}
  \item{A2}{Dimension for the second axis. 2 is the default.}
  \item{ShowAxes}{Should the axes be shown?}
  \item{margin}{Margin around the plot limits (as a proportion).}
  \item{ShowBox}{Should a box around the plot be drawn?}
  \item{PlotGroups}{Should the group centers be plotted?}
  \item{LabelGroups}{Should the groups be labelled?}
  \item{CexGroup}{Size of the points for groups.}
  \item{PchGroup}{Marker of the points for groups.}
  \item{ColorGroup}{User colors for the groups. Default colors will be used if NULL.}
  \item{PlotInd}{Should the individuals be plotted?}
  \item{LabelInd}{Should the individuals be labelled?}
  \item{CexInd}{Size of the points for individuals.}
  \item{PchInd}{Marker of the points for individuals.}
  \item{ColorInd}{User colors for the individuals. Default colors will be used if NULL.}
  \item{WhatInds}{What individuals should be plotted? A logical or index vector.}
  \item{IndLabels}{Labels for the individuals.}
  \item{SmartLabels}{Should the labels be plotted using non-overlapping smart positioning?}
  \item{AddLegend}{Should a legend be added?}
  \item{LegendPos}{Position of the legend.}
  \item{ConvexHulls}{Should the convex hulls containing the individuals for each group be plotted?}
  \item{PlotClus}{Should the clusters be plotted?}
  \item{TypeClus}{Type of plot for the clusters (\code{"ch"} - Convex Hull, \code{"el"} - Ellipse or \code{"st"} - Star).}
  \item{ClustConf}{Confidence level for cluster ellipses or convex hulls.}
  \item{CexClustCenters}{Size of the cluster centers.}
  \item{ClustCenters}{Should the cluster centers be plotted?}
  \item{UseClusterColors}{Should the cluster colors be used for individuals?}
  \item{\dots}{Any other graphical parameter.}
}
\details{
Plots the principal coordinates of the group centers and, optionally, the projected individuals and cluster regions from a \code{PERMANOVA} analysis.
}
\value{
No return value, called for side effects (generates a plot).
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
\examples{
\donttest{
data(wine)
X <- wine[, 4:21]
X <- IniTransform(X)
D <- DistContinuous(X)
perwine <- PERMANOVA(D, wine$Group, CoordPrinc = TRUE)
plot(perwine)
}
}
