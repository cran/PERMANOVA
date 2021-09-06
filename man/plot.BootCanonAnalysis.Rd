\name{plot.BootCanonAnalysis}
\alias{plot.BootCanonAnalysis}

\title{
Plots the principal coordinates of the group centers and the bootstrap confidence regions
}
\description{
Plots the principal coordinates of the group centers and the bootstrap confidence regions.
}
\usage{
\method{plot}{BootCanonAnalysis}(x, A1 = 1, A2 = 2, centred = FALSE, confidence = 0.9,
PlotReplicates = TRUE, MeanCex = 1.5, MeanPch = 16, 
Title = "Bootstrap Canonical Analysis based on Distances", 
LabelMeans = TRUE, MeanLabels = NULL, MeanColors = NULL, SmartLabels = TRUE, 
BootstrapPlot = "el", PlotIndiv = FALSE, LabelInd = FALSE, IndLabels = NULL, 
IndColors = NULL, CexInd = 0.5, PchInd = 1, ConvexHullsInd = FALSE, ...)
}

\arguments{
  \item{x}{
The object to plot.
}
  \item{A1}{
Dimension to plot on the X axis.
}
  \item{A2}{
Dimension to plot on the Y axis.
}
  \item{centred}{
Plot centers of the bootstrap regions.
}
  \item{confidence}{
Confidence for the bootstrap regions.
}
  \item{PlotReplicates}{
Should all the bootstrap replicates be plotted.
}
  \item{MeanCex}{
Size of the points and text representing the centroids of the groups.
}
  \item{MeanPch}{
Marks of the points and text representing the centroids of the groups.
}
  \item{Title}{
Main title of the graph.
}
  \item{LabelMeans}{
Should the groups be labelled?
}
  \item{MeanLabels}{
Labels for the group means.
}
  \item{MeanColors}{
Colors for the group means.
}
  \item{SmartLabels}{
Should the labels be placed in a smart way?
}
  \item{BootstrapPlot}{
Should the bootstrap regions be plotted?
}
  \item{PlotIndiv}{
Should the individuals be plotted?
}
  \item{LabelInd}{
Should individuals be labelled?
}
  \item{IndLabels}{
Labels for the individuals.
}
  \item{IndColors}{
Colors for the rows (individuals).
}
  \item{CexInd}{
Sizes for the rows (individuals).
}
  \item{PchInd}{
Marks for the rows (individuals).
}
  \item{ConvexHullsInd}{
Convex hulls containing the individuals of each group.
}
  \item{\dots}{
Any other graphical parameters.
}
}
\details{
Plots the principal coordinates of the means of the groups with the confidence bootstrap region.
}
\value{
The plot.
}

\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}
