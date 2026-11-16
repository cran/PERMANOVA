\name{plot.BootDisMANOVA}
\alias{plot.BootDisMANOVA}
\title{Plot method for BootDisMANOVA objects}
\description{
Plots the principal coordinates representation of group centroids and individuals from a \code{BootDisMANOVA} object (when computed with \code{CoordPrinc = TRUE}).
}
\usage{
\method{plot}{BootDisMANOVA}(x, A1 = 1, A2 = 2, ShowAxes = FALSE, margin = 0.1, ShowBox = TRUE,
PlotGroups = TRUE, LabelGroups = TRUE, CexGroup = 1.5, PchGroup = 16,
ColorGroup = NULL, PlotInd = TRUE, LabelInd = TRUE, CexInd = 0.8,
PchInd = 3, ColorInd = NULL, WhatInds = NULL, IndLabels = NULL,
SmartLabels = TRUE, AddLegend = TRUE, LegendPos = "topright",
ConvexHulls = FALSE, PlotClus = TRUE, TypeClus = "ch",
ClustConf = 1, CexClustCenters = 1, ClustCenters = FALSE,
UseClusterColors = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{"BootDisMANOVA"} obtained with \code{CoordPrinc = TRUE}.}
  \item{A1}{Dimension to be plotted on the horizontal axis (default is 1).}
  \item{A2}{Dimension to be plotted on the vertical axis (default is 2).}
  \item{ShowAxes}{Logical; if \code{TRUE}, coordinate axes are drawn.}
  \item{margin}{Proportion of margin added around the plot limits.}
  \item{ShowBox}{Logical; if \code{TRUE}, a box is drawn around the plot.}
  \item{PlotGroups}{Logical; if \code{TRUE}, group centroids are plotted.}
  \item{LabelGroups}{Logical; if \code{TRUE}, group centroids are labeled.}
  \item{CexGroup}{Character expansion size for group centroids.}
  \item{PchGroup}{Plotting symbol for group centroids.}
  \item{ColorGroup}{Vector of colors for group centroids.}
  \item{PlotInd}{Logical; if \code{TRUE}, individual points are plotted.}
  \item{LabelInd}{Logical; if \code{TRUE}, individual points are labeled.}
  \item{CexInd}{Character expansion size for individual points.}
  \item{PchInd}{Plotting symbol for individual points.}
  \item{ColorInd}{Vector of colors for individual points.}
  \item{WhatInds}{Logical or numeric vector indicating which individuals to display.}
  \item{IndLabels}{Optional character vector of labels for individuals.}
  \item{SmartLabels}{Logical; if \code{TRUE}, non-overlapping smart labels are used.}
  \item{AddLegend}{Logical; if \code{TRUE}, a legend for the groups is added.}
  \item{LegendPos}{Position of the legend (e.g., \code{"topright"}).}
  \item{ConvexHulls}{Logical; if \code{TRUE}, convex hulls are drawn around each group.}
  \item{PlotClus}{Logical; if \code{TRUE}, cluster structures are plotted.}
  \item{TypeClus}{Type of cluster representation (\code{"ch"} for convex hull, \code{"el"} for ellipse, \code{"st"} for star).}
  \item{ClustConf}{Confidence level for cluster ellipses or convex hulls.}
  \item{CexClustCenters}{Character expansion size for cluster centers.}
  \item{ClustCenters}{Logical; if \code{TRUE}, cluster centers are drawn.}
  \item{UseClusterColors}{Logical; if \code{TRUE}, cluster colors are applied to individuals.}
  \item{\dots}{Additional graphical parameters passed to \code{plot}, \code{points}, and \code{text}.}
}
\value{
No return value, called for side effects (generates a plot).
}
\seealso{
\code{\link{BootDisMANOVA}}, \code{\link{PERMANOVA}}
}
\examples{
D <- DistContinuous(iris[, 1:4])
res <- BootDisMANOVA(D, iris$Species, nB = 50, seed = 123, CoordPrinc = TRUE)
plot(res, LabelInd = FALSE)
}
