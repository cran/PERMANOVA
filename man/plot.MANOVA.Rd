\encoding{UTF-8}
\name{plot.MANOVA}
\alias{plot.MANOVA}
\title{
Plots the results of a MANOVA Biplot
}
\description{
Plots the results of a MANOVA Biplot.
}
\usage{
\method{plot}{MANOVA}(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE,
PlotVars = TRUE, PlotInd = TRUE, WhatInds = NULL, WhatVars = NULL, WhatGroups = NULL,
IndLabels = NULL, VarLabels = NULL, GroupLabels = NULL, AbbreviateLabels = FALSE,
LabelInd = TRUE, LabelVars = TRUE, LabelGroups = TRUE, CexGroup = 1, PchGroup = 16,
CexInd = 0.8, PchInd = 1, CexVar = 0.8, PchVar = 16, ColorGroups = NULL,
ColorVars = NULL, ColorInd = NULL, margin = 0.1, AddLegend = FALSE,
LegendPos = "topright", ShowAxes = FALSE, LabelAxes = TRUE, ShowBox = FALSE,
ShowTitle = TRUE, PlotCircle = TRUE, TypeCircle = "M", ConvexHulls = FALSE,
voronoi = TRUE, VoronoiColor = "black", mode = "a", TypeScale = "Complete",
ValuesScale = "Original", MinQualityVars = 0, dpg = 0, dp = 0, PredPoints = 0,
SmartLabels = FALSE, PlotClus = FALSE, TypeClus = "ch", ClustConf = 1,
ClustCenters = FALSE, UseClusterColors = TRUE, CexClustCenters = 1, ...)
}
\arguments{
  \item{x}{An object of class \code{"MANOVA"}.}
  \item{A1}{Dimension for the first axis. 1 is the default.}
  \item{A2}{Dimension for the second axis. 2 is the default.}
  \item{ScaleGraph}{Rescale the coordinates to optimal matching.}
  \item{PlotGroups}{Should the group centers be plotted?}
  \item{PlotVars}{Should the variables be plotted?}
  \item{PlotInd}{Should the individuals be plotted?}
  \item{WhatInds}{What individuals should be plotted? A logical or index vector.}
  \item{WhatVars}{What variables should be plotted? A logical or index vector.}
  \item{WhatGroups}{What groups should be plotted? A logical or index vector.}
  \item{IndLabels}{Labels for the individuals.}
  \item{VarLabels}{Labels for the variables.}
  \item{GroupLabels}{Labels for the groups.}
  \item{AbbreviateLabels}{Should the labels be abbreviated for simplicity of the graph?}
  \item{LabelInd}{Should the individuals be labelled?}
  \item{LabelVars}{Should the variables be labelled?}
  \item{LabelGroups}{Should the groups be labelled?}
  \item{CexGroup}{Size of the points for groups.}
  \item{PchGroup}{Marker of the points for groups.}
  \item{CexInd}{Size of the points for individuals.}
  \item{PchInd}{Marker of the points for individuals.}
  \item{CexVar}{Size of the points for variables.}
  \item{PchVar}{Marker of the points for variables.}
  \item{ColorGroups}{User colors for the groups. Default colors will be used if NULL.}
  \item{ColorVars}{User colors for the variables. Default colors will be used if NULL.}
  \item{ColorInd}{User colors for the individuals. Default colors will be used if NULL.}
  \item{margin}{Margin for labels (as a proportion).}
  \item{AddLegend}{Should a legend be added?}
  \item{LegendPos}{Position of the legend.}
  \item{ShowAxes}{Should the axes be shown?}
  \item{LabelAxes}{Should the axes be labelled?}
  \item{ShowBox}{Should a box around the points be plotted?}
  \item{ShowTitle}{Should title be shown?}
  \item{PlotCircle}{Should the confidence regions for the groups be plotted?}
  \item{TypeCircle}{Type of confidence region: Univariate (U), Bonferroni (B), Multivariate (M) or Classical (C). Default is "M".}
  \item{ConvexHulls}{Should the convex hulls containing the individuals for each group be plotted?}
  \item{voronoi}{Should the Voronoi diagram with the prediction regions for each group be plotted?}
  \item{VoronoiColor}{Color for the Voronoi diagram.}
  \item{mode}{Mode of the biplot: "p", "a", "b", "h", "ah" and "s". See details.}
  \item{TypeScale}{Type of scale to use: "Complete", "StdDev" or "BoxPlot". See details.}
  \item{ValuesScale}{Values to show on the scale: "Original" or "Transformed".}
  \item{MinQualityVars}{Minimum quality of representation for a variable to be plotted.}
  \item{dpg}{A set of indices with the variables that will show projections of the group centers onto the variables.}
  \item{dp}{A set of indices with the variables that will show the projections of the individuals.}
  \item{PredPoints}{A vector with integers. The group centers listed in the vector are projected onto all the variables.}
  \item{SmartLabels}{Logical; if \code{TRUE}, non-overlapping smart labels are used.}
  \item{PlotClus}{Should the clusters be plotted?}
  \item{TypeClus}{Type of plot for the clusters (\code{"ch"} - Convex Hull, \code{"el"} - Ellipse or \code{"st"} - Star).}
  \item{ClustConf}{Percent of points included in the cluster. Only the \code{ClustConf} percent of the points nearest to the center will be used to calculate the cluster.}
  \item{ClustCenters}{Should the cluster centers be plotted?}
  \item{UseClusterColors}{Should the cluster colors be used in the plot?}
  \item{CexClustCenters}{Size of the cluster centers.}
  \item{\dots}{Any other graphical parameter.}
}
\details{
The function plots the results of a MANOVA Biplot. The coordinates for Groups, Individuals and Variables can be shown or not on the plot; each of the three can also be labeled separately. There are parameters to control the way each different set of coordinates is plotted and labeled.

There are several modes for plotting the biplot:

"p" .- Points (Rows and Columns are represented by points).

"a" .- Arrows (The traditional representation with points for rows and arrows for columns).

"b" .- The arrows for the columns are extended to both extremes of the plot and labeled outside the plot area.

"h" .- The arrows for the columns are extended to the positive extreme of the plot and labeled outside the plot area.

"ah" .- Same as arrows but labeled outside the plot area.

"s" .- The directions (or biplot axes) have a graded scale for prediction of the original values.

The \emph{TypeScale} argument applies only to the \emph{"s"} mode. There are three types:

"Complete" .- An equally spaced scale covering the whole range of the data is calculated.

"StdDev" .- Mean with one, two and three standard deviations.

"BoxPlot" .- Box-Plot like scale (Median, 25 and 75 percentiles, maximum and minimum values).

The \emph{ValuesScale} argument applies only to the \emph{"s"} mode and controls if the labels show the \emph{Original} or \emph{Transformed} values.

Some of the initial transformations are not compatible with some of the types of biplots and scales. For example, it is not possible to recover by projection the original values when you double centre the data. In that case you have the residuals for interaction and only the transformed values make sense.
}
\value{
A plot is returned.
}
\references{
Gabriel, K. R. (1972). Analysis of meteorological data by means of canonical decomposition and biplots. Journal of Applied Meteorology, 11(7), 1071-1077.

Amaro, I. R., Vicente-Villardon, J. L., & Galindo Villardon, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Sierra, C., Ruiz-Barzola, O., Menendez, M., Demey, J. R., & Vicente-Villardon, J. L. (2017). Geochemical interactions study in surface river sediments at an artisanal mining area by means of Canonical (MANOVA)-Biplot. Journal of Geochemical Exploration, 175, 72-81.
}
\author{
Jose Luis Vicente-Villardon, Laura Vicente-Gonzalez
}
\examples{
\donttest{
data(wine)
X <- wine[, 4:21]
manbip <- MANOVA(X, Group = wine$Group)
plot(manbip, TypeCircle = "B", voronoi = FALSE)
}
}
\keyword{Discriminant Analysis}
\keyword{Canonical Variate Analysis}
\keyword{MANOVA}
