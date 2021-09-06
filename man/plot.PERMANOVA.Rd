\name{plot.PERMANOVA}
\alias{plot.PERMANOVA}

\title{
Plots the results of the PERMANOVA function
}
\description{
Plots the principal coordinates of the group centers a the bootstrap confidence regions.
}
\usage{
\method{plot}{PERMANOVA}(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, ShowAxis = FALSE,
ShowAxes = FALSE, LabelAxes = TRUE, margin = 0.1, ShowBox = TRUE, PlotGroups = TRUE,
LabelGroups = TRUE, CexGroup = 1.5, PchGroup = 16, ColorGroup = NULL, voronoi = TRUE,
VoronoiColor = "black", PlotInd = FALSE, LabelInd = TRUE, CexInd = 0.8, PchInd = 3,
ColorInd = NULL, WhatInds = NULL, IndLabels = NULL, PlotVars = TRUE, LabelVar = TRUE,
CexVar = NULL, PchVar = NULL, ColorVar = NULL, WhatVars = NULL, VarLabels = NULL,
mode = "a", TypeScale = "Complete", ValuesScale = "Original", SmartLabels = TRUE,
AddLegend = TRUE, LegendPos = "topright", PlotCircle = TRUE, ConvexHulls = FALSE,
TypeCircle = "M", MinQualityVars = 0, dpg = 0, dpi = 0, PredPoints = 0,
PlotClus = TRUE, TypeClus = "ch", ClustConf = 1, CexClustCenters=1,ClustCenters = FALSE, 
UseClusterColors = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class "PERMANOVA"
}
  \item{A1}{
Dimension for the first axis. 1 is the default.
}
  \item{A2}{
Dimension for the second axis. 2 is the default.
}
  \item{ScaleGraph}{
Resale the coordinates to optimal matching.
}
  \item{ShowAxis}{
Should the axis passing trough the origin be plotted?
}
  \item{ShowAxes}{
Should the axes be shown?
}
  \item{LabelAxes}{
Should the axes be labelled?
}
  \item{margin}{
Margin for labels (in percentage).
}
  \item{ShowBox}{
Should a box around the points be plotted?
}
  \item{PlotGroups}{
Should the group centers be plotted?
}
  \item{LabelGroups}{
Should the groups be labelled?
}
  \item{CexGroup}{
Size of the points for groups.
}
  \item{PchGroup}{
Mark of the points for groups.
}
  \item{ColorGroup}{
User colors for the groups. Default colors will be used if NULL.
}
  \item{voronoi}{
Should the Voronoi diagram with the prediction regions for each group be plotted?
}
  \item{VoronoiColor}{
Color for the Voronoi diagram.
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelInd}{
Should the individuals be labelled?
}
  \item{CexInd}{
Size of the points for individuals.
}
  \item{PchInd}{
Markers of the points for individuals.
}
  \item{ColorInd}{
User colors for the individuals. Default colors will be used if NULL.
}
  \item{WhatInds}{
What individuals should be plotted? A vector.
}
  \item{IndLabels}{
Labels for the individuals.
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{LabelVar}{
Should the variables be labelled?
}
  \item{CexVar}{
Size of the points for variables.
}
  \item{PchVar}{
Markers of the points for variables.
}
  \item{ColorVar}{
Colors of the points for variables.
}
  \item{WhatVars}{
What variables should be plotted? A vector.
}
  \item{VarLabels}{
Labels for the variables.
}
  \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot".
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed".
}
  \item{SmartLabels}{
Should the labels be plotted in a smart way?
}
  \item{AddLegend}{
Should a legend be added?
}
  \item{LegendPos}{
Position of the legend.
}
  \item{PlotCircle}{
Should the confidence regions for the groups be plotted?
}
  \item{ConvexHulls}{
Should the convex hulls containing the individuals for each group be plotted?
}
  \item{TypeCircle}{
Type of confidence region: Univariate (U), Bonferroni(B), Multivariate (M) or Classical (C). By default is "M".
}
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be plotted.
}
  \item{dpg}{
A set of indices with the variables that will show projections onto the variables.
}
  \item{dpi}{
A set of indices with the individuals that will show the projections.
}
  \item{PredPoints}{
A vector with integers. The group centers listed in the vector are projected onto all the variables.
}
  \item{PlotClus}{
Should the clusters be plotted?
}
  \item{TypeClus}{
	
Type of plot for the clusters. ("ch"- Convex Hull, "el"- Ellipse or "st"- Star)
}
  \item{ClustConf}{
Percent of points included in the cluster. only the ClusConf percent of the points nearest to the center will be used to calculate the cluster.
}
  \item{CexClustCenters}{
Size of the cluster centers.
}
  \item{ClustCenters}{
Should the cluster centers be plotted.
}
  \item{UseClusterColors}{
Should the cluster colors be used in the plot. 
}
  \item{\dots}{
Any other graphical parameter.
}
}
\details{
Plots the principal coordinates of the group centers a the bootstrap confidence regions.
}
\value{
The plot.
}

\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}


\examples{
\donttest{
data(wine)
X = wine[,4:21]
X=IniTransform(X)
D = DistContinuous (X)
perwine=PERMANOVA(D, wine$Group, CoordPrinc=TRUE)
plot(perwine)

}
}
