\name{plot.MANOVA.Biplot}
\alias{plot.MANOVA.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the results of a MANOVA Biplot
}
\description{
Plots the results of a MANOVA Biplot
}
\usage{
\method{plot}{MANOVA.Biplot}(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE, 
PlotVars = TRUE, PlotInd = TRUE, WhatInds = NULL, WhatVars = NULL, WhatGroups = NULL, 
IndLabels = NULL, VarLabels = NULL, GroupLabels = NULL, AbbreviateLabels = FALSE, 
LabelInd = TRUE, LabelVars = TRUE, CexGroup = 1, PchGroup = 16, margin = 0.1, 
AddLegend = FALSE, ShowAxes = FALSE, LabelAxes = FALSE, LabelGroups = TRUE, 
PlotCircle = TRUE, ConvexHulls = FALSE, TypeCircle = "M", ColorGroups = NULL, 
ColorVars = NULL, LegendPos = "topright", ColorInd = NULL, voronoi = TRUE, 
mode = "a", TypeScale = "Complete", ValuesScale = "Original", MinQualityVars = 0, 
dpg = 0, dpi = 0, dp = 0, PredPoints = 0, PlotAxis = FALSE, CexInd = NULL, 
CexVar = NULL, PchInd = NULL, PchVar = NULL, ColorVar = NULL, ShowAxis = FALSE, 
VoronoiColor = "black", ShowBox = FALSE, ShowTitle = TRUE, PlotClus = FALSE, 
TypeClus = "ch", ClustConf = 1, ClustCenters = FALSE, UseClusterColors = TRUE, 
CexClustCenters = 1, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "MANOVA.Biplot"
}
  \item{A1}{
Dimension for the first axis. 1 is the default.
}
  \item{A2}{
Dimension for the second axis. 2 is the default.
}
  \item{ScaleGraph}{
Reescale the coordinates to optimal matching.
}
  \item{PlotGroups}{
Shoud the group centers be plotted?
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{WhatInds}{
What individuals should be plooted? A vector.
}
  \item{WhatVars}{
What variables should be plooted? A vector.
}
  \item{WhatGroups}{
What groups should be plooted? A vector.
}
  \item{IndLabels}{
Labels for the individuals
}
  \item{VarLabels}{
Labels for the variables
}
  \item{GroupLabels}{
Labels for the groups
}
  \item{AbbreviateLabels}{
Should the labels be Abbreviated for simplicity of the graph?
}
  \item{LabelInd}{
Should the individuaqls be labelled?
}
  \item{LabelVars}{
Should the variables be labelled?
}
  \item{CexGroup}{
Size of the points for groups
}
  \item{PchGroup}{
Mark of the points for groups
}
  \item{margin}{
Margin for labels (in percentage)
}
  \item{AddLegend}{
Should a legend be added?
}
  \item{ShowAxes}{
Should the axes be shown?
}
  \item{LabelAxes}{
Should the axes be labelled?
}
  \item{LabelGroups}{
Should the groups be labelled?
}
  \item{PlotCircle}{
Should the confidence regions for the groups be plotted?
}
  \item{ConvexHulls}{
Should the convex hulls containing the individuals for each group be plotted?
}
  \item{TypeCircle}{
Type of confidence region: Univariate (U), Bonferroni(B), Multivariate (M) or Classical (C)
}
  \item{ColorGroups}{
User colors for the groups. Default colors will be used if NULL.
}
  \item{ColorVars}{
User colors for the variables. Default colors will be used if NULL.
}
  \item{LegendPos}{
Position of the legend.
}
  \item{ColorInd}{
User colors for the individuals. Default colors will be used if NULL.
}
  \item{voronoi}{
Should the voronoi diagram with the prediction regións for each group be plotted?
}
  \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot"
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed"
}
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be plotted
}
  \item{dpg}{
A set of indices with the variables that will show projections onto the variables 
}
  \item{dpi}{
A set of indices with the individuals that will show the projections 
}
  \item{dp}{
A set of indices with the variables that will show the projections of the individuals
}
\item{PredPoints}{
A vector with integers. The group centers listed in the vector are projected onto all the variables.
}
  \item{PlotAxis}{
Not Used
}
  \item{CexInd}{
Size of the points for individuals.
}
  \item{CexVar}{
Size of the points for variables.
}
  \item{PchInd}{
Marhers of the points for individuals.
}
  \item{PchVar}{
Markers of the points for variables.
}
  \item{ColorVar}{
Colors of the points for variables.
}
  \item{ShowAxis}{
Should axis scales be shown?
}
  \item{VoronoiColor}{
Color for the Voronoi diagram
}
  \item{ShowBox}{
Should a box around the poitns be plotted?
} 
  \item{ShowTitle}{
Should Title be shown?
}
  \item{PlotClus}{
Should the clusters be plotted?
}
  \item{TypeClus}{
Type of cluster to plot (convex hull, ellipse or star)
}
  \item{ClustConf}{
Confidence level for the clusters (percentage plotted)
}
  \item{ClustCenters}{
Should the cluster centres be plotted?
}
  \item{UseClusterColors}{
Use cluster colors for individuals
}
  \item{CexClustCenters}{
Size of the cluster centers
}
  \item{\dots}{
Any other graphical parameter
}
}
\details{
The function plots the results of a MANOVA Biplot. The coordinates for Groups, Individuals and Variables can be shown or not on the plot, each of the three can also be labeled separately.
The are parameters to control the way each different set of coordinates is plotted and labeled.

There are several modes for plotting the biplot.

"p".- Points (Rows and Columns are represented by points)

"a" .- Arrows (The traditional representation with points for rows and arrows for columns)

"b" .- The arrows for the columns are extended to both extremes of the plot and labeled outside the plot area.

"h" .- The arrows for the columns are extended to the positive extreme of the plot and labeled outside the plot area.

"ah" .- Same as arrows but labeled outside the plot area.

"s" .- The directions (or biplot axes) have a graded scale for prediction of the original values.

The \emph{TypeScale} argument applies only to the \emph{"s"} mode. There are three types:

"Complete" .- An equally spaced scale covering the whole range of the data is calculates.

"StdDev" .- Mean with one, two and three stadard deviations 

"BoxPlot" .- Box-Plot like Scale (Median, 25 and 75 percentiles, maximum and minimum values.)


The \emph{ValuesScale} argument applies only to the \emph{"s"} mode and controls if the labels show the  \emph{Original} ot \emph{Transformed} values.

Some of the initial transformations are not compatible with some of the types of biplots and scales. For example, It is not possible to recover by projection the original values when you double centre de data. In that case you have the residuals for interaction and only the transformed values make sense.
}
\value{
A plot is returned
}
\references{
Gabriel, K. R. (1972). Analysis of meteorological data by means of canonical decomposition and biplots. Journal of Applied Meteorology, 11(7), 1071-1077.

Amaro, I. R., Vicente-Villardón, J. L., & Galindo Villardón, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Sierra, C., Ruíz-Barzola, O., Menéndez, M., Demey, J. R., & Vicente-Villardón, J. L. (2017). Geochemical interactions study in surface river sediments at an artisanal mining area by means of Canonical (MANOVA)-Biplot. Journal of Geochemical Exploration, 175, 72-81.
}

\author{
Jose Luis Vicente Villardon, Laura Vicente Gonzalez
}

\examples{
\donttest{data(wine)
X=wine[,4:21]
manbip=MANOVA(X, Group=wine$Group)
plot(manbip, TypeCircle="U", Voronoi=FALSE)
}}
\keyword{ Discriminant Analysis }
\keyword{ Canonical Variate Analysis}
\keyword{ MANOVA}
