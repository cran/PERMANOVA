plot.PERMANOVA <- function(x, A1 = 1, A2 = 2, ShowAxes = FALSE, margin = 0.1, ShowBox = TRUE,
                           PlotGroups = TRUE, LabelGroups = TRUE, CexGroup = 1.5, PchGroup = 16,
                           ColorGroup = NULL, PlotInd = FALSE, LabelInd = TRUE, CexInd = 0.8,
                           PchInd = 3, ColorInd = NULL, WhatInds = NULL, IndLabels = NULL,
                           SmartLabels = TRUE, AddLegend = TRUE, LegendPos = "topright",
                           ConvexHulls = FALSE, PlotClus = TRUE, TypeClus = "ch",
                           ClustConf = 1, CexClustCenters = 1, ClustCenters = FALSE,
                           UseClusterColors = TRUE, ...) {

  if (!isTRUE(x$CoordPrinc)) {
    stop("You have to calculate the principal coordinates when running PERMANOVA (CoordPrinc = TRUE).")
  }

  plot.BootDisMANOVA(
    x = x, A1 = A1, A2 = A2, ShowAxes = ShowAxes, margin = margin, ShowBox = ShowBox,
    PlotGroups = PlotGroups, LabelGroups = LabelGroups, CexGroup = CexGroup,
    PchGroup = PchGroup, ColorGroup = ColorGroup, PlotInd = PlotInd, LabelInd = LabelInd,
    CexInd = CexInd, PchInd = PchInd, ColorInd = ColorInd, WhatInds = WhatInds,
    IndLabels = IndLabels, SmartLabels = SmartLabels, AddLegend = AddLegend,
    LegendPos = LegendPos, ConvexHulls = ConvexHulls, PlotClus = PlotClus,
    TypeClus = TypeClus, ClustConf = ClustConf, CexClustCenters = CexClustCenters,
    ClustCenters = ClustCenters, UseClusterColors = UseClusterColors, ...
  )
}
