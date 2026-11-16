plot.MANOVA <- function(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE,
                        PlotVars = TRUE, PlotInd = TRUE, WhatInds = NULL,
                        WhatVars = NULL, WhatGroups = NULL, IndLabels = NULL,
                        VarLabels = NULL, GroupLabels = NULL, AbbreviateLabels = FALSE,
                        LabelInd = TRUE, LabelVars = TRUE, LabelGroups = TRUE,
                        CexGroup = 1, PchGroup = 16, CexInd = 0.8, PchInd = 1,
                        CexVar = 0.8, PchVar = 16, ColorGroups = NULL, ColorVars = NULL,
                        ColorInd = NULL, margin = 0.1, AddLegend = FALSE,
                        LegendPos = "topright", ShowAxes = FALSE, LabelAxes = TRUE,
                        ShowBox = FALSE, ShowTitle = TRUE, PlotCircle = TRUE,
                        TypeCircle = "M", ConvexHulls = FALSE, voronoi = TRUE,
                        VoronoiColor = "black", mode = "a", TypeScale = "Complete",
                        ValuesScale = "Original", MinQualityVars = 0, dpg = 0,
                        dp = 0, PredPoints = 0, SmartLabels = FALSE,
                        PlotClus = FALSE, TypeClus = "ch", ClustConf = 1,
                        ClustCenters = FALSE, UseClusterColors = TRUE,
                        CexClustCenters = 1, ...) {

  plot.MANOVA.Biplot(
    x = x, A1 = A1, A2 = A2, ScaleGraph = ScaleGraph, PlotGroups = PlotGroups,
    PlotVars = PlotVars, PlotInd = PlotInd, WhatInds = WhatInds,
    WhatVars = WhatVars, WhatGroups = WhatGroups, IndLabels = IndLabels,
    VarLabels = VarLabels, GroupLabels = GroupLabels, AbbreviateLabels = AbbreviateLabels,
    LabelInd = LabelInd, LabelVars = LabelVars, LabelGroups = LabelGroups,
    CexGroup = CexGroup, PchGroup = PchGroup, CexInd = CexInd, PchInd = PchInd,
    CexVar = CexVar, PchVar = PchVar, ColorGroups = ColorGroups, ColorVars = ColorVars,
    ColorInd = ColorInd, margin = margin, AddLegend = AddLegend,
    LegendPos = LegendPos, ShowAxes = ShowAxes, LabelAxes = LabelAxes,
    ShowBox = ShowBox, ShowTitle = ShowTitle, PlotCircle = PlotCircle,
    TypeCircle = TypeCircle, ConvexHulls = ConvexHulls, voronoi = voronoi,
    VoronoiColor = VoronoiColor, mode = mode, TypeScale = TypeScale,
    ValuesScale = ValuesScale, MinQualityVars = MinQualityVars, dpg = dpg,
    dp = dp, PredPoints = PredPoints, SmartLabels = SmartLabels,
    PlotClus = PlotClus, TypeClus = TypeClus, ClustConf = ClustConf,
    ClustCenters = ClustCenters, UseClusterColors = UseClusterColors,
    CexClustCenters = CexClustCenters, ...
  )
}

