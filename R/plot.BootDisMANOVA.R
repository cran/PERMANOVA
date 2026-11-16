plot.BootDisMANOVA <- function(x, A1 = 1, A2 = 2, ShowAxes = FALSE, margin = 0.1, ShowBox = TRUE,
                               PlotGroups = TRUE, LabelGroups = TRUE, CexGroup = 1.5, PchGroup = 16,
                               ColorGroup = NULL, PlotInd = TRUE, LabelInd = TRUE, CexInd = 0.8,
                               PchInd = 3, ColorInd = NULL, WhatInds = NULL, IndLabels = NULL,
                               SmartLabels = TRUE, AddLegend = TRUE, LegendPos = "topright",
                               ConvexHulls = FALSE, PlotClus = TRUE, TypeClus = "ch",
                               ClustConf = 1, CexClustCenters = 1, ClustCenters = FALSE,
                               UseClusterColors = TRUE, ...) {

  if (is.null(x$MeanCoordinates)) {
    stop("Principal coordinates were not computed. Re-run with CoordPrinc = TRUE.")
  }

  grp <- if (!is.null(x$Clusters)) x$Clusters else x$Groups
  if (is.null(x$ClusterType) && !is.null(grp)) {
    x <- MultBiplotR::AddCluster2Biplot(x, ClusterType = "us", Groups = grp)
  }

  A <- x$RowCoordinates[, c(A1, A2), drop = FALSE]
  n <- dim(A)[1]
  M <- x$MeanCoordinates[, c(A1, A2), drop = FALSE]
  g <- dim(M)[1]

  if (is.null(rownames(A))) rownames(A) <- paste("i", 1:n, sep = "")
  if (is.null(IndLabels)) IndLabels <- rownames(A)

  # Selección de individuos a representar
  if (is.null(WhatInds)) {
    WhatInds <- rep(TRUE, n)
  } else if (!MultBiplotR::CheckBinaryVector(WhatInds)) {
    AllRows <- rep(FALSE, n)
    AllRows[WhatInds] <- TRUE
    WhatInds <- AllRows
  } else {
    WhatInds <- as.logical(WhatInds)
  }

  clus_cols <- if (!is.null(x$ClusterColors)) x$ClusterColors else grDevices::rainbow(g)

  if (is.null(CexInd)) CexInd <- rep(0.8, n) else if (length(CexInd) == 1) CexInd <- rep(CexInd, n)
  if (is.null(PchInd)) PchInd <- rep(3, n) else if (length(PchInd) == 1) PchInd <- rep(PchInd, n)
  if (is.null(ColorInd)) {
    ColorInd <- if (!is.null(grp)) clus_cols[as.integer(grp)] else rep("black", n)
  } else if (length(ColorInd) == 1) {
    ColorInd <- rep(ColorInd, n)
  }

  if (is.null(CexGroup)) CexGroup <- rep(1.5, g) else if (length(CexGroup) == 1) CexGroup <- rep(CexGroup, g)
  if (is.null(PchGroup)) PchGroup <- rep(16, g) else if (length(PchGroup) == 1) PchGroup <- rep(PchGroup, g)
  if (is.null(ColorGroup)) {
    ColorGroup <- clus_cols
  } else if (length(ColorGroup) == 1) {
    ColorGroup <- rep(ColorGroup, g)
  }

  if (ShowAxes) {
    xaxt <- "s"
    yaxt <- "s"
  } else {
    xaxt <- "n"
    yaxt <- "n"
  }

  # Incluir tanto centroides (M) como individuos visibles (A) en los límites del gráfico
  P <- M
  if (PlotInd && any(WhatInds)) {
    P <- rbind(P, A[WhatInds, , drop = FALSE])
  }

  xmin <- min(P[, 1], na.rm = TRUE)
  xmax <- max(P[, 1], na.rm = TRUE)
  ymin <- min(P[, 2], na.rm = TRUE)
  ymax <- max(P[, 2], na.rm = TRUE)

  P <- rbind(P, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  P <- rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))

  XLabel <- paste("Dimension", A1, " (", round(x$ExplainedVariance[A1], digits = 2), "%)", sep = "")
  YLabel <- paste("Dimension", A2, " (", round(x$ExplainedVariance[A2], digits = 2), "%)", sep = "")

  if (identical(x$Type, "PERMANOVA")) {
    Main <- "PERMANOVA: Graphical representation"
  } else {
    Main <- paste("Canonical Distance Analysis (p-value =", round(x$pvalue, digits = 5), ")")
  }

  plot(P[, 1], P[, 2], cex = 0, asp = 1, xlab = XLabel, ylab = YLabel,
       xaxt = xaxt, yaxt = yaxt, main = Main, axes = ShowAxes, ...)

  if (ShowBox) graphics::box()

  RowColors <- ColorInd
  if (PlotClus && !is.null(x$Clusters)) {
    RowColors <- MultBiplotR::PlotBiplotClusters(
      A, x$Clusters, TypeClus = TypeClus, ClusterColors = clus_cols,
      ClusterNames = x$ClusterNames, centers = ClustCenters, ClustConf = ClustConf,
      CexClustCenters = CexClustCenters, ...
    )
    if (!UseClusterColors) RowColors <- ColorInd
  }

  if (PlotInd && any(WhatInds)) {
    points(A[WhatInds, 1], A[WhatInds, 2],
           cex = CexInd[WhatInds], col = RowColors[WhatInds], pch = PchInd[WhatInds], ...)
  }

  if (LabelInd && any(WhatInds)) {
    if (SmartLabels) {
      MultBiplotR::textsmart(
        A[WhatInds, , drop = FALSE], Labels = IndLabels[WhatInds],
        CexPoints = CexInd[WhatInds], ColorPoints = RowColors[WhatInds], ...
      )
    } else {
      text(A[WhatInds, 1], A[WhatInds, 2], labels = IndLabels[WhatInds],
           cex = CexInd[WhatInds], col = RowColors[WhatInds], pos = 1, ...)
    }
  }

  if (PlotGroups) {
    points(M[, 1], M[, 2], cex = CexGroup, col = ColorGroup, pch = PchGroup, ...)
  }

  if (LabelGroups) {
    if (SmartLabels) {
      MultBiplotR::textsmart(M, Labels = rownames(M),
                             CexPoints = CexGroup, ColorPoints = ColorGroup, ...)
    } else {
      text(M[, 1], M[, 2], labels = rownames(M), cex = CexGroup, col = ColorGroup, pos = 1, ...)
    }
  }

  if (ConvexHulls && !is.null(grp)) {
    lev <- levels(grp)
    for (i in seq_along(lev)) {
      XP <- A[which(grp == lev[i]), , drop = FALSE]
      if (nrow(XP) >= 3) {
        hpts <- chull(XP)
        hpts <- c(hpts, hpts[1])
        lines(XP[hpts, , drop = FALSE], col = ColorGroup[i])
      }
    }
  }

  if (AddLegend && !is.null(grp)) {
    legend(LegendPos, legend = levels(grp), col = ColorGroup, pch = PchGroup, bty = "n")
  }
}
