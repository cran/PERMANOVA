plot.MANOVA.Biplot <- function(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE,
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

  if (ncol(x$GroupCoordinates) < 2) {
    stop("At least 2 canonical dimensions (3 groups and 2 variables) are required for a 2D biplot.")
  }

  modes <- c("p", "a", "b", "h", "ah", "s")
  if (is.numeric(mode)) mode <- modes[mode]
  TypeScales <- c("Complete", "StdDev", "BoxPlot")
  if (is.numeric(TypeScale)) TypeScale <- TypeScales[TypeScale]
  ValuesScales <- c("Original", "Transformed")
  if (is.numeric(ValuesScale)) ValuesScale <- ValuesScales[ValuesScale]

  grp <- if (!is.null(x$Clusters)) x$Clusters else x$groups
  clus_cols <- if (!is.null(x$ClusterColors)) x$ClusterColors else grDevices::rainbow(x$g)

  if (is.null(ColorGroups)) ColorGroups <- clus_cols
  if (is.null(ColorInd)) {
    ColorInd <- if (!is.null(grp)) clus_cols[as.integer(grp)] else rep("black", nrow(x$RowCoordinates))
  }

  if (LabelAxes) {
    xlabel <- paste("CV", A1, "(", round(x$Inertia[A1], digits = 2), "%)")
    ylabel <- paste("CV", A2, "(", round(x$Inertia[A2], digits = 2), "%)")
  } else {
    xlabel <- ""
    ylabel <- ""
  }

  if (ShowTitle) {
    Title <- paste(x$Title, "/", A1, "-", A2, "(", round(x$Inertia[A1] + x$Inertia[A2], digits = 2), "%)")
  } else {
    Title <- ""
  }

  scf <- 1
  if (ScaleGraph) {
    sca <- sum(x$GroupCoordinates^2) / x$g
    scb <- sum(x$ColCoordinates^2) / x$p
    scf <- sqrt(sqrt(scb / sca))
  }

  J <- x$GroupCoordinates[, c(A1, A2), drop = FALSE] * scf
  H <- x$ColCoordinates[, c(A1, A2), drop = FALSE] / scf
  V <- x$RowCoordinates[, c(A1, A2), drop = FALSE] * scf

  g <- dim(J)[1]
  p <- dim(H)[1]
  n <- dim(V)[1]

  if (is.null(GroupLabels)) GroupLabels <- rownames(J)
  if (is.null(IndLabels)) IndLabels <- rownames(V)
  if (is.null(VarLabels)) VarLabels <- rownames(H)

  if (AbbreviateLabels) {
    GroupLabels <- abbreviate(GroupLabels, minlength = 5L)
    IndLabels <- abbreviate(IndLabels, minlength = 5L)
    VarLabels <- abbreviate(VarLabels, minlength = 5L)
  }

  # Filtros de individuos, variables y grupos
  if (is.null(WhatInds)) {
    WhatInds <- rep(TRUE, n)
  } else if (!MultBiplotR::CheckBinaryVector(WhatInds)) {
    AllRows <- rep(FALSE, n)
    AllRows[WhatInds] <- TRUE
    WhatInds <- AllRows
  } else {
    WhatInds <- as.logical(WhatInds)
  }

  if (is.null(WhatVars)) {
    WhatVars <- rep(TRUE, p)
  } else if (!MultBiplotR::CheckBinaryVector(WhatVars)) {
    AllCols <- rep(FALSE, p)
    AllCols[WhatVars] <- TRUE
    WhatVars <- AllCols
  } else {
    WhatVars <- as.logical(WhatVars)
  }

  if (is.null(WhatGroups)) {
    WhatGroups <- rep(TRUE, g)
  } else if (!MultBiplotR::CheckBinaryVector(WhatGroups)) {
    AllGrps <- rep(FALSE, g)
    AllGrps[WhatGroups] <- TRUE
    WhatGroups <- AllGrps
  } else {
    WhatGroups <- as.logical(WhatGroups)
  }

  if (length(CexInd) == 1) CexInd <- rep(CexInd, n)
  if (length(PchInd) == 1) PchInd <- rep(PchInd, n)
  if (length(ColorInd) == 1) ColorInd <- rep(ColorInd, n)

  if (length(CexGroup) == 1) CexGroup <- rep(CexGroup, g)
  if (length(PchGroup) == 1) PchGroup <- rep(PchGroup, g)
  if (length(ColorGroups) == 1) ColorGroups <- rep(ColorGroups, g)

  if (is.null(CexVar)) CexVar <- rep(0.8, p) else if (length(CexVar) == 1) CexVar <- rep(CexVar, p)
  if (is.null(PchVar)) PchVar <- rep(16, p) else if (length(PchVar) == 1) PchVar <- rep(PchVar, p)
  if (is.null(ColorVars)) ColorVars <- rep("black", p) else if (length(ColorVars) == 1) ColorVars <- rep(ColorVars, p)

  if (!LabelVars) VarLabels <- rep(" ", p)
  if (!is.null(x$Sup_Individual_Coord)) VS <- x$Sup_Individual_Coord * scf

  if (ShowAxes) {
    xaxt <- "s"
    yaxt <- "s"
  } else {
    xaxt <- "n"
    yaxt <- "n"
  }

  PP <- rbind(J[WhatGroups, , drop = FALSE], H[WhatVars, , drop = FALSE])
  if ((PlotInd || PlotClus) && any(WhatInds)) {
    PP <- rbind(PP, V[WhatInds, , drop = FALSE])
  }

  xmin <- min(PP[, 1], na.rm = TRUE)
  xmax <- max(PP[, 1], na.rm = TRUE)
  ymin <- min(PP[, 2], na.rm = TRUE)
  ymax <- max(PP[, 2], na.rm = TRUE)
  xrang <- abs(xmax - xmin)
  yrang <- abs(ymax - ymin)
  fact <- 0.03
  xmin <- xmin - fact * xrang
  xmax <- xmax + fact * xrang
  ymin <- ymin - fact * yrang
  ymax <- ymax + fact * yrang

  PP <- rbind(PP, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  PP <- rbind(PP, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))

  plot(PP[, 1], PP[, 2], cex = 0, asp = 1, main = Title, xlab = xlabel, ylab = ylabel,
       xaxt = xaxt, yaxt = yaxt, axes = ShowAxes, ...)

  if (ShowBox) rect(xmin, ymin, xmax, ymax)

  if (PlotClus && !is.null(grp)) {
    ClusCols <- MultBiplotR::PlotBiplotClusters(
      V, grp, TypeClus = TypeClus, ClusterColors = clus_cols,
      ClusterNames = x$ClusterNames, centers = ClustCenters, ClustConf = ClustConf,
      CexClustCenters = CexClustCenters, ...
    )
    if (!is.null(x$ClusterType) && x$ClusterType == "gm") {
      ColorInd2 <- rgb((x$P %*% t(col2rgb(clus_cols))) / 255)
      if (UseClusterColors) ColorInd <- ColorInd2
      PchInd <- rep(16, n)
    } else if (UseClusterColors) {
      ColorInd <- ClusCols
    }
  }

  if (PlotInd && any(WhatInds)) {
    points(V[WhatInds, 1], V[WhatInds, 2], col = ColorInd[WhatInds],
           pch = PchInd[WhatInds], cex = CexInd[WhatInds])
    if (LabelInd) {
      if (SmartLabels) {
        MultBiplotR::textsmart(V[WhatInds, , drop = FALSE], Labels = IndLabels[WhatInds],
                               CexPoints = CexInd[WhatInds], ColorPoints = ColorInd[WhatInds])
      } else {
        text(V[WhatInds, 1], V[WhatInds, 2], IndLabels[WhatInds],
             adj = 0, col = ColorInd[WhatInds], cex = CexInd[WhatInds])
      }
    }
    if (!is.null(x$Sup_Individual_Coord)) {
      points(VS[, A1], VS[, A2], cex = 0.5, col = "black")
      text(VS[, A1], VS[, A2], paste(rownames(VS), "*", sep = ""), adj = 0, col = "black", cex = 0.8)
    }
  }

  if (PlotGroups && any(WhatGroups)) {
    points(J[WhatGroups, 1], J[WhatGroups, 2], cex = CexGroup[WhatGroups],
           pch = PchGroup[WhatGroups], col = ColorGroups[WhatGroups])
    if (LabelGroups) {
      if (SmartLabels) {
        MultBiplotR::textsmart(J[WhatGroups, , drop = FALSE], Labels = GroupLabels[WhatGroups],
                               CexPoints = CexGroup[WhatGroups], ColorPoints = ColorGroups[WhatGroups])
      } else {
        text(J[WhatGroups, 1], J[WhatGroups, 2], GroupLabels[WhatGroups],
             adj = 0, col = ColorGroups[WhatGroups], cex = CexGroup[WhatGroups] * 1.2)
      }
    }
  }

  if (AddLegend && any(WhatGroups)) {
    legend(LegendPos, legend = GroupLabels[WhatGroups], col = ColorGroups[WhatGroups],
           pch = PchGroup[WhatGroups], bty = "n")
  }

  if (PlotVars) {
    qlrcols <- x$ColContributions[, A1] + x$ColContributions[, A2]
    Scales <- NULL
    if (mode == "s") {
      Scales <- MultBiplotR::GetBiplotScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)
    }

    for (j in seq_len(p)) {
      if (WhatVars[j] && qlrcols[j] > MinQualityVars) {
        MultBiplotR::VarBiplot(
          H[j, 1], H[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          label = VarLabels[j], mode = mode, CexPoint = CexVar[j], Color = ColorVars[j],
          ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale,
          PchPoint = PchVar[j], ...
        )
      }
    }

    for (idp in dpg) {
      if ((idp > 0) && (idp <= p)) {
        gv <- H[idp, ]
        nn <- sum(gv^2)
        scal <- as.vector(J %*% gv) / nn
        Fpr <- outer(scal, as.vector(gv))
        MultBiplotR::dlines(J, Fpr, color = ColorGroups)
      }
    }

    if (PlotInd) {
      for (idp in dp) {
        if ((idp > 0) && (idp <= p)) {
          gv <- H[idp, ]
          nn <- sum(gv^2)
          scal <- as.vector(V %*% gv) / nn
          Fpr <- outer(scal, as.vector(gv))
          MultBiplotR::dlines(V, Fpr, color = ColorInd)
        }
      }
    }

    for (idp in PredPoints) {
      if ((idp > 0) && (idp <= g)) {
        for (j in seq_len(p)) {
          if (WhatVars[j]) {
            gv <- H[j, ]
            nn <- sum(gv^2)
            scal <- sum(J[idp, ] * gv) / nn
            Fpr <- matrix(scal * gv, nrow = 1, ncol = 2)
            MultBiplotR::dlines(matrix(J[idp, ], 1, 2), Fpr, color = ColorVars[j])
          }
        }
      }
    }
  }

  if (PlotCircle) {
    radius <- switch(TypeCircle,
                     U = x$UnivRad,
                     B = x$BonfRad,
                     M = x$MultRad,
                     C = x$ChisRad
    )
    radius <- radius * scf
    for (i in seq_len(x$g)) {
      if (WhatGroups[i]) {
        MultBiplotR::Circle(radius[i], c(J[i, 1], J[i, 2]), col = ColorGroups[i])
      }
    }
  }

  if (ConvexHulls && !is.null(grp)) {
    lev <- levels(grp)
    for (i in seq_along(lev)) {
      XP <- V[which(grp == lev[i]), , drop = FALSE]
      if (nrow(XP) >= 3) {
        hpts <- chull(XP)
        hpts <- c(hpts, hpts[1])
        lines(XP[hpts, , drop = FALSE], col = ColorGroups[i])
      }
    }
  }

  if (voronoi && g >= 2) {
    tv <- deldir::deldir(J[, 1], J[, 2], rw = c(xmin, xmax, ymin, ymax))
    plot(tv, add = TRUE, wlines = "tess", xlim = c(xmin, xmax), ylim = c(ymin, ymax), cmpnt_col = VoronoiColor)
  }

  if (!is.null(x$ContSupVarsBiplot)) {
    plot(x$ContSupVarsBiplot, F1 = A1, F2 = A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
         mode = mode, TypeScale = TypeScale)
  }
}

