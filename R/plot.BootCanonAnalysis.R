plot.BootCanonAnalysis <- function(x, A1 = 1, A2 = 2, centred = FALSE, confidence = 0.90,
                                   PlotReplicates = TRUE, MeanCex = 1.5, MeanPch = 16,
                                   Title = "Bootstrap Canonical Analysis based on Distances",
                                   LabelMeans = TRUE, MeanLabels = NULL, MeanColors = NULL,
                                   SmartLabels = TRUE, BootstrapPlot = "el", PlotIndiv = FALSE,
                                   LabelInd = FALSE, IndLabels = NULL, IndColors = NULL,
                                   CexInd = 0.5, PchInd = 1, ConvexHullsInd = FALSE, ...) {

  L <- dim(x$MeanCoordinates)[1]
  if (is.null(MeanColors)) MeanColors <- grDevices::rainbow(L)
  if (is.null(MeanLabels)) MeanLabels <- x$GroupNames
  if (is.null(IndColors)) IndColors <- MeanColors[as.integer(x$Groups)]
  if (is.null(IndLabels)) IndLabels <- rownames(x$RowCoordinates)

  A <- x$RowCoordinates[, c(A1, A2), drop = FALSE]

  # Calcular límites globales incluyendo medias, réplicas bootstrap e individuos
  all_x <- c(x$MeanCoordinates[, A1], as.vector(x$CoordBoot[, A1, ]))
  all_y <- c(x$MeanCoordinates[, A2], as.vector(x$CoordBoot[, A2, ]))
  if (PlotIndiv || ConvexHullsInd) {
    all_x <- c(all_x, A[, 1])
    all_y <- c(all_y, A[, 2])
  }
  rx <- range(all_x, na.rm = TRUE)
  ry <- range(all_y, na.rm = TRUE)
  mx <- diff(rx) * 0.1
  my <- diff(ry) * 0.1
  xmin <- rx[1] - mx
  xmax <- rx[2] + mx
  ymin <- ry[1] - my
  ymax <- ry[2] + my

  xlabel <- paste("Dim", A1, "(", round(x$Inertia[A1], digits = 2), "%)")
  ylabel <- paste("Dim", A2, "(", round(x$Inertia[A2], digits = 2), "%)")

  extrem <- matrix(c(xmin, ymin, xmax, ymin, xmin, ymax, xmax, ymax), ncol = 2, byrow = TRUE)
  plot(extrem[, 1], extrem[, 2], cex = 0, asp = 1, main = Title, xlab = xlabel, ylab = ylabel)
  points(x$MeanCoordinates[, A1], x$MeanCoordinates[, A2], col = MeanColors, pch = MeanPch, cex = MeanCex)

  CoordinatesMeans <- x$MeanCoordinates[, c(A1, A2), drop = FALSE]
  for (i in 1:L) {
    CoordinatesMeans[i, ] <- rowMeans(matrix(x$CoordBoot[i, c(A1, A2), ], nrow = 2, ncol = x$nB))
  }
  points(CoordinatesMeans[, 1], CoordinatesMeans[, 2], col = MeanColors, pch = 1, cex = 1)

  if (LabelMeans) {
    if (SmartLabels) {
      MultBiplotR::textsmart(
        cbind(x$MeanCoordinates[, A1], x$MeanCoordinates[, A2]),
        MeanLabels, CexPoints = MeanCex, ColorPoints = MeanColors
      )
    } else {
      text(x$MeanCoordinates[, A1], x$MeanCoordinates[, A2],
           labels = MeanLabels, col = MeanColors, cex = MeanCex, pos = 1)
    }
  }

  if (BootstrapPlot == "el") {
    for (i in 1:L) {
      boot_i <- t(matrix(x$CoordBoot[i, c(A1, A2), ], nrow = 2, ncol = x$nB))
      if (centred) {
        ellipse <- MultBiplotR::ConcEllipse(
          boot_i, center = c(x$MeanCoordinates[i, A1], x$MeanCoordinates[i, A2]),
          confidence = confidence
        )
      } else {
        ellipse <- MultBiplotR::ConcEllipse(boot_i, confidence = confidence)
      }
      plot(ellipse, col = MeanColors[i])
    }
  } else {
    for (i in 1:L) {
      boot_i <- t(matrix(x$CoordBoot[i, c(A1, A2), ], nrow = 2, ncol = x$nB))
      if (centred) {
        fraction <- MultBiplotR::Fraction(
          boot_i, center = c(x$MeanCoordinates[i, A1], x$MeanCoordinates[i, A2]),
          confidence = confidence
        )
      } else {
        fraction <- MultBiplotR::Fraction(boot_i, confidence = confidence)
      }
      plot(fraction, type = BootstrapPlot, col = MeanColors[i])
    }
  }

  points(CoordinatesMeans[, 1], CoordinatesMeans[, 2], col = MeanColors, pch = MeanPch, cex = MeanCex)
  difs <- CoordinatesMeans - x$MeanCoordinates[, c(A1, A2), drop = FALSE]

  for (i in 1:L) {
    segments(x$MeanCoordinates[i, A1], x$MeanCoordinates[i, A2],
             CoordinatesMeans[i, 1], CoordinatesMeans[i, 2], col = MeanColors[i])
  }

  if (PlotReplicates) {
    for (i in 1:L) {
      rep_pts <- t(matrix(x$CoordBoot[i, c(A1, A2), ], nrow = 2, ncol = x$nB))
      if (centred) {
        rep_pts <- rep_pts - matrix(1, nrow = x$nB, ncol = 1) %*% difs[i, ]
      }
      points(rep_pts[, 1], rep_pts[, 2], col = MeanColors[i], cex = 0.7, pch = 16)
    }
  }

  if (PlotIndiv) {
    points(A[, 1], A[, 2], cex = CexInd, col = IndColors, pch = PchInd, ...)
  }
  if (LabelInd) {
    if (SmartLabels) {
      MultBiplotR::textsmart(cbind(A[, 1], A[, 2]), Labels = IndLabels,
                             CexPoints = CexInd, ColorPoints = IndColors, ...)
    } else {
      text(A[, 1], A[, 2], labels = IndLabels, cex = CexInd, col = IndColors, pos = 1, ...)
    }
  }

  if (ConvexHullsInd) {
    lev <- levels(x$Groups)
    for (i in seq_along(lev)) {
      XP <- A[which(x$Groups == lev[i]), , drop = FALSE]
      if (nrow(XP) >= 3) {
        hpts <- chull(XP)
        hpts <- c(hpts, hpts[1])
        lines(XP[hpts, , drop = FALSE], col = MeanColors[i])
      }
    }
  }
}




