plot.BootCanonAnalysis <- function(x, A1=1, A2=2, centred=FALSE, confidence=0.90, PlotReplicates=TRUE, MeanCex=1.5, MeanPch=16,
                                   Title="Bootstrap Canonical Analysis based on Distances",
                                   LabelMeans=TRUE, MeanLabels = NULL, MeanColors = NULL, SmartLabels = TRUE,
                                   BootstrapPlot="el", PlotIndiv=FALSE, LabelInd = FALSE, IndLabels=NULL, IndColors = NULL, CexInd=0.5,
                                   PchInd=1, ConvexHullsInd=FALSE, ...){

  # BootstrapPlot=c("el", "ch", "st")
  # ColColors=NULL, ColLabels=NULL, SizeQualInd = FALSE
  # ColorQualInd = FALSE, ColorQual="black", PlotSup=TRUE, Bootstrap=FALSE,
  #  margin=0,
  # PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, CexClustCenters=1, LegendClust=TRUE,
  # ClustCenters = FALSE, UseClusterColors = TRUE, ShowAxis=FALSE, PlotBinaryMeans=FALSE,
  # MinIncidence=0, ShowBox=FALSE,  ColorSupContVars=NULL, ColorSupBinVars=NULL, ColorSupOrdVars=NULL,
  # TypeScale = "Complete", SupMode="s", PlotSupVars=FALSE, ...

  L=dim(x$MeanCoordinates)[1]
  if (is.null(MeanColors)) MeanColors=1:L
  if (is.null(MeanLabels)) MeanLabels=x$GroupNames
  if (is.null(IndColors)) IndColors= MeanColors[as.integer(x$Groups)]
  if (is.null(IndLabels)) IndLabels = rownames(x$RowCoordinates)
  A = x$RowCoordinates[, c(A1, A2)]
  xmax=0
  xmin=0
  ymax=0
  ymin=0
  for (i in 1:x$nB){
    xmax=max(xmax,max(x$CoordBoot[,A1,i]))
    xmin=min(xmin,min(x$CoordBoot[,A1,i]))
    ymax=max(ymax,max(x$CoordBoot[,A2,i]))
    ymin=min(ymin,min(x$CoordBoot[,A2,i]))
  }

  xlabel = paste("Dim", A1, "(", round(x$Inertia[A1], digits=2),"%)")
  ylabel = paste("Dim", A2, "(", round(x$Inertia[A2], digits=2),"%)")

  extrem=matrix(c(xmin, ymin, xmax, ymin, xmin, ymax, xmax, ymax), ncol=2, byrow=TRUE)
  plot(extrem[,1], extrem[,2], cex=0, asp=1, main=Title, xlab = xlabel, ylab = ylabel)
  points(x$MeanCoordinates[,A1], x$MeanCoordinates[,A2], col=MeanColors, pch=MeanPch, cex=MeanCex)

  CoordinatesMeans=x$MeanCoordinates[,c(A1,A2)]

  for (i in 1:L){
    CoordinatesMeans[i,]=apply(t(x$CoordBoot[i,c(A1,A2),]), 2, mean)
  }
  points(CoordinatesMeans[,A1], CoordinatesMeans[,A2], col=MeanColors, pch=1, cex=1)

  if (LabelMeans)
    if (SmartLabels)
      TextSmart(cbind(x$MeanCoordinates[,A1], x$MeanCoordinates[,A2]), MeanLabels, CexPoints = MeanCex, ColorPoints = MeanColors)
  else text(x$MeanCoordinates[,A1], x$MeanCoordinates[,A2], labels = MeanLabels, col=MeanColors, cex=MeanCex, pos=1)

  if (BootstrapPlot=="el"){
    for (i in 1:L){
      if (centred)
        ellipse=ConcEllipse( t(x$CoordBoot[i,c(A1,A2),]), center= c(x$MeanCoordinates[i,A1], x$MeanCoordinates[i,A2]), confidence=confidence)
      else
        ellipse=ConcEllipse( t(x$CoordBoot[i,c(A1,A2),]), confidence=confidence)
      plot(ellipse , col=MeanColors[i])
      CoordinatesMeans[i,]=apply(t(x$CoordBoot[i,c(A1,A2),]), 2, mean)
      }
  }
  else{
    for (i in 1:L){
      if (centred)
        fraction=Fractions(t(x$CoordBoot[i,c(A1,A2),]), center= c(x$MeanCoordinates[i,A1], x$MeanCoordinates[i,A2]), confidence=confidence)
      else
        fraction=Fractions(t(x$CoordBoot[i,c(A1,A2),]), confidence=confidence)
      plot(fraction, type=BootstrapPlot , col=MeanColors[i])
      CoordinatesMeans[i,]=apply(t(x$CoordBoot[i,c(A1,A2),]), 2, mean)}
  }

  points(CoordinatesMeans[,1], CoordinatesMeans[,2], col=MeanColors, pch=MeanPch, cex=MeanCex)
  difs=CoordinatesMeans-x$MeanCoordinates

  for (i in 1:L)
    arrows(x$MeanCoordinates[i,A1], x$MeanCoordinates[i,A2],CoordinatesMeans[i,1], CoordinatesMeans[i,2],  col=MeanColors[i])

  if (PlotReplicates)
    for (i in 1:L){
      if (centred)
        rep=t(x$CoordBoot[i,c(A1,A2),])- matrix(1, nrow=x$nB, ncol=1)%*%difs[i,]
      else
        rep=t(x$CoordBoot[i,c(A1,A2),])
      points(rep[,A1], rep[,A2], col=MeanColors[i], cex=0.7, pch=16)
    }


  if (PlotIndiv)
    points(x$RowCoordinates[, A1], x$RowCoordinates[, A2], cex = CexInd, col = IndColors, pch = PchInd, ...)
  if (LabelInd)
    if (SmartLabels)
      TextSmart(cbind(A[, A1], A[, A2]), CexPoints = CexInd, ColorPoints = IndColors, ...)
  else text(A[, A1], A[, A2], rownames(A), cex = CexInd, col = IndColors, pos = 1, ...)

  if (ConvexHullsInd) {
    lev = levels(x$Groups)
    for (i in 1:nlevels(x$Groups)) {
      XP = x$RowCoordinates[which(x$Groups == lev[i]), ]
      XP = cbind(XP[, 1], XP[, 2])
      hpts <- chull(XP)
      hpts <- c(hpts, hpts[1])
      lines(XP[hpts, ], col = MeanColors[i])
    }
  }


}




