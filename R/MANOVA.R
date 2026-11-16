MANOVA <- function(Y, Group, C = NULL, M = NULL, Effects = NULL,
                   InitialTransform = 5, AddOnes = FALSE, Contrasts = TRUE) {
  if (!is.factor(Group)) Group <- as.factor(Group)

  ContinuousDataTransform <- c(
    "Raw Data", "Substract the global mean", "Double centering",
    "Column centering", "Standardize columns", "Row centering",
    "Standardize rows", "Divide by the column means and center",
    "Normalized residuals from independence", "Divide by the range"
  )
  if (is.numeric(InitialTransform)) {
    InitialTransform <- ContinuousDataTransform[InitialTransform]
  }

  Bip <- list()
  if (is.null(rownames(Y))) {
    rownames(Y) <- rownames(Y, do.NULL = FALSE, prefix = "I")
  }
  RowNames <- rownames(Y)
  if (is.null(colnames(Y))) {
    colnames(Y) <- colnames(Y, do.NULL = FALSE, prefix = "V")
  }

  Bip$Title <- "MANOVA Biplot"
  Bip$Type <- "MANOVA"
  Bip$Non_Scaled_Data <- Y
  Bip$Means <- apply(Y, 2, mean)
  Bip$Medians <- apply(Y, 2, median)
  Bip$Deviations <- apply(Y, 2, sd)
  Bip$Minima <- apply(Y, 2, min)
  Bip$Maxima <- apply(Y, 2, max)
  Bip$P25 <- apply(Y, 2, quantile)[2, ]
  Bip$P75 <- apply(Y, 2, quantile)[4, ]
  Bip$GMean <- mean(as.matrix(Y))
  Bip$Initial_Transformation <- InitialTransform

  Y <- IniTransform(as.matrix(Y), transform = InitialTransform)
  rownames(Y) <- RowNames

  GroupNames <- levels(Group)
  L <- length(GroupNames)
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  r <- min(c(L - 1, J))
  X <- MultBiplotR::Factor2Binary(Group)
  colnames(X) <- GroupNames

  if (AddOnes) {
    X <- cbind(rep(1, I), X)
    rownames(X) <- RowNames
    colnames(X) <- c("Constante", GroupNames)
  }

  if (is.null(C)) {
    C <- diag(L)
    rownames(C) <- paste("C", GroupNames)
    colnames(C) <- GroupNames
  }
  if (is.null(M)) {
    M <- diag(J)
    colnames(M) <- colnames(Y)
  }
  VarNames <- colnames(M)
  k_m <- ncol(M)

  Bip$ncols <- J
  Bip$nrows <- I
  Bip$g <- L
  Bip$p <- J
  Bip$C <- C
  Bip$Combinaciones <- M
  DimNames <- paste("Dim", seq_len(r))

  Man <- MANOVA.Estimation(Y, X, C, M)

  Bip$MANOVA <- Man$MANOVA
  Bip$EigenValues <- Man$VP
  Bip$Inertia <- Man$Inercia
  Bip$CumInertia <- cumsum(Man$Inercia)

  Bip$GroupCoordinates <- Man$B %*% M %*% Man$A[, seq_len(r), drop = FALSE]
  rownames(Bip$GroupCoordinates) <- GroupNames
  colnames(Bip$GroupCoordinates) <- DimNames

  YM <- Y %*% M
  Bip$RowCoordinates <- YM %*% Man$A[, seq_len(r), drop = FALSE]
  rownames(Bip$RowCoordinates) <- RowNames
  colnames(Bip$RowCoordinates) <- DimNames

  Bip$ColCoordinates <- MultBiplotR::matrixsqrt(Man$E) %*% Man$V[, seq_len(r), drop = FALSE]
  rownames(Bip$ColCoordinates) <- VarNames
  colnames(Bip$ColCoordinates) <- DimNames

  Bip$Structure_Correlations <- cor(YM, Bip$RowCoordinates)
  rownames(Bip$Structure_Correlations) <- VarNames
  colnames(Bip$Structure_Correlations) <- DimNames

  Bip$Canonical_Weights <- Man$A[, seq_len(r), drop = FALSE]
  rownames(Bip$Canonical_Weights) <- VarNames
  colnames(Bip$Canonical_Weights) <- DimNames

  Bip$GroupContributions <- diag(1 / rowSums(Bip$GroupCoordinates^2), nrow = L) %*% Bip$GroupCoordinates^2
  Bip$ColContributions <- diag(1 / rowSums(Bip$ColCoordinates^2), nrow = k_m) %*% Bip$ColCoordinates^2

  Bip$ExplTotal <- matrix(0, r, 1)
  Bip$RowContributions <- matrix(0, I, r)
  Bip$QLRVars <- matrix(0, k_m, r)

  SCT <- sum(YM^2)
  SCRows <- rowSums(YM^2)
  SCCols <- colSums(YM^2)

  for (j in seq_len(r)) {
    Fitted <- Bip$RowCoordinates[, seq_len(j), drop = FALSE] %*% t(Bip$ColCoordinates[, seq_len(j), drop = FALSE])
    residuals <- YM - Fitted
    Bip$ExplTotal[j] <- 1 - sum(residuals^2) / SCT
    Bip$RowContributions[, j] <- 1 - rowSums(residuals^2) / SCRows
    Bip$QLRVars[, j] <- 1 - colSums(residuals^2) / SCCols
  }

  Bip$Clusters <- Group
  Bip$ClusterNames <- GroupNames
  Bip$ClusterType <- "us"
  Bip$ClusterColors <- rainbow(L)

  falfau <- qt(1 - (0.025), (I - L))
  falfab <- qt(1 - (0.025 / (I - L)), (I - L))
  falfam <- sqrt(qf(1 - 0.05, L, (I - L - J + 1)) * (((I - L) * J) / (I - L - J + 1)))
  falfac <- sqrt(qchisq(0.95, 2))
  IL <- as.numeric(table(Group))
  Bip$UnivRad <- falfau * (1 / sqrt(IL)) / sqrt(I - L)
  Bip$BonfRad <- falfab * (1 / sqrt(IL)) / sqrt(I - L)
  Bip$MultRad <- falfam * (1 / sqrt(IL)) / sqrt(I - L)
  Bip$ChisRad <- falfac * (1 / sqrt(IL)) / sqrt(I - L)

  nContrasts <- dim(C)[1]

  if (!is.null(Effects)) {
    if (!is.factor(Effects)) Effects <- as.factor(Effects)
    EstimEfectos <- list()
    efectos <- levels(Effects)
    nefec <- length(efectos)
    for (i in seq_len(nefec)) {
      cuales <- which(Effects == efectos[i])
      Man <- MANOVA.Estimation(Y, X, C[cuales, , drop = FALSE], M)
      EstimEfectos[[i]] <- Man$MANOVA
    }
    names(EstimEfectos) <- efectos
    Bip$Effects <- EstimEfectos
  }

  if (Contrasts) {
    EstimContrastes <- list()
    for (i in seq_len(nContrasts)) {
      Man <- MANOVA.Estimation(Y, X, C[i, , drop = FALSE], M)
      EstimContrastes[[i]] <- Man$MANOVA
    }
    names(EstimContrastes) <- rownames(C)
    Bip$Contrasts <- EstimContrastes
  }

  class(Bip) <- c("MANOVA", "MANOVA.Biplot")
  return(Bip)
}

