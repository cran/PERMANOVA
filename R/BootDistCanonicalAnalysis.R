BootDistCanonicalAnalysis <- function(Distance, groups, dimens = NULL, nB = 100,
                                      seed = NULL, PCoA = "Standard",
                                      ProcrustesRot = TRUE, DatosIni = TRUE,
                                      tol = 0.0001) {
  if (!is.factor(groups)) groups <- as.factor(groups)

  result <- list()
  D <- Distance$D

  if (is.null(rownames(D))) {
    rownames(D) <- rownames(D, do.NULL = FALSE, prefix = "I")
  }

  result$D <- D
  result$Coefficient <- Distance$Coefficient
  result$nB <- nB
  result$Groups <- groups
  GroupNames <- levels(groups)
  result$GroupNames <- GroupNames
  L <- length(GroupNames)
  I <- dim(D)[1]
  X <- MultBiplotR::Factor2Binary(groups)
  colnames(X) <- GroupNames
  rownames(X) <- rownames(D)

  N <- t(X) %*% X
  Ns <- diag(N)
  invN <- diag(1 / Ns, nrow = L)

  if (is.null(dimens)) dimens <- L - 1

  F_mat <- invN %*% (t(X) %*% (0.5 * D^2) %*% X) %*% invN

  switch(PCoA,
         Standard = {
           H <- diag(L) - matrix(1, L, L) / L
         },
         Weighted = {
           H <- diag(L) - matrix(1, L, 1) %*% matrix(Ns, 1, L) / I
         }
  )

  B <- -1 * H %*% F_mat %*% H
  solut <- eigen(B, symmetric = TRUE)
  pos <- which(solut$values > tol)
  if (length(pos) == 0) pos <- 1
  vals <- solut$values[pos]
  vecs <- solut$vectors[, pos, drop = FALSE]

  dimens <- min(dimens, length(vals))
  vp <- vals[1:dimens]
  Inertia <- (vp / sum(vals)) * 100
  result$Inertia <- Inertia
  Inertias <- cbind(vp, Inertia, cumsum(Inertia))
  rownames(Inertias) <- paste("PCo", 1:dimens)
  colnames(Inertias) <- c("Eigenvalue", "Explained Variance", "Cumulative")
  result$Inertias <- Inertias

  Y_full <- vecs %*% diag(sqrt(vals), nrow = length(vals))
  d0 <- rowSums(Y_full^2)
  Y <- Y_full[, 1:dimens, drop = FALSE]
  rownames(Y) <- GroupNames

  result$MeanCoordinates <- Y
  colnames(result$MeanCoordinates) <- paste("Dim", 1:dimens)
  qlr <- diag(1 / d0, nrow = L) %*% (result$MeanCoordinates^2)
  result$Qualities <- round(qlr * 100, digits = 2)
  rownames(result$Qualities) <- GroupNames
  colnames(result$Qualities) <- paste("Dim", 1:dimens)
  result$CummulativeQualities <- t(apply(result$Qualities, 1, cumsum))

  if (!is.null(seed)) set.seed(seed)
  coord <- array(0, c(L, dimens, nB))
  idx_groups <- split(seq_len(I), groups)

  for (i in 1:nB) {
    # Remuestreo estratificado independiente del orden de filas
    muestra <- unlist(lapply(idx_groups, function(idx) {
      idx[sample.int(length(idx), size = length(idx), replace = TRUE)]
    }), use.names = FALSE)

    D2 <- D[muestra, muestra]
    X_boot <- X[muestra, , drop = FALSE]

    F2 <- invN %*% (t(X_boot) %*% (0.5 * D2^2) %*% X_boot) %*% invN
    B2 <- -1 * H %*% F2 %*% H
    solut2 <- eigen(B2, symmetric = TRUE)
    vals2 <- pmax(solut2$values, 0)
    k2 <- min(dimens, length(vals2))
    Y2 <- matrix(0, L, dimens)
    Y2[, 1:k2] <- solut2$vectors[, 1:k2, drop = FALSE] %*%
      diag(sqrt(vals2[1:k2]), nrow = k2)

    if (ProcrustesRot) {
      Y2 <- ProcrustesSimple(Y, Y2, centre = TRUE)$Yrot
    }
    coord[, , i] <- Y2
  }
  result$CoordBoot <- coord

  CoordinatesMeans <- result$MeanCoordinates
  for (i in 1:L) {
    CoordinatesMeans[i, ] <- rowMeans(matrix(coord[i, , ], nrow = dimens, ncol = nB))
  }
  result$BootstrapCenters <- CoordinatesMeans

  # Proyección de individuos
  cont_coefs <- c("euclidean", "manhattan", "minkowski", "canberra", "bray", "kulczynski", "gower")
  can_use_raw <- DatosIni && !is.null(Distance$Data) &&
    (is.matrix(Distance$Data) || all(sapply(Distance$Data, is.numeric))) &&
    (!is.null(Distance$Coefficient) && tolower(Distance$Coefficient) %in% cont_coefs)

  if (can_use_raw) {
    YY <- as.matrix(Distance$Data)
    Means <- invN %*% t(X) %*% YY
    Di <- DistContinuous(Means, y = YY, coef = tolower(Distance$Coefficient))$D^2
  } else {
    d_il_mean <- (D^2) %*% X %*% invN
    d_ll_mean <- diag(invN %*% (t(X) %*% (D^2) %*% X) %*% invN)
    Di <- sweep(d_il_mean, 2, 0.5 * d_ll_mean, "-")
  }

  Yi <- -0.5 * MultBiplotR::ginv(t(Y) %*% Y) %*% t(Y) %*% H %*%
    t(Di - matrix(1, I, 1) %*% matrix(d0, 1, L))
  Yi <- t(Yi)
  rownames(Yi) <- rownames(D)
  result$RowCoordinates <- Yi

  class(result) <- "BootCanonAnalysis"
  return(result)
}


