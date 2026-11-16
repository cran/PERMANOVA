BootDisMANOVA <- function(Distance, groups, C = NULL, Effects = NULL, nB = 1000,
                          seed = NULL, CoordPrinc = FALSE, dimens = 2,
                          PCoA = "Standard", ProjectInd = TRUE, tol = 1e-4,
                          DatosIni = TRUE) {
  D <- Distance$D
  cl <- match.call()

  PCoAs <- c("Standard", "Weighted")
  if (is.numeric(PCoA)) PCoA <- PCoAs[PCoA]
  PCoA <- match.arg(PCoA, PCoAs)

  if (!is.factor(groups)) stop("The grouping variable must be a factor")

  if (!is.null(seed)) set.seed(seed)
  Perm <- list()
  Perm$call <- cl

  if (is.null(rownames(D))) {
    rownames(D) <- rownames(D, do.NULL = FALSE, prefix = "I")
  }

  Perm$Title <- "MANOVA BASED ON DISTANCES"
  Perm$Type <- "BootMANOVA"
  Perm$Distances <- D

  GroupNames <- levels(groups)
  L <- length(GroupNames)
  I <- dim(D)[1]
  X <- MultBiplotR::Factor2Binary(groups)
  colnames(X) <- GroupNames
  rownames(X) <- rownames(D)

  if (is.null(C)) {
    C <- diag(L)
    rownames(C) <- paste("C", GroupNames)
    colnames(C) <- GroupNames
  }

  nc <- dim(C)[1]
  Perm$C <- C

  if (!is.null(Effects)) {
    Perm$Effects <- Effects
  }

  # PERMANOVA INICIAL
  Perm$Initial <- PERMANOVA.Estimation(D, X, C, Effects)

  if (is.null(seed)) set.seed(1) else set.seed(seed)

  # Bootstrap
  Ftotal <- matrix(0, 1, nB)
  FContrastes <- matrix(0, nc, nB)
  if (!is.null(Effects)) {
    ne <- length(levels(Effects))
    FEffects <- matrix(0, ne, nB)
  }

  for (i in 1:nB) {
    muestra <- sample.int(I, size = I, replace = TRUE)
    Man <- PERMANOVA.Estimation(D[muestra, muestra], X, C, Effects)
    Ftotal[i] <- Man$Global[5]
    FContrastes[, i] <- Man$Contrastes[, 5]
    if (!is.null(Effects)) {
      FEffects[, i] <- Man$Effects[, 5]
    }
  }

  Perm$DistMuestral <- Ftotal
  Perm$pvalue <- (sum(Ftotal >= Perm$Initial$Global[1, 5]) + 1) / (nB + 1)
  Perm$Initial$Global <- cbind(Perm$Initial$Global, Perm$pvalue)
  colnames(Perm$Initial$Global) <- c("Explained", "Residual", "G.L. Num", "G.L. Denom", "F-exp", "p-value")
  rownames(Perm$Initial$Global) <- "Total"

  Perm$DistMuestralC <- FContrastes
  pvalC <- (rowSums(FContrastes >= matrix(Perm$Initial$Contrastes[, 5], nrow = nc, ncol = nB)) + 1) / (nB + 1)
  Perm$Initial$Contrastes <- cbind(Perm$Initial$Contrastes, pvalC)
  colnames(Perm$Initial$Contrastes) <- c("Explained", "Residual", "G.L. Num", "G.L. Denom", "F-exp", "p-value")

  if (!is.null(Effects)) {
    Perm$DistMuestralEffects <- FEffects
    pvalE <- (rowSums(FEffects >= matrix(Perm$Initial$Effects[, 5], nrow = ne, ncol = nB)) + 1) / (nB + 1)
    Perm$Initial$Effects <- cbind(Perm$Initial$Effects, pvalE)
    colnames(Perm$Initial$Effects) <- c("Explained", "Residual", "G.L. Num", "G.L. Denom", "F-exp", "p-value")
  }

  Perm$CoordPrinc <- CoordPrinc

  if (CoordPrinc) {
    N <- t(X) %*% X
    invN <- solve(N)
    FS <- invN %*% (t(X) %*% (0.5 * D^2) %*% X) %*% invN
    f <- matrix(diag(FS), L, 1)

    DB <- 2 * FS - f %*% matrix(1, 1, L) - t(f %*% matrix(1, 1, L))
    DB <- 0.5 * DB

    switch(PCoA,
           Standard = {
             H <- diag(L) - matrix(1, L, L) / L
             B <- -H %*% DB %*% H
           },
           Weighted = {
             H <- diag(L) - matrix(1, L, 1) %*% matrix(diag(N), 1, L) / I
             B <- -H %*% DB %*% H
           }
    )

    solut <- eigen(B, symmetric = TRUE)
    pos <- which(solut$values > tol)
    if (length(pos) == 0) pos <- 1
    vals <- solut$values[pos]
    vecs <- solut$vectors[, pos, drop = FALSE]

    dimens_eff <- min(dimens, length(vals))
    vp <- vals[1:dimens_eff]
    Inertia <- (vp / sum(vals)) * 100
    Perm$ExplainedVariance <- Inertia
    Perm$Inertias <- cbind(vp, Inertia, cumsum(Inertia))
    rownames(Perm$Inertias) <- paste("PCo", 1:dimens_eff)
    colnames(Perm$Inertias) <- c("Eigenvalue", "Explained Variance", "Cumulative")

    Y <- vecs %*% diag(sqrt(vals), nrow = length(vals))
    d0 <- rowSums(Y^2)
    rownames(Y) <- GroupNames

    st <- d0
    Perm$MeanCoordinates <- Y[, 1:dimens_eff, drop = FALSE]
    colnames(Perm$MeanCoordinates) <- paste("Dim", 1:dimens_eff)
    qlr <- diag(1 / st, nrow = L) %*% (Perm$MeanCoordinates^2)
    Perm$Qualities <- round(qlr * 100, digits = 2)
    rownames(Perm$Qualities) <- GroupNames
    colnames(Perm$Qualities) <- paste("Dim", 1:dimens_eff)
    Perm$CummulativeQualities <- t(apply(Perm$Qualities, 1, cumsum))

    if (ProjectInd) {
      Y_sub <- Y[, 1:dimens_eff, drop = FALSE]
      cont_coefs <- c("euclidean", "manhattan", "minkowski", "canberra", "bray", "kulczynski", "gower")
      can_use_raw <- DatosIni && !is.null(Distance$Data) &&
        (is.matrix(Distance$Data) || all(sapply(Distance$Data, is.numeric))) &&
        (!is.null(Distance$Coefficient) && tolower(Distance$Coefficient) %in% cont_coefs)

      if (can_use_raw) {
        YY <- as.matrix(Distance$Data)
        Means <- invN %*% t(X) %*% YY
        Di <- DistContinuous(Means, y = YY, coef = tolower(Distance$Coefficient))$D^2
      } else {
        # Proyección exacta de Gower desde la matriz de distancias D^2
        d_il_mean <- (D^2) %*% X %*% invN
        d_ll_mean <- diag(invN %*% (t(X) %*% (D^2) %*% X) %*% invN)
        Di <- sweep(d_il_mean, 2, 0.5 * d_ll_mean, "-")
      }

      Yi <- -0.5 * MultBiplotR::ginv(t(Y_sub) %*% Y_sub) %*% t(Y_sub) %*% H %*%
        t(Di - matrix(1, I, 1) %*% matrix(d0, 1, L))
      Yi <- t(Yi)
      rownames(Yi) <- rownames(D)
      Perm$RowCoordinates <- Yi
    }

    Perm <- MultBiplotR::AddCluster2Biplot(Perm, ClusterType = "us", Groups = groups)
  }

  class(Perm) <- c("BootDisMANOVA")
  return(Perm)
}



