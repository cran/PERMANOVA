PERMANOVA <- function(Distance, group, C = NULL, Effects = NULL, nperm = 1000,
                      seed = NULL, CoordPrinc = FALSE, dimens = 2, PCoA = "Standard",
                      ProjectInd = TRUE, tol = 1e-4, DatosIni = TRUE,
                      PostHoc = "bonferroni", parallel = FALSE, cores = NULL) {

  D <- Distance$D
  cl <- match.call()

  PCoAs <- c("Standard", "Weighted")
  if (is.numeric(PCoA)) PCoA <- PCoAs[PCoA]
  PCoA <- match.arg(PCoA, PCoAs)

  if (!is.factor(group)) stop("The grouping variable must be a factor")

  if (!is.null(seed)) set.seed(seed)
  Perm <- list()
  Perm$call <- cl

  if (is.null(rownames(D))) {
    rownames(D) <- rownames(D, do.NULL = FALSE, prefix = "I")
  }

  Perm$Title <- "MANOVA based on PERMUTATION and DISTANCES"
  Perm$Type <- "PERMANOVA"
  Perm$Distances <- D

  GroupNames <- levels(group)
  L <- length(GroupNames)
  I <- dim(D)[1]
  X <- MultBiplotR::Factor2Binary(group)
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
    if (!is.factor(Effects)) Effects <- as.factor(Effects)
    Perm$Effects <- Effects
  }

  # PERMANOVA INICIAL
  Perm$Initial <- PERMANOVA.Estimation(D, X, C, Effects)

  # Permutaciones
  Ftotal <- matrix(0, 1, nperm)
  FContrastes <- matrix(0, nc, nperm)
  if (!is.null(Effects)) {
    ne <- length(levels(Effects))
    FEffects <- matrix(0, ne, nperm)
  }

  seed_used <- if (is.null(seed)) 1 else seed
  set.seed(seed_used)

  if (parallel) {
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (is.null(cores)) {
      if (nzchar(chk) && chk == "TRUE") {
        cores <- 2L
      } else {
        avail <- parallel::detectCores(logical = FALSE)
        if (is.na(avail)) avail <- 2L
        cores <- max(1L, min(avail - 1L, 2L))
      }
    } else if (nzchar(chk) && chk == "TRUE") {
      cores <- min(as.integer(cores), 2L)
    }

    if (cores > 1L) {
      cl_par <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl_par), add = TRUE)
      parallel::clusterSetRNGStream(cl_par, iseed = seed_used)

      FPer <- parallel::parLapply(
        cl_par,
        seq_len(nperm),
        function(i, D_mat, X_mat, C_mat, Eff, n_ind, est_fn) {
          muestra <- sample.int(n_ind)
          est_fn(D_mat[muestra, muestra], X_mat, C_mat, Eff)
        },
        D_mat = D, X_mat = X, C_mat = C, Eff = Effects,
        n_ind = I, est_fn = PERMANOVA.Estimation
      )

      for (i in seq_len(nperm)) {
        Ftotal[i] <- FPer[[i]]$Global[5]
        FContrastes[, i] <- FPer[[i]]$Contrastes[, 5]
        if (!is.null(Effects)) {
          FEffects[, i] <- FPer[[i]]$Effects[, 5]
        }
      }
    } else {
      parallel <- FALSE
    }
  }

  if (!parallel) {
    for (i in seq_len(nperm)) {
      muestra <- sample.int(I)
      Man <- PERMANOVA.Estimation(D[muestra, muestra], X, C, Effects)
      Ftotal[i] <- Man$Global[5]
      FContrastes[, i] <- Man$Contrastes[, 5]
      if (!is.null(Effects)) {
        FEffects[, i] <- Man$Effects[, 5]
      }
    }
  }

  Perm$DistMuestral <- Ftotal
  Perm$pvalue <- (sum(Ftotal >= Perm$Initial$Global[1, 5]) + 1) / (nperm + 1)
  Perm$Initial$Global <- cbind(Perm$Initial$Global, Perm$pvalue, Perm$pvalue)

  colnames(Perm$Initial$Global) <- c("Explained", "Residual", "df Num", "df Denom", "F-exp", "p-value", "p-value adj.")
  rownames(Perm$Initial$Global) <- "Total"

  Perm$DistMuestralC <- FContrastes
  pvalC <- (rowSums(FContrastes >= matrix(Perm$Initial$Contrastes[, 5], nrow = nc, ncol = nperm)) + 1) / (nperm + 1)
  pvalCadjust <- round(stats::p.adjust(pvalC, PostHoc), 3)
  Perm$Initial$Contrastes <- cbind(Perm$Initial$Contrastes, pvalC, pvalCadjust)
  colnames(Perm$Initial$Contrastes) <- c("Explained", "Residual", "df Num", "df Denom", "F-exp", "p-value", "p-value adj.")

  if (!is.null(Effects)) {
    Perm$DistMuestralEffects <- FEffects
    pvalE <- (rowSums(FEffects >= matrix(Perm$Initial$Effects[, 5], nrow = ne, ncol = nperm)) + 1) / (nperm + 1)
    pvalEadjust <- round(stats::p.adjust(pvalE, PostHoc), 3)
    Perm$Initial$Effects <- cbind(Perm$Initial$Effects, pvalE, pvalEadjust)
    colnames(Perm$Initial$Effects) <- c("Explained", "Residual", "df Num", "df Denom", "F-exp", "p-value", "p-value adj.")
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

    Perm$MeanCoordinates <- Y[, 1:dimens_eff, drop = FALSE]
    colnames(Perm$MeanCoordinates) <- paste("Dim", 1:dimens_eff)
    qlr <- diag(1 / d0, nrow = L) %*% (Perm$MeanCoordinates^2)
    Perm$Qualities <- round(qlr * 100, digits = 2)
    rownames(Perm$Qualities) <- GroupNames
    colnames(Perm$Qualities) <- paste("Dim", 1:dimens_eff)
    Perm$CummulativeQualities <- t(apply(Perm$Qualities, 1, cumsum))

    if (ProjectInd) {
      Y_sub <- Y[, 1:dimens_eff, drop = FALSE]
      cont_coefs <- c("Pythagorean", "Taxonomic", "City", "Minkowski", "Divergence",
                      "dif_sum", "Camberra", "Bray_Curtis", "Soergel", "Ware_Hedges", "Gower")
      can_use_datos_ini <- DatosIni && !is.null(Distance$Data) &&
        (is.matrix(Distance$Data) || all(vapply(Distance$Data, is.numeric, logical(1)))) &&
        !is.null(Distance$Coefficient) && (Distance$Coefficient %in% cont_coefs)

      if (can_use_datos_ini) {
        YY <- as.matrix(Distance$Data)
        Means <- invN %*% t(X) %*% YY
        Di <- DistContinuous(Means, y = YY, coef = Distance$Coefficient)$D^2
      } else {
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

    Perm <- MultBiplotR::AddCluster2Biplot(Perm, ClusterType = "us", Groups = group)
  }

  class(Perm) <- c("PERMANOVA")
  return(Perm)
}




