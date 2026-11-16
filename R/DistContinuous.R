DistContinuous <- function(x, y = NULL, coef = "Pythagorean", r = 1,
                           alpha = 0.1, euclidean = FALSE) {
  distances <- c("Pythagorean", "Taxonomic", "City", "Minkowski", "Divergence",
                 "dif_sum", "Camberra", "Bray_Curtis", "Soergel", "Ware_Hedges",
                 "Gower", "Euclidean_Standardized", "Maximum", "Cosine",
                 "Correlation", "Mahalanobis", "Robust_GV", "Robust_MCD")
  if (is.numeric(coef)) coef <- distances[coef]
  coef <- match.arg(coef, distances)

  if (!is.matrix(x)) x <- as.matrix(x)
  if (is.null(y)) {
    same_xy <- TRUE
    y <- x
  } else {
    same_xy <- identical(x, y)
    if (!is.matrix(y)) y <- as.matrix(y)
  }

  n <- nrow(x)
  p <- ncol(x)
  s <- nrow(y)
  q <- ncol(y)
  if (p != q) stop("The matrices should have the same number of columns")

  NamesX <- rownames(x)
  NamesY <- rownames(y)

  dbrobust_methods <- c("Euclidean_Standardized", "Maximum", "Cosine",
                        "Correlation", "Mahalanobis", "Robust_GV", "Robust_MCD")

  if (coef %in% dbrobust_methods) {
    if (!same_xy) {
      stop("The selected dbrobust coefficient currently supports only a single data matrix (y = NULL).")
    }
    if (coef %in% c("Robust_GV", "Robust_MCD")) {
      cov_m <- if (coef == "Robust_GV") "gv" else "mcd"
      df_x <- as.data.frame(x)
      if (is.null(colnames(df_x))) colnames(df_x) <- paste0("V", seq_len(p))
      d_sq <- dbrobust::robust_distances(
        data = df_x,
        cont_vars = colnames(df_x),
        method = "ggower",
        cov_method = cov_m,
        alpha = alpha,
        return_dist = FALSE
      )
      dis <- sqrt(pmax(d_sq, 0))
    } else {
      method_map <- c(
        Euclidean_Standardized = "euclidean_standardized",
        Maximum = "maximum",
        Cosine = "cosine",
        Correlation = "correlation",
        Mahalanobis = "mahalanobis"
      )
      dis <- dbrobust::calculate_distances(
        x,
        method = method_map[[coef]],
        output_format = "matrix"
      )
    }
  } else if (coef == "Pythagorean" && same_xy) {
    dis <- as.matrix(stats::dist(x, method = "euclidean"))
  } else if (coef == "City" && same_xy) {
    dis <- as.matrix(stats::dist(x, method = "manhattan"))
  } else if (coef == "Minkowski" && same_xy) {
    dis <- as.matrix(stats::dist(x, method = "minkowski", p = r))
  } else {
    if (coef == "Gower") {
      comb <- rbind(y, x)
      rank <- apply(comb, 2, max) - apply(comb, 2, min)
      rank[rank == 0] <- 1
    }
    dis <- matrix(0, s, n)
    for (i in seq_len(s)) {
      for (j in seq_len(n)) {
        yi <- y[i, ]
        xj <- x[j, ]
        switch(coef,
               Pythagorean = { dis[i, j] <- sqrt(sum((yi - xj)^2)) },
               Taxonomic   = { dis[i, j] <- sqrt(sum(((yi - xj)^2) / r^2)) },
               City        = { dis[i, j] <- sum(abs(yi - xj)) },
               Minkowski   = { dis[i, j] <- (sum(abs(yi - xj)^r))^(1 / r) },
               Divergence  = {
                 den <- (yi + xj)^2
                 dis[i, j] <- sqrt(sum(ifelse(den == 0, 0, (yi - xj)^2 / den)))
               },
               dif_sum     = {
                 den <- abs(yi + xj)
                 dis[i, j] <- sum(ifelse(den == 0, 0, abs(yi - xj) / den))
               },
               Camberra    = {
                 den <- abs(yi) + abs(xj)
                 dis[i, j] <- sum(ifelse(den == 0, 0, abs(yi - xj) / den))
               },
               Bray_Curtis = {
                 den <- sum(yi + xj)
                 dis[i, j] <- if (den == 0) 0 else sum(abs(yi - xj)) / den
               },
               Soergel     = {
                 den <- sum(pmax(yi, xj))
                 dis[i, j] <- if (den == 0) 0 else sum(abs(yi - xj)) / den
               },
               Ware_Hedges = {
                 den <- pmax(yi, xj)
                 dis[i, j] <- sum(ifelse(den == 0, 0, 1 - pmin(yi, xj) / den))
               },
               Gower       = { dis[i, j] <- sum(abs(yi - xj) / rank) }
        )
      }
    }
  }

  if (euclidean && same_xy) {
    euc_res <- dbrobust::make_euclidean(dis^2)
    dis <- sqrt(pmax(euc_res$D_euc, 0))
  }

  rownames(dis) <- NamesY
  colnames(dis) <- NamesX
  distances <- list(Data = x, D = dis, Coefficient = coef, Type = "continuous")
  class(distances) <- "proximities"
  return(distances)
}



