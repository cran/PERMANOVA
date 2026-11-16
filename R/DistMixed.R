DistMixed <- function(x, coef = "gower", cont_vars = NULL, bin_vars = NULL,
                      cat_vars = NULL, w = NULL, type = list(),
                      cov_method = "gv", alpha = 0.1, euclidean = FALSE) {
  coefficients <- c("gower", "ggower", "relms")
  if (is.numeric(coef)) coef <- coefficients[coef]
  coef <- match.arg(coef, coefficients)

  df_x <- as.data.frame(x)
  if (is.null(colnames(df_x))) colnames(df_x) <- paste0("V", seq_len(ncol(df_x)))

  if (coef == "gower") {
    dis <- dbrobust::calculate_distances(
      df_x,
      method = "gower",
      output_format = "matrix",
      type = type,
      weights = w
    )
  } else {
    if (is.null(cont_vars) && is.null(bin_vars) && is.null(cat_vars)) {
      is_bin <- vapply(df_x, function(col) {
        u <- unique(stats::na.omit(col))
        length(u) <= 2 && all(u %in% c(0, 1, TRUE, FALSE))
      }, logical(1))
      is_cat <- vapply(df_x, function(col) is.factor(col) || is.character(col), logical(1)) & !is_bin
      is_cont <- !is_bin & !is_cat

      cont_vars <- names(df_x)[is_cont]
      bin_vars  <- names(df_x)[is_bin]
      cat_vars  <- names(df_x)[is_cat]
      if (length(cont_vars) == 0) cont_vars <- NULL
      if (length(bin_vars) == 0) bin_vars <- NULL
      if (length(cat_vars) == 0) cat_vars <- NULL
    }

    d_sq <- dbrobust::robust_distances(
      data = df_x,
      cont_vars = cont_vars,
      bin_vars = bin_vars,
      cat_vars = cat_vars,
      w = w,
      method = coef,
      cov_method = cov_method,
      alpha = alpha,
      return_dist = FALSE
    )
    dis <- sqrt(pmax(d_sq, 0))
  }

  if (euclidean) {
    euc_res <- dbrobust::make_euclidean(dis^2)
    dis <- sqrt(pmax(euc_res$D_euc, 0))
  }

  rownames(dis) <- rownames(df_x)
  colnames(dis) <- rownames(df_x)
  distances <- list(Data = df_x, D = dis, Coefficient = coef, Type = "mixed")
  class(distances) <- "proximities"
  return(distances)
}
