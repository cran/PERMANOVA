DistCategorical <- function(x, coef = "cramer", weights = NULL, euclidean = FALSE) {
  coefficients <- c("matching_coefficient", "cramer")
  if (is.numeric(coef)) coef <- coefficients[coef]
  coef <- match.arg(coef, coefficients)

  df_x <- as.data.frame(x)
  df_x[] <- lapply(df_x, as.factor)

  dis <- dbrobust::calculate_distances(
    df_x,
    method = coef,
    output_format = "matrix",
    weights = weights
  )

  if (euclidean) {
    euc_res <- dbrobust::make_euclidean(dis^2)
    dis <- sqrt(pmax(euc_res$D_euc, 0))
  }

  rownames(dis) <- rownames(df_x)
  colnames(dis) <- rownames(df_x)
  distances <- list(Data = df_x, D = dis, Coefficient = coef, Type = "categorical")
  class(distances) <- "proximities"
  return(distances)
}
