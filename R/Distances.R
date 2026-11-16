Distances <- function(x, y = NULL, type = c("auto", "continuous", "binary", "categorical", "mixed"),
                      coef = NULL, ...) {
  type <- match.arg(type)

  if (type == "auto") {
    df_x <- as.data.frame(x)
    is_bin <- vapply(df_x, function(col) {
      u <- unique(stats::na.omit(col))
      length(u) <= 2 && all(u %in% c(0, 1, TRUE, FALSE))
    }, logical(1))
    is_cat <- vapply(df_x, function(col) is.factor(col) || is.character(col), logical(1)) & !is_bin
    is_cont <- vapply(df_x, is.numeric, logical(1)) & !is_bin

    if (all(is_bin)) {
      type <- "binary"
    } else if (all(is_cont)) {
      type <- "continuous"
    } else if (all(is_cat)) {
      type <- "categorical"
    } else {
      type <- "mixed"
    }
  }

  switch(type,
         continuous = DistContinuous(x = x, y = y,
                                     coef = if (is.null(coef)) "Pythagorean" else coef, ...),
         binary     = DistBinary(x = x, y = y,
                                 coefficient = if (is.null(coef)) "Simple_Matching" else coef, ...),
         categorical = {
           if (!is.null(y)) stop("y must be NULL for categorical distances.")
           DistCategorical(x = x, coef = if (is.null(coef)) "cramer" else coef, ...)
         },
         mixed = {
           if (!is.null(y)) stop("y must be NULL for mixed distances.")
           DistMixed(x = x, coef = if (is.null(coef)) "gower" else coef, ...)
         }
  )
}
