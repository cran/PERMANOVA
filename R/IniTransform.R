IniTransform <- function(X, InitTransform = "None", transform = "Standardize columns") {
  n <- nrow(X)
  p <- ncol(X)
  RowNames <- rownames(X)
  ColNames <- colnames(X)

  InitTransforms <- c("None", "Log", "Logit")
  if (is.numeric(InitTransform)) {
    InitTransform <- InitTransforms[InitTransform]
  }
  InitTransform <- match.arg(InitTransform, InitTransforms)

  ContinuousDataTransform <- c(
    "Raw Data", "Substract the global mean", "Double centering",
    "Column centering", "Standardize columns", "Row centering",
    "Standardize rows", "Divide by the column means and center",
    "Normalized residuals from independence", "Divide by the range"
  )
  if (is.numeric(transform)) {
    transform <- ContinuousDataTransform[transform]
  }
  transform <- match.arg(transform, ContinuousDataTransform)

  switch(InitTransform,
         Log = {
           if (any(X <= 0, na.rm = TRUE)) {
             stop("Initial log transformation is not compatible with negative or zero values")
           }
           X <- log(X)
         },
         Logit = {
           if (any(X < 0 | X > 1, na.rm = TRUE)) {
             stop("Initial logit transformation requires values between 0 and 1")
           }
           X <- X + 0.01 * (X == 0) - 0.01 * (X == 1)
           X <- log(X / (1 - X))
         }
  )

  switch(transform,
         `Substract the global mean` = {
           gmean <- mean(X)
           X <- X - gmean
         },
         `Double centering` = {
           X <- (diag(n) - matrix(1, n, n) / n) %*% X %*% (diag(p) - matrix(1, p, p) / p)
         },
         `Column centering` = {
           means <- apply(X, 2, mean)
           X <- X - matrix(1, n, 1) %*% matrix(means, 1, p)
         },
         `Standardize columns` = {
           means <- apply(X, 2, mean)
           stdDevs <- apply(X, 2, sd)
           X <- (X - matrix(1, n, 1) %*% matrix(means, 1, p)) /
             (matrix(1, n, 1) %*% matrix(stdDevs, 1, p))
         },
         `Row centering` = {
           X <- X %*% (diag(p) - matrix(1, p, p) / p)
         },
         `Standardize rows` = {
           stdDevs <- apply(X, 1, sd)
           X <- solve(diag(stdDevs)) %*% X %*% (diag(p) - matrix(1, p, p) / p)
         },
         `Divide by the column means and center` = {
           means <- apply(X, 2, mean)
           for (i in seq_len(p)) X[, i] <- X[, i] / means[i]
           X <- (diag(n) - matrix(1, n, n) / n) %*% X
         },
         `Normalized residuals from independence` = {
           nt <- sum(X)
           dr <- rowSums(X)
           dc <- colSums(X)
           esp <- (matrix(dr, ncol = 1) %*% matrix(dc, nrow = 1)) / nt
           var_res <- (matrix(1 - dr / nt, ncol = 1)) %*% matrix(1 - dc / nt, nrow = 1)
           X <- ((X - esp) / sqrt(esp)) / sqrt(var_res)
         },
         `Divide by the range` = {
           Rangos <- apply(X, 2, max) - apply(X, 2, min)
           X <- X %*% diag(1 / Rangos, nrow = p)
         }
  )

  rownames(X) <- RowNames
  colnames(X) <- ColNames
  return(X)
}
