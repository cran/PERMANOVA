ProcrustesSimple <- function(X, Y, centre = FALSE) {
  nx <- nrow(X)
  px <- ncol(X)
  ny <- nrow(Y)
  py <- ncol(Y)
  if (nx != ny) {
    stop("Matrices should have the same number of rows: ", nrow(X), " or ", nrow(Y))
  }
  if (px != py) {
    warning("Matrices should have the same number of columns: Adjusting to conform")
    add <- abs(px - py)
    if (px < py) X <- cbind(X, matrix(0, nx, add))
    if (px > py) Y <- cbind(Y, matrix(0, ny, add))
    px <- ncol(X)
    py <- ncol(Y)
  }

  J <- diag(nx) - matrix(1, nx, nx) / nx
  if (centre) {
    X <- J %*% X
    Y <- J %*% Y
  }

  Procrustes <- list()
  Procrustes$xmean <- apply(X, 2, mean)
  Procrustes$X <- X
  Procrustes$Y <- Y
  C <- t(X) %*% Y
  SVD <- svd(C)
  TT <- SVD$u %*% t(SVD$v)

  denom_s <- tr(t(Y) %*% Y)
  s <- if (denom_s == 0) 1 else tr(TT %*% t(Y) %*% X) / denom_s
  t_vec <- colSums(X - s * Y %*% TT) / ny
  Z <- s * Y %*% TT + matrix(1, nx, 1) %*% matrix(t_vec, nrow = 1)

  Procrustes$Yrot <- Z
  Procrustes$rotation <- TT
  Procrustes$translation <- t_vec
  Procrustes$scale <- s
  Procrustes$rss <- sum((X - Z)^2)
  denom_x <- sum(X^2)
  Procrustes$fit <- if (denom_x == 0) NA_real_ else 1 - Procrustes$rss / denom_x
  Procrustes$correlations <- cor(X, Z)
  class(Procrustes) <- "Procrustes"
  return(Procrustes)
}

