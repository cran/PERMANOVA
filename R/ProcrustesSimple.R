ProcrustesSimple <- function (X, Y, centre = FALSE)
{
  nx = nrow(X)
  px = ncol(X)
  ny = nrow(Y)
  py = ncol(Y)
  if (nx != ny)
    stop("Matrices should have the same number of rows: ",
         nrow(X), " or ", nrow(Y))
  if (px != py) {
    warning("Matrices should have the same number of columns: Adjusting to conform")
    add = abs(px - py)
    if (px < py)
      X = cbind(X, matrix(0, nx, add))
    if (px > py)
      Y = cbind(Y, matrix(0, ny, add))
    px = ncol(X)
    py = ncol(Y)
  }
  nx = nrow(X)
  px = ncol(X)
  ny = nrow(Y)
  py = ncol(Y)
  J = diag(nx) - matrix(1, nx, nx)/nx
  if (centre) {
    X = J %*% X
    Y = J %*% Y
  }
  Procrustes = list()
  Procrustes$xmean = apply(X, 2, mean)
  Procrustes$X = X
  Procrustes$Y = Y
  C = t(X) %*% Y
  SVD = svd(C)
  TT = SVD$u %*% t(SVD$v)
  s=tr(TT %*% t(Y) %*% X)/tr(t(Y) %*% Y)
  t = apply((X - s * Y %*% TT), 2, sum)/ny
  Z = s * Y %*% TT + matrix(1, nx, 1) %*% t
  Procrustes$Yrot = Z
  Procrustes$rotation = TT
  Procrustes$translation = t
  Procrustes$scale = s
  Procrustes$rss = sum((X - Z)^2)
  Procrustes$fit = 1 - sum((X - Z)^2)/sum(X^2)
  Procrustes$correlations = cor(X, Z)
  class(Procrustes) = "Procrustes"
  return(Procrustes)
}
