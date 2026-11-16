PERMANOVA.Estimation <- function(D, X, C = NULL, Effects = NULL) {
  Result <- list()
  L <- dim(X)[2]
  I <- dim(D)[1]
  if (is.null(C)) C <- diag(L)
  if (!is.matrix(C)) C <- matrix(C, ncol = L)
  nc <- dim(C)[1]
  if (!is.null(Effects) && !is.factor(Effects)) Effects <- as.factor(Effects)
  ne <- length(levels(Effects))

  S <- min(c(rankMatrix(C)[1]), L - 1)
  G <- (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I) %*% (-0.5 * D^2) %*%
    (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I)
  H <- X %*% ginv(t(X) %*% X) %*% t(X)
  R <- C %*% ginv(t(X) %*% X) %*% t(C)
  R12 <- MultBiplotR::matrixsqrtinv(R)
  A <- R12 %*% C %*% ginv(t(X) %*% X) %*% t(X)
  SCET <- sum(diag(H %*% G %*% H))
  SCEC <- diag(A %*% G %*% t(A))
  SCR <- sum(diag((diag(I) - H) %*% G %*% (diag(I) - H)))
  Fexp <- (SCET / S) / (SCR / (I - L))
  Result$Global <- cbind(SCET, SCR, S, I - L, Fexp)

  FCexp <- SCEC / (SCR / (I - L))
  SCEC <- matrix(SCEC, ncol = 1)

  Result$Contrastes <- cbind(SCEC, rep(SCR, nc), rep(1, nc), rep(I - L, nc), FCexp)
  rownames(Result$Contrastes) <- rownames(C)

  if (!is.null(Effects)) {
    XX <- t(MultBiplotR::Factor2Binary(Effects))
    tam <- XX %*% t(XX)
    SCE <- XX %*% SCEC
    FEexp <- (solve(tam) %*% XX %*% SCEC) / (SCR / (I - L))
    rownames(FEexp) <- levels(Effects)
    Result$Effects <- cbind(SCE, rep(SCR, ne), as.numeric(diag(tam)), rep(I - L, ne), FEexp)
    colnames(Result$Effects) <- c("Explained", "Residual", "df Num", "df Denom", "F-exp")
    rownames(Result$Effects) <- levels(Effects)
  }

  return(Result)
}
