MANOVA.Estimation <- function(Y, X, C, M) {
  Man <- list()
  L <- dim(X)[2]
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  r <- min(c(L - 1, J))
  C <- matrix(C, ncol = L)

  XtX_inv <- ginv(t(X) %*% X)
  B <- XtX_inv %*% t(X) %*% Y
  O <- C %*% B %*% M
  R <- C %*% XtX_inv %*% t(C)

  E <- t(M) %*% t(Y) %*% (diag(I) - (X %*% XtX_inv %*% t(X))) %*% Y %*% M
  H <- t(O) %*% ginv(R) %*% O

  E2 <- MultBiplotR::matrixsqrtinv(E)
  YB2 <- E2 %*% H %*% E2

  dvs <- svd(YB2, 0)
  VP <- dvs$d[seq_len(r)]
  Inercia <- (VP / sum(VP)) * 100
  V <- dvs$v
  A <- E2 %*% V

  Man$B <- B
  Man$E <- E
  Man$H <- H
  Man$VP <- VP
  Man$Inercia <- Inercia
  Man$V <- V
  Man$A <- A

  Wilks <- prod(1 / (1 + VP))
  glh <- rankMatrix(H)[1]
  J_eff <- rankMatrix(H + E)[1]

  gle <- I - L
  denom_t <- J_eff^2 + glh^2 - 5
  t_stat <- if (denom_t <= 0) 1 else sqrt((glh^2 * J_eff^2 - 4) / denom_t)
  w <- gle + glh - 0.5 * (J_eff + glh + 1)
  df1w <- J_eff * glh
  df2w <- w * t_stat - 0.5 * (J_eff * glh - 2)
  s <- min(c(J_eff, glh))
  Wilksf <- ((1 - Wilks^(1 / t_stat)) / (Wilks^(1 / t_stat))) * (df2w / df1w)
  Wilksp <- round(1 - pf(Wilksf, df1w, df2w), digits = 5)
  Wilsr <- c(Wilks, Wilksf, df1w, df2w, Wilksp)

  Pillai <- sum(VP / (1 + VP))
  m <- (abs(J_eff - glh) - 1) / 2
  n <- (gle - J_eff - 1) / 2
  Pillaif <- ((2 * n + s + 1) * Pillai) / ((2 * m + s + 1) * (s - Pillai))
  df1p <- (2 * m + s + 1) * s
  df2p <- (2 * n + s + 1) * s
  Pillaip <- round(1 - pf(Pillaif, df1p, df2p), digits = 5)
  Pillair <- c(Pillai, Pillaif, df1p, df2p, Pillaip)

  Hotteling <- sum(VP)
  Hottelingf <- Hotteling * ((2 * (s * n + 1)) / (s^2 * (2 * m + s + 1)))
  df1h <- s * (2 * m + s + 1)
  df2h <- 2 * (s * n + 1)
  Hottelingp <- round(1 - pf(Hottelingf, df1h, df2h), digits = 5)
  Hottelingr <- c(Hotteling, Hottelingf, df1h, df2h, Hottelingp)

  Roy <- VP[1]
  v1 <- (abs(J_eff - glh) - 1) / 2
  v2 <- (gle - J_eff - 1) / 2
  Royf <- Roy * (2 * v2 + 2) / (2 * v1 + 2)
  df1r <- 2 * v1 + 2
  df2r <- 2 * v2 + 2
  Royp <- round(1 - pf(Royf, df1r, df2r), digits = 5)
  Royr <- c(Roy, Royf, df1r, df2r, Royp)

  ManovaR <- rbind(Wilsr, Pillair, Hottelingr, Royr)
  rownames(ManovaR) <- c("Wilk's lambda", "Pillai Trace", "Lawley-Hotteling", "Roy's greatest root")
  colnames(ManovaR) <- c("Estadistico", "F aprox", "df Num", "df Den", "p-valor")
  Man$MANOVA <- ManovaR
  return(Man)
}
