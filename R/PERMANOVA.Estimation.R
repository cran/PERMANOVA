PERMANOVA.Estimation <- function(D, X, C, Effects=NULL){
  Result=list()
  L = dim(X)[2] #Número de grupos
  I = dim(D)[1] #Número de de individuos
  if (!is.matrix(C)) C=matrix(C, ncol=L)
  nc=dim(C)[1]
  ne=length(levels(Effects))

  S=min(c(rankMatrix(C)[1]),L-1)
  G <- (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I) %*% (-0.5 * D^2) %*% (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I)
  H=X %*% Ginv(t(X) %*% X) %*% t(X)
  R=C %*% Ginv(t(X) %*% X) %*% t(C)
  R12=MatrixSqrtInv(R)
  A=R12 %*% C %*% Ginv(t(X) %*% X) %*% t(X)
  SCET= sum(diag(H %*% G %*% H))
  SCEC=diag(A %*% G %*% t(A))
  SCR=sum(diag((diag(I)-H) %*% G %*% (diag(I) - H)))
  Fexp= (SCET/S)/(SCR/(I-L))
  Result$Global=cbind(SCET, SCR, S, I-L, Fexp)

  FCexp=SCEC/(SCR/(I-L))

  SCEC=matrix(SCEC, ncol=1)

  Result$Contrastes=cbind(SCEC, rep(SCR, nc), rep(1, nc), rep(I-L, nc), FCexp)
  rownames(Result$Contrastes)=rownames(C)

  if (!is.null(Effects)){
    XX=t(FactorToBinary(Effects))
    tam= XX %*% t(XX)
    SCE= XX %*% SCEC
    FEexp= (solve(tam) %*% XX %*% SCEC) / (SCR/(I-L))
    rownames(FEexp)=levels(Effects)
    Result$Effects=cbind(SCE, rep(SCR, ne), diag(tam), rep(I-L, ne), FEexp)
    colnames(Result$Effects)=c("Explicada", "Residual","G.L. Num", "G.L. Denom", "F-exp")
    rownames(Result$Effects)=levels(Effects)
  }

  return(Result)
}
