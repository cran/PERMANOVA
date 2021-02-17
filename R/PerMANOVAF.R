PerMANOVAF <- function(D, grupo, C, Efectos){
  g = length(levels(grupo))
  n = dim(D)[1]
  G = FactorToBinary(grupo) # Matrix of indicators
  Eij = G %*% t(G)
  TSS=0.5*sum(D^2)/n
  ng=diag(t(G) %*% G)
  WSS=sum(0.5*apply((Eij * D^2)%*%G, 2, sum)/ng)
  BSS=TSS-WSS
  Fexp=(BSS/(g-1))/(WSS/(n-g))
  Result=list()
  Result$TSS= TSS
  Result$BSS = BSS
  Result$WSS = WSS
  Result$glt=n-1
  Result$glb=g-1
  Result$glw=n-g
  Result$Fexp=Fexp
  return(Result)
}
