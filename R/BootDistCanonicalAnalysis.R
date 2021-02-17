BootDistCanonicalAnalysis <- function(Distance, groups,  dimens=NULL, nB = 100, seed=NULL,
                                      PCoA="Standard", ProcrustesRot=TRUE, DatosIni=TRUE,
                                      tol=0.0001){
  result=list()

  D=Distance$D

  result$D=D
  result$Coefficient=Distance$Coefficient
  result$nB=nB
  result$Groups=groups
  result$GroupNames = levels(groups)
  L = length(levels(groups)) #Número de grupos
  I = dim(D)[1] #Número de de individuos
  X = FactorToBinary(groups)

  colnames(X)=levels(groups)
  rownames(X)=rownames(D)

  N=t(X) %*% X
  Ns=diag(N)
  NsC=cumsum(Ns)

  if (is.null(dimens)) dimens=(L-1)


  FS=solve(N)%*% (t(X) %*% (0.5*D^2) %*% X) %*% solve(N)
  f=matrix(diag(FS), L, 1)

  DB=2*FS - f %*% matrix(1, 1, L) - t(f %*% matrix(1, 1, L))
  DB=0.5*DB

  switch(PCoA, Standard = {
    H=(diag(L) - matrix(1, L, L)/L)
    B =  H %*% DB %*% H
  },Weighted = {
    H=(diag(L) - matrix(1, L, 1) %*% matrix(diag(N), 1, L)/I)
    B = H %*% DB  %*% H
  })

  solut <- svd(B)
  b=diag(B)

  vp=solut$d[1:dimens]
  Inertia=(solut$d[1:dimens]/sum(solut$d)) * 100
  result$Inertia=Inertia
  Inertias= cbind(vp, Inertia, cumsum(Inertia))
  rownames(Inertias)=paste("PCo", 1:dimens)
  colnames(Inertias)=c("Eigenvalue", "Explained Variance", "Cumulative")
  result$Inertias=Inertias
  Y = solut$u %*% diag(sqrt(solut$d))
  Y=Y[,1:dimens]

  d0=apply(Y^2,1, sum)
  rownames(Y)=levels(groups)

  st <- apply(Y^2, 1, sum)
  result$MeanCoordinates=Y
  colnames(result$MeanCoordinates)=paste("Dim",1:dimens)
  colnames(result$MeanCoordinates)=paste("Dim",1:dimens)
  qlr <- diag(1/st) %*% (result$MeanCoordinates^2)
  result$Qualities=round(qlr[, 1:dimens]*100, digits=2)
  rownames(result$Qualities)=levels(groups)
  colnames(result$Qualities)=paste("Dim",1:dim(result$Qualities)[2])
  result$CummulativeQualities=t(apply(result$Qualities,1,cumsum))

  coord=array(0, c(L,dimens,nB))
  print("Calculating Bootstrap")
  for (i in 1:nB){
    # Obtención de las muestras bootstrap dentro de los grupos
    muestra=sample.int(Ns[1], size=Ns[1], replace=TRUE)
    for (j in 2:L){
      muestra2=sample.int(Ns[j], size=Ns[j], replace=TRUE)+NsC[j-1]
      muestra=c(muestra, muestra2)
    }

    D2=D[muestra, muestra]
    FS=solve(N)%*% (t(X) %*% (0.5*D2^2) %*% X) %*% solve(N)
    f=matrix(diag(FS), L, 1)

    DB=2*FS - f %*% matrix(1, 1, L) - t(f %*% matrix(1, 1, L))
    DB=0.5*DB

    switch(PCoA, Standard = {
      H=(diag(L) - matrix(1, L, L)/L)
      B =  H %*% DB %*% H
    },Weighted = {
      H=(diag(L) - matrix(1, L, 1) %*% matrix(diag(N), 1, L)/I)
      B = H %*% DB  %*% H
    })

    solut <- svd(B)
    vp2=solut$d[1:dimens]
    Inertia2=(solut$d[1:dimens]/sum(solut$d)) * 100
    Y2 = solut$u %*% diag(sqrt(solut$d))
    Y2=Y2[,1:dimens]
    #Y2= Y2[,1:dimens]
    if (ProcrustesRot)
      Y2=ProcrustesSimple(Y,Y2)$Yrot
    coord[,,i]=Y2[,1:dimens]
  }
  result$CoordBoot=coord


  # Coordenadas de los individuos
  if (DatosIni){
    YY=as.matrix(Distance$Data)
    Means= solve(N) %*% t(X) %*% YY
    Di=DistContinuous(Means, y=YY, coef = Distance$Coefficient)$D^2
    Yi=-0.5 * ginv(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,I,1) %*% matrix(d0,1,L))
    Yi=t(Yi)}
  else{
    G <- (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I) %*% (-0.5 * D^2) %*% (diag(I) - matrix(1, I, 1) %*% matrix(1, 1, I) / I)
    soluc <- svd(G)
    dimefec=sum(soluc$d > tol)
    YY=soluc$u[,1:dimefec] %*% diag(sqrt(soluc$d[1:dimefec]))
    Means= solve(N) %*% t(X) %*% YY
    Di=DistContinuous(Means, y=YY, coef = 1)$D
    Yi=-0.5 * solve(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,I,1) %*% matrix(d0,1,L))
    Yi=t(Yi)
  }

  rownames(Yi)=rownames(D)
  result$RowCoordinates=Yi

  class(result)="BootCanonAnalysis"
  return(result)
}



