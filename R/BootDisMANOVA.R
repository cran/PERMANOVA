BootDisMANOVA <- function(Distance, groups, C=NULL, Effects=NULL, nB = 1000, seed=NULL, CoordPrinc=FALSE, dimens=2, PCoA="Standard", ProjectInd=TRUE, tol=1e-4, DatosIni=TRUE) {
  D=Distance$D
  Coefficient=Distance$Coefficient
  cl <- match.call()

  PCoAs= c("Standard", "Weighted")
  if (is.numeric(PCoA)) PcoA=PCoAs(PCoA)

  if (!is.factor(groups)) stop("The grouping variable must be a factor")

  if (!is.null(seed)) set.seed(seed)
  Perm = list() #Container for the solution
  Perm$call=cl
  # Setting the properties of data
  if (is.null(rownames(D)))
    rownames(D) <- rownames(D, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(D)

  Perm$Title = "MANOVA BASED ON  based DISTANCES"
  Perm$Type = "BootMANOVA"
  Perm$Distances = D

  if (is.factor(groups)) {
    GroupNames = levels(groups)
  }

  L = length(levels(groups)) #Número de groupss
  I = dim(D)[1] #Número de de individuos
  X = FactorToBinary(groups)
  colnames(X)=levels(groups)
  rownames(X)=rownames(D)

  if (is.null(C)){
    C=diag(L)
    rownames(C)=paste("C",GroupNames)
    colnames(C)=GroupNames
    Cwasnull=TRUE
  }

  nc = dim(C)[1]

  Perm$C=C

  if (!is.null(Effects)){
    Perm$Effects=Effects
  }

  # PERMANOVA INICIAL
  Perm$Initial=PERMANOVA.Estimation(D, X, C, Effects)

  if (is.null(seed))
    set.seed(1)
  else
    set.seed(seed)

  # Permutaciones
  Ftotal=matrix(0, 1,  nB)
  FContrastes= matrix(0, nc, nB)
  if (!is.null(Effects)){
    ne=length(levels(Effects))
    FEffects=matrix(0, ne, nB)
  }

  for (i in 1:nB){
#    if (i%%100==0) print(i)
    muestra=sample.int(I, size=I, replace=TRUE)
    Man=PERMANOVA.Estimation(D[muestra, muestra], X, C, Effects)
    Ftotal[i]=Man$Global[5]
    FContrastes[,i]=Man$Contrastes[,5]
    if (!is.null(Effects)){
      FEffects[,i]=Man$Effects[,5]
    }
  }

  Perm$DistMuestral=Ftotal
  Perm$pvalue=((sum(Ftotal >= Perm$Initial$Global[1,5]) +1 )/ (nB+1))
  Perm$Initial$Global= cbind(Perm$Initial$Global, Perm$pvalue)

  colnames(Perm$Initial$Global)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")
  rownames(Perm$Initial$Global)="Total"

  Perm$DistMuestralC=FContrastes
  pvalC=(apply(FContrastes > matrix(Perm$Initial$Contrastes[,5], ncol=1) %*% matrix(1, 1, nB), 1, sum)+1)/ (nB+1)
  Perm$Initial$Contrastes=cbind(Perm$Initial$Contrastes, pvalC)
  colnames(Perm$Initial$Contrastes)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")

  if (!is.null(Effects)){
    Perm$DistMuestralEffects=FEffects
    pvalE=(apply(FEffects > matrix(Perm$Initial$Effects[,5], ncol=1) %*% matrix(1, 1, nB), 1, sum)+1)/ (nB+1)
    Perm$Initial$Effects=cbind(Perm$Initial$Effects, pvalE)
    colnames(Perm$Initial$Effects)=c("Explained", "Residual","G.L. Num", "G.L. Denom", "F-exp", "p-value")
    print(Perm$Initial$Effects)
  }

  Perm$CoordPrinc=CoordPrinc

  if (CoordPrinc){

    N=t(X) %*% X
    #D=0.5*D^2
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
    Perm$ExplainedVariance = Inertia
    Perm$Inertias= cbind(vp, Inertia, cumsum(Inertia))
    rownames(Perm$Inertias)=paste("PCo", 1:dimens)
    colnames(Perm$Inertias)=c("Eigenvalue", "Explained Variance", "Cumulative")
    Y = solut$u %*% diag(sqrt(solut$d))
    d0=apply(Y^2,1, sum)

    rownames(Y)=levels(groups)

    st <- apply(Y^2, 1, sum)
    Perm$MeanCoordinates=Y[,1:dimens]
    colnames(Perm$MeanCoordinates)=paste("Dim",1:dimens)
    colnames(Perm$MeanCoordinates)=paste("Dim",1:dimens)
    qlr <- diag(1/st) %*% (Perm$MeanCoordinates^2)
    Perm$Qualities=round(qlr[, 1:dimens]*100, digits=2)
    rownames(Perm$Qualities)=levels(groups)
    colnames(Perm$Qualities)=paste("Dim",1:dim(Perm$Qualities)[2])

    Perm$CummulativeQualities=t(apply(Perm$Qualities,1,cumsum))

    if (ProjectInd){
      Y = Y[,1:dimens]
      if (DatosIni){
        YY=as.matrix(Distance$Data)
        Means= solve(N) %*% t(X) %*% YY
        Di=DistContinuous(Means, y=YY, coef = Distance$Coefficient)$D^2
        Yi=-0.5 * Ginv(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,I,1) %*% matrix(d0,1,L))
        Yi=t(Yi)
      }
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
      Perm$RowCoordinates=Yi

    }

    Perm=AddClusterToBiplot(Perm, ClusterType="us", Groups=groups )
  }
  class(Perm)=c("BootDisMANOVA")
  return(Perm)
}



