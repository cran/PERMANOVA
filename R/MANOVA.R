MANOVA <- function(Y, Group, C=NULL, M=NULL, Effects=NULL, InitialTransform = 5, AddOnes=FALSE, Contrasts=TRUE) {

  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", "Column centering", "Standardize columns", "Row centering",
                              "Standardize rows", "Divide by the column means and center", "Normalized residuals from independence")
  if (is.numeric(InitialTransform))
    InitialTransform = ContinuousDataTransform[InitialTransform]
  Bip = list() #Container for the solution
  # Setting the properties of data
  if (is.null(rownames(Y)))
    rownames(Y) <- rownames(Y, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(Y)
  if (is.null(colnames(Y)))
    colnames(Y) <- colnames(X, do.NULL = FALSE, prefix = "V")


  Bip$Title = "MANOVA Biplot"
  Bip$Type = "MANOVA"
  Bip$Non_Scaled_Data = Y
  Bip$Means = apply(Y, 2, mean)
  Bip$Medians = apply(Y, 2, median)
  Bip$Deviations = apply(Y, 2, sd)
  Bip$Minima = apply(Y, 2, min)
  Bip$Maxima = apply(Y, 2, max)
  Bip$P25 = apply(Y, 2, quantile)[2, ]
  Bip$P75 = apply(Y, 2, quantile)[4, ]
  Bip$GMean = mean(as.matrix(Y))
  Bip$Initial_Transformation=InitialTransform
  Y = IniTransform(as.matrix(Y), transform = InitialTransform) # Initial transformation
  rownames(Y) <- RowNames
  if (is.factor(Group)) {
    GroupNames = levels(Group)
  }
  L = length(levels(Group)) #Número de grupos
  I = dim(Y)[1] #Número de de individuos
  J = dim(Y)[2] #Número de variables
  r = min(c(L - 1, J)) # Rango de la solucion
  X = FactorToBinary(Group)
  colnames(X)=levels(Group)

  if (AddOnes){
    X=cbind(rep(1,I),X)
    rownames(X)=RowNames
    colnames(X)=c("Constante",levels(Group) )
  }

  if (is.null(C)){
    C=diag(L)
    rownames(C)=paste("C",GroupNames)
    colnames(C)=GroupNames
    Cwasnull=TRUE
  }
  if (is.null(M)) {
    M=diag(J)
    colnames(M)=colnames(Y)
  }
  VarNames = colnames(M)

  Bip$ncols=J
  Bip$nrows=I

  Bip$g=L
  Bip$p=J

  Bip$C=C
  Bip$Combinaciones=M
  DimNames = paste("Dim", 1:r)

  # Manova Completo
  Man=MANOVA.Estimation(Y, X, C, M)

  Bip$MANOVA=Man$MANOVA
  Bip$EigenValues = Man$VP
  Bip$Inertia = Man$Inercia
  Bip$CumInertia = cumsum(Man$Inercia)

  Bip$GroupCoordinates = Man$B %*% M %*% Man$A[,1:r]
  rownames(Bip$GroupCoordinates) = GroupNames
  colnames(Bip$GroupCoordinates) = paste("Dim", 1:r)


  Bip$RowCoordinates = Y %*% M %*% Man$A[,1:r]
  rownames(Bip$RowCoordinates) = RowNames
  colnames(Bip$RowCoordinates) = paste("Dim", 1:r)

  Bip$ColCoordinates = MatrixSqrt(Man$E) %*% Man$V[,1:r]
  rownames(Bip$ColCoordinates) = VarNames
  colnames(Bip$ColCoordinates) = paste("Dim", 1:r)

  Bip$Structure_Correlations = cor(Y %*% M, Bip$RowCoordinates)
  rownames(Bip$Structure_Correlations) = VarNames
  colnames(Bip$Structure_Correlations) = paste("Dim", 1:r)

  Bip$Canonical_Weights = Man$A[,1:r]
  rownames(Bip$Canonical_Weights) = VarNames
  colnames(Bip$Canonical_Weights) = paste("Dim", 1:r)

  Bip$GroupContributions = diag(1/rowSums(Bip$GroupCoordinates^2)) %*% Bip$GroupCoordinates^2
  Bip$ColContributions = diag(1/rowSums(Bip$ColCoordinates^2)) %*% Bip$ColCoordinates^2

  Bip$ExplTotal = matrix(0, r, 1)
  Bip$RowContributions = matrix(0, I, r)
  Bip$QLRVars = matrix(0, J, r)

  SCT = sum(Y^2)
  SCRows = rowSums(Y^2)
  SCCols = colSums(Y^2)

  for (j in 1:r) {
    Fitted = Bip$RowCoordinates[, 1:j] %*% t(Bip$ColCoordinates[, 1:j])
    residuals = Y %*%M - Fitted
    Bip$ExplTotal[j] = 1 - sum(residuals^2)/SCT
    Bip$RowContributions[, j] = 1 - apply(residuals^2,1,sum)/SCRows
    Bip$QLRVars[, j] = 1 - colSums(residuals^2)/SCCols
  }


  Ngrupos=length(levels(Group))
  Bip$Clusters = Group
  Bip$ClusterNames = levels(Group)

  palette(rainbow(Ngrupos))
  ClusterColors = palette()
  Bip$ClusterType="us"
  Bip$ClusterColors=ClusterColors

  falfau = qt(1 - (0.025), (I - L))
  falfab = qt(1 - (0.025/(I - L)), (I - L))
  falfam = sqrt(qf(1 - 0.05, L, (I - L - J + 1)) * (((I - L) * J)/(I - L - J + 1)))
  falfac = sqrt(qchisq(0.95, 2))
  IL=table(Group)
  Bip$UnivRad = falfau * (1/sqrt(IL))/sqrt(I - L)
  Bip$BonfRad = falfab * (1/sqrt(IL))/sqrt(I - L)
  Bip$MultRad = falfam * (1/sqrt(IL))/sqrt(I - L)
  Bip$ChisRad = falfac * (1/sqrt(IL))/sqrt(I - L)

  # Todo lo anterior corresponde con el modelo completo
  # si no hay matriz de contrastes sería el biplot canónico habitual
  # Efectos es un factor con tantas observaciones como contrastes
  # El factor contiene los contrastes que se corresponden con el mismo efecto
  # Ahora separaremos cada uno de los efectos

  nContrasts=dim(C)[1]


  if (!is.null(Effects)){
    EstimEfectos=list()
    S=dim(C)[2]
    efectos=levels(Effects)
    nefec=length(efectos)
    for (i in 1:nefec){
      cuales=which(Effects==efectos[i])
      Man=MANOVA.Estimation(Y, X, C[cuales,], M)
      EstimEfectos[[i]]=Man$MANOVA
    }
    names(EstimEfectos)=efectos
    Bip$Effects=EstimEfectos
  }

  if (Contrasts){
    EstimContrastes=list()
    S=dim(C)[2]
    for (i in 1:nContrasts){
      Man=MANOVA.Estimation(Y, X, C[i,], M)
      EstimContrastes[[i]]=Man$MANOVA
    }
    names(EstimContrastes)=rownames(C)
    Bip$Contrasts=EstimContrastes
  }

  class(Bip)=c("MANOVA", "MANOVA.Biplot")
  return(Bip)
}

