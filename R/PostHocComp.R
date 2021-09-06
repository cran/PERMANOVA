PostHocComp <- function(grupo){
  if (!is.factor(grupo)) stop("Debe proporcionarme un factor")
    niveles=levels(grupo)
    g=length(niveles)
    r=g*(g-1)/2
    C=matrix(0,r,g)
    colnames(C)=niveles
    filas=character(r)
    k=0
    for (i in 1:(g-1))
      for (j in (i+1):g){
        k=k+1
        C[k,i]=1
        C[k,j]=-1
        filas[k]=paste(niveles[i], niveles[j], sep="-")
      }
    rownames(C)=filas
    return(C)
}
