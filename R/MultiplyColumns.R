MultiplyColumns <- function(a,b){
  if (nrow(a)!=nrow(b)) stop("Matrices don't have the same number of rows")
  prod=matrix(0, nrow(a), ncol(a)*ncol(b))
  k=0
  for (i in 1:ncol(a)){
    for (j in 1:ncol(b)){
      k=k+1
      prod[,k]=a[,i]*b[,j]
    }
  }
  return(prod)
}
