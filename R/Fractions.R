Fractions <- function(data, center=NULL, confidence=1){
  if (is.null(center))
    center=apply(data,2,mean)

  n=dim(data)[1]
  datacen=(diag(n)-matrix(1, n, n)/n) %*% data
  DE=svd(datacen)
  K=DE$u
  d=matrix(0,n,1)
  for (i in 1:n)
    d[i]=sqrt(sum(K[i,]^2))
  radius = quantile(d,confidence)

  data2=datacen[which(d<=radius),]
  data2[,1]=data2[,1]+center[1]
  data2[,2]=data2[,2]+center[2]
  res=list(data=data, fraction=data2, confidence=confidence, center=center)
  class(res)<-"fraction"
  return(res)
}
