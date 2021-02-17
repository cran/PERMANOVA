ConstructContrasts <- function(Factors, MaxOrderIter=2){
  Classes=sapply(Factors,class)
  Factors=Factors[,Classes=="factor"]
  nf=dim(Factors)[2] # Number of factors
  if (MaxOrderIter>nf) {
    print(paste("The maximum order of the iterations has been set to", nf))
    MaxOrderIter=nf}

  nlev=rep(0, nf) #Number of levels per factor
  for (i in 1:nf)
    nlev[i]=length(levels(Factors[[i]]))
  nt=prod(nlev) # Number of treatments

  # The complete set of groups
  Groups=interaction(Factors[[1]],Factors[[2]])
  if (nf>2){
    for (i in 3:nf)
      Groups=interaction(Groups,Factors[[i]])
  }

  contr=list() # This is a list with all the  contrasts for each effect

  # Contrasts for main effects
  k=1
  a1=contr.helmert(nlev[nf])
  contr[[nf]]=kronecker(a1, matrix(1, nrow=prod(nlev[1:(nf-1)]), ncol=1))

  for (i in (nf-1):1){
    k=k+1
    if (i>1){
      a2=contr.helmert(nlev[i])
      c2=kronecker(a2, matrix(1, nrow=prod(nlev[1:(i-1)]), ncol=1))
      contr[[nf-k+1]]=kronecker(matrix(1, nrow=prod(nlev[(i+1):nf]), ncol=1), c2)
    }
    else{
      a2=contr.helmert(nlev[i])
      contr[[nf-k+1]]=kronecker(matrix(1, nrow=prod(nlev[(i+1):nf]), ncol=1), a2)
    }
  }
  names(contr)=colnames(Factors)

  #  Contrasts for two-way interactions

  for (i in 1:(nf-1)){
    for (j in (i+1):nf){
      k=k+1
      contr[[k]]=MultiplyColumns(contr[[i]],contr[[j]])
      names(contr)[k]=paste(names(contr)[i], names(contr)[j], sep="*")
    }
  }

  eff=c()

  for (l in 1:k){
    eff=c(eff, rep(l, ncol(contr[[l]])))
    if (l==1) Cont=contr[[l]]
    else Cont=cbind(Cont, contr[[l]])
  }

  eff=factor(eff)
  levels(eff)=names(contr)

  return(list(Groups=Groups, Contrasts=t(Cont), Effects=eff))
}

