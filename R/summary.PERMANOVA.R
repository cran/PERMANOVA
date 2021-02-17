summary.PERMANOVA <- function(object, Latex=FALSE, ...){
  cat(" ###### PERMANOVA Analysis #######\n\n")
  cat("Call\n")
  print(object$call)
  cat("________________________________________________\n\n")
  cat("Contrast Matrix\n")
  print(object$C)
  cat("________________________________________________\n\n")
  if (!is.null(object$Effects)){
    cat("Effects vector\n")
    print(object$Effects)
    cat("________________________________________________\n\n")
  }
  cat(paste("PerMANOVA\n"))
  print(object$Initial$Global)
  cat("________________________________________________\n\n")
  cat(paste("Contrasts\n"))
  print(rbind(object$Initial$Contrastes, object$Initial$Global))
  cat("________________________________________________\n\n")
  if (!is.null(object$Effects)){
    cat("Effects \n")
    print(rbind(object$Initial$Effects, object$Initial$Global))
    cat("________________________________________________\n\n")
  }
  if (Latex){
    xtable(round(object$C, digits=0), caption="Contrast Matrix")
    cat("________________________________________________\n\n")
    xtable(object$Initial$Global, caption="PerMANOVA")
    cat("________________________________________________\n\n")
    xtable(rbind(object$Initial$Contrastes, object$Initial$Global), caption="perMANOVA with contrasts")

  }
}

