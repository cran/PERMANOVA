summary.BootDisMANOVA <- function(object, Latex = FALSE, ...){
  cat(" ###### Bootstrap Distance Based MANOVA Analysis #######\n\n")
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
  cat(paste("MANOVA\n"))
  print(object$Initial$Global)
  cat("________________________________________________\n\n")
  cat(paste("Contrasts\n"))
  print(rbind(object$Initial$Contrastes, object$Initial$Global))
  cat("________________________________________________\n\n")
  if (!is.null(object$Effects)){
    cat("Effects \n")
    print(rbind(object$Initial$Efectos, object$Initial$Global))
    cat("________________________________________________\n\n")
  }

  if (Latex){
    xtable(object$Initial$Global)
    xtable(rbind(object$Initial$Contrastes, object$Initial$Global))
    xtable(rbind(object$Initial$Efectos, object$Initial$Global))
  }
}
