summary.BootDisMANOVA <- function(object, Latex = FALSE, ...) {
  if (!Latex) {
    cat(" ###### Bootstrap Distance Based MANOVA Analysis #######\n\n")
    cat("Call\n")
    print(object$call)
    cat("________________________________________________\n\n")
    cat("Contrast Matrix\n")
    print(object$C)
    cat("________________________________________________\n\n")
    if (!is.null(object$Effects)) {
      cat("Effects vector\n")
      print(object$Effects)
      cat("________________________________________________\n\n")
    }
    cat("MANOVA\n")
    print(object$Initial$Global)
    cat("________________________________________________\n\n")
    cat("Contrasts\n")
    print(rbind(object$Initial$Contrastes, object$Initial$Global))
    cat("________________________________________________\n\n")
    if (!is.null(object$Effects)) {
      cat("Effects\n")
      print(rbind(object$Initial$Effects, object$Initial$Global))
      cat("________________________________________________\n\n")
    }
  } else {
    print(xtable::xtable(object$Initial$Global, caption = "Bootstrap Distance MANOVA"))
    print(xtable::xtable(rbind(object$Initial$Contrastes, object$Initial$Global), caption = "Contrasts"))
    if (!is.null(object$Effects)) {
      print(xtable::xtable(rbind(object$Initial$Effects, object$Initial$Global), caption = "Effects"))
    }
  }
}

