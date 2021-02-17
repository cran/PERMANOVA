print.BootDisMANOVA <- function(x, ...){
  cat(" ###### PERMANOVA Analysis #######\n\n")
  cat("Call\n")
  print(x$call)
  cat("________________________________________________\n\n")
  cat("Contrast Matrix\n")
  print(x$C)
  cat("________________________________________________\n\n")
  if (!is.null(x$Effects)){
    cat("Effects\n")
    print(x$Effects)
    cat("________________________________________________\n\n")
  }

  cat(paste("MANOVA\n"))
  print(x$Initial$Global)

}
