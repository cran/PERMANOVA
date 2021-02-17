summary.spermanova <- function(object, Latex=FALSE, ...){
  cat(" ###### PERMANOVA : MANOVA based on distances #######\n")
  cat("________________________________________________\n")
  print(object$call)
  cat("Number of permutations:", object$nperm, "\n\n")
  SC<- c(round(object$Inicial$BSS,3), round(object$Inicial$WSS,2), round(object$Inicial$TSS,3))
  grados <- c(object$Inicial$glb, object$Inicial$glw, object$Inicial$glt)
  Fexp <- c(round(object$Inicial$Fexp,3), "", "")
  pvalor <- c(round(object$pval,10), "", "")
  res <- data.frame(SC, grados, Fexp, pvalor, row.names = c("Between Groups", "Within Groups", "TOTAL"))
  print(res)

  if (Latex){
    print(xtable(res, caption="PERMANOVA results"))
  }
}
