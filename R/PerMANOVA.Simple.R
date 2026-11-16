PerMANOVA.Simple <- function(D, grupo, nperm = 999, seed = NULL, C = NULL) {
  cl <- match.call()
  if (!is.factor(grupo)) stop("The grouping variable must be a factor")
  if (!is.null(seed)) set.seed(seed) else set.seed(1)

  Result <- list()
  Result$call <- cl
  Result$nperm <- nperm
  Result$Inicial <- PerMANOVAF(D, grupo)
  n <- length(grupo)
  Result$Fvals <- rep(0, nperm)

  for (i in seq_len(nperm)) {
    muestra <- sample.int(n)
    Result$Fvals[i] <- PerMANOVAF(D[muestra, muestra], grupo)$Fexp
  }
  Result$pval <- (sum(Result$Fvals >= Result$Inicial$Fexp) + 1) / (nperm + 1)
  class(Result) <- "permanova"
  return(Result)
}

