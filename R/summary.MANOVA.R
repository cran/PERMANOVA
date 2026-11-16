summary.MANOVA <- function(object, ShowCanonical = FALSE, latex = FALSE, ...) {
  row_cumsum <- function(mat) {
    if (ncol(mat) == 1) return(mat)
    res <- t(apply(mat, 1, cumsum))
    colnames(res) <- colnames(mat)
    rownames(res) <- rownames(mat)
    res
  }

  if (!latex) {
    cat(" ###### MULTIVARIATE ANALYSIS OF VARIANCE : MANOVA  #######\n\n")
    cat("________________________________________________\n\n")
    cat("Contrasts matrix\n\n")
    print(object$C)
    cat("________________________________________________\n\n")
    cat("Global Test\n\n")
    print(object$MANOVA)
    cat("________________________________________________\n\n")

    if (!is.null(object$ANOVAS)) {
      cat("Analysis of Variance for each Variable\n")
      print(object$ANOVAS)
      cat("________________________________________________\n\n")
    }

    if (!is.null(object$Effects)) {
      cat("Test for each Effect\n")
      print(object$Effects)
      cat("________________________________________________\n\n")
    }

    if (!is.null(object$Contrasts)) {
      cat("Test for each separate contrast\n")
      print(object$Contrasts)
      cat("________________________________________________\n\n")
    }

    cat("Eigenvalues and explained variance\n")
    pp <- cbind(seq_along(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
    colnames(pp) <- c("Axis", "Eigenvalue", "Explained Variance", "Cumulative")
    print(pp)
    cat("________________________________________________\n\n")

    if (ShowCanonical) {
      cat("Correlations between the canonical axis and the original variables\n")
      print(round(object$Structure_Correlations, digits = 3))
      cat("________________________________________________\n\n")
      cat("Squared Correlations\n")
      print(round(object$Structure_Correlations^2, digits = 3))
      cat("________________________________________________\n\n")
      cat("Quality of representation of the group means\n")
      print(round(object$GroupContributions * 100, digits = 3))
      cat("________________________________________________\n\n")
      cat("Quality of representation of the group means - Cumulative\n")
      print(round(row_cumsum(object$GroupContributions) * 100, digits = 3))
      cat("________________________________________________\n\n")
      cat("Goodness of fit of the variables (to explain the group means)\n")
      print(round(object$ColContributions * 100, digits = 3))
      cat("________________________________________________\n\n")
      cat("Goodness of fit of the variables (to explain the group means) - Cumulative\n")
      print(round(row_cumsum(object$ColContributions) * 100, digits = 3))
      cat("________________________________________________\n\n")
    }
  } else {
    print(xtable::xtable(object$C, digits = 0, caption = "Contrast Matrix"))
    print(xtable::xtable(object$MANOVA, digits = c(0, 3, 3, 0, 0, -3), caption = "MANOVA"))

    if (!is.null(object$ANOVAS)) {
      print(xtable::xtable(object$ANOVAS, caption = "Analysis of Variance for each Variable"))
    }

    if (!is.null(object$Effects)) {
      eff_list <- object$Effects
      attr(eff_list, "subheadings") <- names(eff_list)
      xList <- xtable::xtableList(eff_list, digits = c(0, 3, 3, 0, 0, -3), caption = "Test for each Effect")
      print(xList, colnames.format = "multiple")
    }

    if (!is.null(object$Contrasts)) {
      cont_list <- object$Contrasts
      attr(cont_list, "subheadings") <- names(cont_list)
      xList <- xtable::xtableList(cont_list, digits = c(0, 3, 3, 0, 0, -3), caption = "Test for each separate contrast")
      print(xList, colnames.format = "multiple")
    }

    pp <- cbind(seq_along(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
    colnames(pp) <- c("Axis", "Eigenvalue", "Explained Variance", "Cumulative")
    print(xtable::xtable(pp, caption = "Eigenvalues and explained variance"))

    if (ShowCanonical) {
      print(xtable::xtable(round(object$Structure_Correlations, digits = 3),
                           caption = "Correlations between the canonical axis and the original variables"))
      print(xtable::xtable(round(object$Structure_Correlations^2, digits = 3),
                           caption = "Squared Correlations"))
      print(xtable::xtable(round(object$GroupContributions * 100, digits = 3),
                           caption = "Quality of representation of the group means"))
      print(xtable::xtable(round(row_cumsum(object$GroupContributions) * 100, digits = 3),
                           caption = "Quality of representation of the group means - Cumulative"))
      print(xtable::xtable(round(object$ColContributions * 100, digits = 3),
                           caption = "Goodness of fit of the variables (to explain the group means)"))
      print(xtable::xtable(round(row_cumsum(object$ColContributions) * 100, digits = 3),
                           caption = "Goodness of fit of the variables (to explain the group means) - Cumulative"))
    }
  }
}

