summary.MANOVA <- function(object, ShowCanonical=FALSE, latex=FALSE, ...) {
  if (!latex){
    cat(" ###### MULTIVARIATE ANALYSIS OF VARIANCE : MANOVA  #######\n\n")
    cat("________________________________________________\n\n")
    cat("Contrasts matrix\n\n")
    print(object$C)
    cat("________________________________________________\n\n")
    cat("Global Test\n\n")
    print(object$MANOVA)
    cat("________________________________________________\n\n")

    if (!is.null(object$ANOVAS)){
      cat("Analysis of Variance for each Variable\n")
      print(object$ANOVAS)
      cat("________________________________________________\n\n")
    }

    if (!is.null(object$Effects)){
      cat("Test for each Effect\n")
      print(object$Effects)
      cat("________________________________________________\n\n")
    }

    if (!is.null(object$Contrasts)){
      cat("Test for each separate contrast\n")
      print(object$Contrasts)
      cat("________________________________________________\n\n")
    }


    cat("Eigenvalues and explained variance\n")
    pp = cbind(1:length(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
    colnames(pp) = c("Axis", "Eigenvalue", "Explained Variance", "Cummulative")
    print(pp)
    cat("________________________________________________\n\n")

    if (ShowCanonical){
      cat("Correlations between the canonical axis and the original variables\n")
      print(round(object$Structure_Correlations, digits=3))
      cat("________________________________________________\n\n")
      cat("Squared Correlations \n")
      print(round(object$Structure_Correlations^2, digits = 3))
      cat("________________________________________________\n\n")
      cat("Quality of representation of the group means \n")
      print(round(object$GroupContributions*100, digits=3))
      cat("________________________________________________\n\n")
      cat("Quality of representation of the group means - Cummulative\n")
      print(round(cumsum(object$GroupContributions)*100, digits=3))
      cat("________________________________________________\n\n")
      cat("Goodness of fit if the variables (to explain the group means) \n")
      print(round(object$ColContributions*100, digits=3))
      cat("________________________________________________\n\n")
      cat("Goodness of fit if the variables (to explain the group means) - Cummulative \n")
      print(round(cumsum(object$ColContributions*100), digits=3))
      cat("________________________________________________\n\n")}
  }


  if (latex){
    print(xtable(object$C, digits=0, caption="Contrast Matrix"))

    print(xtable(object$MANOVA,  digits=c(0, 3, 3, 0, 0, -3), caption="MANOVA"))

    if (!is.null(object$ANOVAS)){
      print(xtable(object$ANOVAS, caption="Analysis of Variance for each Variable"))
    }

    if (!is.null(object$Effects)){
      attr(object$Efectos, "subheadings") <- names(object$Effects)
      xList=xtableList(object$Effects,  digits=c(0, 3, 3, 0, 0, -3), caption="Test for each Effect")
      print(xList, colnames.format = "multiple")
    }

    if (!is.null(object$Contrasts)){
      attr(object$Contrasts, "subheadings") <- names(object$Contrasts)
      xList=xtableList(object$Contrasts,  digits=c(0, 3, 3, 0, 0, -3), caption="Test for each separate contrast")
      print(xList, colnames.format = "multiple")
    }




    pp = cbind(1:length(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
    colnames(pp) = c("Axis", "Eigenvalue", "Explained Variance", "Cummulative")
    print(xtable(pp, caption="Eigenvalues and explained variance"))


    if (ShowCanonical){
      cat("Correlations between the canonical axis and the original variables\n")
      print(xtable(round(object$Structure_Correlations, digits=3), caption="Correlations between the canonical axis and the original variables"))

      print(xtable(round(object$Structure_Correlations^2, digits=3), caption="Squared Correlations"))

      print(xtable(round(object$GroupContributions*100, digits=3), caption="Quality of representation of the group means"))

      print(xtable(round(cumsum(object$GroupContributions)*100, digits=3), caption="Quality of representation of the group means - Cummulative"))

      print(xtable(round(object$ColContributions*100, digits=3), caption="Goodness of fit if the variables (to explain the group means)"))
      cat("________________________________________________\n\n")
      cat("Goodness of fit of the variables (to explain the group means) - Cummulative \n")
      print(round(cumsum(object$ColContributions*100), digits=3))
      print(xtable(round(cumsum(object$ColContributions*100), digits=3), caption="Goodness of fit if the variables (to explain the group means) - Cummulative"))
      cat("________________________________________________\n\n")}

  }
}
