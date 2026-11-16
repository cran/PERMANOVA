# PERMANOVA 0.2.1

## New features

* Added unified distance wrapper `Distances()` with automatic variable type detection.
* Added `DistCategorical()` (Cramer's V and matching coefficient) and `DistMixed()` (Gower, generalized Gower, and robust RELMS distances) powered by `dbrobust`.
* Expanded `DistContinuous()` with Mahalanobis and robust covariance-based distances (`Robust_GV`, `Robust_MCD`), as well as optional Euclidean correction (`euclidean = TRUE`) across distance functions.
* Added parallel permutation support (`parallel` and `cores` arguments) in `PERMANOVA()`, compliant with CRAN core limits.

## Bug fixes and improvements

* **Graphical methods (`plot.BootDisMANOVA`, `plot.PERMANOVA`, `plot.BootCanonAnalysis`, `plot.MANOVA.Biplot`, `plot.MANOVA`)**: Fixed plot boundaries so projected individuals and confidence regions are no longer clipped; fixed cluster/group factor references (`x$Clusters`) and axis indexing when `A1`/`A2` differ from `1` and `2`; removed non-functional redundant arguments.
* **Principal coordinates & bootstrap**: Corrected sign and `eigen()` decomposition in Gower double-centering within `PERMANOVA()`, `BootDisMANOVA()`, and `BootDistCanonicalAnalysis()`, fixed individual projection scaling when `DatosIni = FALSE`, and made stratified bootstrap resampling independent of row ordering.
* **MANOVA & contrasts**: Fixed Wilks' lambda calculation and dimension dropping in two-group designs (`r = 1`), corrected degrees of freedom for Roy's greatest root in `MANOVA.Estimation()`, enabled `MaxOrderIter` in `ConstructContrasts()`, and fixed `Logit` and normalized residuals in `IniTransform()`.
* **Summary & LaTeX output**: Fixed field references (`Effects`) and row-wise cumulative contributions in `summary.MANOVA()`, `summary.BootDisMANOVA()`, and `summary.PERMANOVA()`.
