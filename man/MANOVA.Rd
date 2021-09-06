\name{MANOVA}
\alias{MANOVA}

\title{
Multivariate Analysis of Variance (MANOVA)
}
\description{
Performs a Multivariate Analysis of Variance (MANOVA) based on matrix calculations. Is an extension of the function in the base package of R.

The function can take care of matrices of contrasts and matrices of linear combinations of variables in order to cope with complex designs.
}
\usage{
MANOVA(Y, Group, C = NULL, M = NULL, Effects = NULL, InitialTransform = 5,
AddOnes = FALSE, Contrasts = TRUE)
}

\arguments{
  \item{Y}{
A matrix containing the response variables.
}
  \item{Group}{
A factor containing the groups to compare (or the treatments).
}
  \item{C}{
A matrix of contrasts on the groups.
}
  \item{M}{
A matrix with the linear combinations of the variables to test.
}
  \item{Effects}{
A factor with the effects in the rows of C.
}
  \item{InitialTransform}{
Initial transformation of the response variables.
}
  \item{AddOnes}{
Add a column of ones to the design matrix. By default is false.
}
  \item{Contrasts}{
Should each contrast be tested separately?. By default is TRUE.
}
}
\details{
Performs a general MANOVA to compare several groups or treatments. Additional contrasts can be tested using the contrasts matrix \eqn{\bf{C}}. This can also be separated into several effects.

Here we use an approach with matrices rather than te usual multivariate general linear model. Several designs can be obtained using the appropriate matrix of contrasts.
}

\value{
An object of class "MANOVA"  and "Canonical.Biplot".
}
\references{
Amaro, I. R., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Gabriel KR (1971) The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3):453-467.

Gabriel, K. R. (1995). MANOVA biplots for two-way contingency tables. WJ Krzanowski (Ed.), Recent advances in descriptive multivariate analysis, Oxford University Press, Toronto. 227-268.

Gower y Hand (1996): Biplots. Chapman & Hall.

Varas, M. J., Vicente-Tavera, S., Molina, E., & Vicente-Villardon, J. L. (2005). Role of canonical biplot method in the study of building stones: an example from Spanish monumental heritage. Environmetrics, 16(4), 405-419.

Santana, M. A., Romay, G., Matehus, J., Vcente-Villardon, J. L., & Demey, J. R. (2009). simple and low-cost strategy for micropropagation of cassava (Manihot esculenta Crantz). African Journal of Biotechnology, 8(16).
}
\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
\donttest{
data(wine)
Y=as.matrix(wine[,4:21])
group=wine$Group
manvin=MANOVA(Y, Group=group)
summary(manvin)
}
}
