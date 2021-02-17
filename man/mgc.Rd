\name{mgc}
\alias{mgc}

\title{
Mixture Gaussian Clustering
}
\description{
Model based clustering using mixtures of gaussian distributions.
}
\usage{
mgc(x, NG = 2, init = "km", RemoveOutliers = FALSE, ConfidOutliers = 0.995,
tolerance = 1e-07, maxiter = 100, show = TRUE, ...)
}

\arguments{
  \item{x}{
The data matrix.
}
  \item{NG}{
Number of groups or clusters to obtain.
}
  \item{init}{
Initial centers can be obtained from k-means ("km") or at random ("rd").
}
  \item{RemoveOutliers}{
Should the extreme values be removed to calculate the clusters?
}
  \item{ConfidOutliers}{
Percentage of the points to keep for the calculations when RemoveOutliers is true.
}
  \item{tolerance}{
Tolerance for convergence.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{show}{
Should the likelihood at each iteration be shown?
}
  \item{\dots}{
Any other parameter that can affect k-means if that is the initial configuration.
}
}
\details{
A basic algorithm for clustering with mixtures of gaussians with no restrictions on the covariance matrices.
}
\value{
Clusters.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
X=as.matrix(iris[,1:4])
mod1=mgc(X,NG=3)
plot(iris[,1:4], col=mod1$Classification)
table(iris[,5],mod1$Classification)
}
