\name{densityDRM}
\alias{densityDRM}
\title{Estimate the density of the populations under DRM}
\description{
  Suppose we have m+1 samples, labeled as
  \eqn{0, \, 1, \, \ldots, \, m}{0, 1, ..., m},
  whose population distributions satisfy the density ratio model
  (DRM) (see \code{\link{drmdel}} for the definition of DRM). The
  \code{\link{densityDRM}} function estimates the density,
  \eqn{dF_k(x)}{dF_k(x)}, of the k\eqn{\textsuperscript{th}}{^th},
  \eqn{k=0, \, 1, \, \ldots, \, m}{k = 0, 1, ..., m},
  population using a kernel density estimator with weights been
  the estimated \eqn{dF_k(x)}{dF_k(x)}'s at the observed data
  points.
}
\usage{
densityDRM(k, drmfit, interpolation=TRUE, adjust=FALSE,
           adj_factor=NULL, ...)
}
\arguments{
  \item{k}{the label of the population whose density is to be
    estimated, k = 0, 1, ..., m. It must be a single value, not a
    vector.}
  \item{drmfit}{a fitted DRM object (an output from the
    \code{\link{drmdel}} function). See \code{\link{drmdel}} for
    details.}
  \item{interpolation}{an argument passed to the function
    \code{\link{quantileDRM}} for calculating the default
    bandwidth for the kernel density estimator, if bandwidth (bw)
    is not specfied as a ... argument passing to the R functin
    \code{density}. The default is TRUE. See
    \code{\link{quantileDRM}} for details.}
  \item{adjust}{an argument passed to the function
    \code{\link{quantileDRM}} for calculating the default
    bandwidth for the kernel density estimator, if bandwidth (bw)
    is not specfied as a ... argument passing to R functin
    \code{density}. The default is FALSE. See
    \code{\link{quantileDRM}} for details.}
  \item{adj_factor}{an argument passed to the function
    \code{\link{quantileDRM}} for calculating the default
    bandwidth for the kernel density estimator, if bandwidth (bw)
    is not specfied as a ... argument passing to R functin
    \code{density}. The default is NULL. See
    \code{\link{quantileDRM}} for details.}
  \item{...}{further arguments to be passed to the \code{density}
    function, which performs kernel density estimation in R. See
    \code{help(density)} for details. One should not pass "x" and
    "weights" arguments to the \code{density} function, since
    those are supposed to be extracted automatically from the
    fitted DRM object, drmfit.}
}
\value{
  An output from the \code{density} function, usually an object
  with class '"density"'. See \code{help(density)} for details.
}
\details{
  Note that the default bandwidth for this density estimator is
  set as that described in Chen and Liu (2013):
  \deqn{1.06 n_{total}^{-1/5} \min(\sigma_k, \, \mbox{IQR}_k/1.34)}{
      1.06 * n_total^(-1/5) * min(sigma_k, IQR_k/1.34)},
  where \eqn{n_{total}}{n_total} is the total sample size, and
  \eqn{sigma_k}{sigma_k} and \eqn{IQR_k}{IQR_k} are the standard
  deviation and inter-quartile range of the estimated CDF
  \eqn{F_k}{F_k}, respectively.

  If bandwidth (bw) is not specfied as a ... argument passing to
  the R functin \code{density}, the default bandwidth (as
  described above) will be calculated. That calculation involves
  the estimation of population quartiles. In this situation, the
  arguments "interpolation", "adjust", and "adj_factor" are passed
  to \code{\link{quantileDRM}} for quartile estimations. See
  \code{\link{quantileDRM}} for details.
}
\references{
  K. Fokianos (2004),
  Merging information for semiparametric density estimation.
  \emph{Journal of the Royal Statistical Society, Series B
    (Statistical Methodology)}, 66(4):941-958.

  J. Chen and Y. Liu (2013),
  Quantile and quantile-function estimations under density ratio
  model.
  To appear in \emph{The Annals of Statistics}, 2013.
}
\examples{
# Data generation
set.seed(25)
n_samples <- c(100, 200, 180, 150, 175)  # sample sizes
x0 <- rgamma(n_samples[1], shape=5, rate=1.8)
x1 <- rgamma(n_samples[2], shape=12, rate=1.2)
x2 <- rgamma(n_samples[3], shape=12, rate=1.2)
x3 <- rgamma(n_samples[4], shape=18, rate=5)
x4 <- rgamma(n_samples[5], shape=25, rate=2.6)
x <- c(x0, x1, x2, x3, x4)

# Fit a DRM with the basis function q(x) = (x, log(abs(x))), which
# is the basis function for gamma family. This basis function is
# the built-in basis function 6.
drmfit <- drmdel(x=x, n_samples=n_samples, basis_func=6)

# Estimate the density of population 3 under the DRM
dens_pop3 <- densityDRM(k=3, drmfit=drmfit)

# Plot the estimated density
plot(dens_pop3, main=bquote(F[3]), ylim=range(c(0, 0.5)))

# Add the empirical kernel density estimation curve of F_3 based
# on the third sample on the above density plot
lines(density(x3), col="blue", lty="28F8")

# Add the true density curve of F_3 on the above density plot
lines(seq(min(dens_pop3$y), max(dens_pop3$x), 0.01),
      dgamma(seq(min(dens_pop3$y), max(dens_pop3$x), 0.01), 18,
             5),
      type="l", col="red", lty="dotted")

legend(9, 0.5,
       legend=c("DRM density estimator",
                "Empirical kernel density estimator",
                "True density"),
       col=c("black", "blue", "red"),
       lty=c("solid", "28F8", "dotted"))
}
\keyword{methods}
