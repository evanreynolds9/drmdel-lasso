\name{cdfDRM}
\alias{cdfDRM}
\title{Estimate the CDF of the populations under the DRM}
\description{
  Suppose we have m+1 samples, labeled as 0, 1, ..., m, whose
  population distributions satisfy the density ratio model (DRM)
  (see \code{\link{drmdel}} for the definition of DRM). The
  \code{\link{cdfDRM}} function estimates the cumulative
  distribution functions (CDF), \eqn{F_k}{F_k}, of the
  populations.
}
\usage{
cdfDRM(k, q=NULL, drmfit, interpolation=TRUE)
}
\arguments{
  \item{k}{a vector of labels of populations whose CDFs are to be
    estimated, with k[i] = 0, 1, ..., m.}
  \item{q}{It can be:

       (1) a list whose length is the same as the argument "k".
       The i\eqn{\textsuperscript{th}}{^th} component of this list
       must be a vector of quantiles at which the
       k[i]\eqn{\textsuperscript{th}}{^th} CDF is estimated.

       (2) a single vector of quantiles, in which case, for each
       population k[i], we estimate the CDF at the same quantiles.

       (3) NULL (default), in which case , for each population
       k[i], we estimate the CDF at the all the observed data
       points.
  }
  \item{drmfit}{a fitted DRM object (an output from the
    \code{\link{drmdel}} function). See \code{\link{drmdel}} for
    details.}
  \item{interpolation}{a logical variable specifying whether to
    linearly interpolate the estimated CDF to make a smooth CDF
    estimate; The default is TRUE.}
}
\value{
  A list whose length is the same as the argument "k". The
  i\eqn{\textsuperscript{th}}{^th} component of the list is a data
  frame with the following two columns:

     x -- quantiles at which the CDF of population k[i] is
     estimated.

     cdf_est -- the corresponding estimated CDF values of
     population k[i].
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

# Estimate the CDF of population 1 at c(0.5, 1.8, 3) and that of
# population 3 at c(1.0, 2.2).
cdf_est1 <- cdfDRM(k=c(1, 3), q=list(c(0.5, 1.8, 3), c(1.0, 2.2)),
                   drmfit=drmfit)
# Show the output.
names(cdf_est1)
cdf_est1$F1
cdf_est1$F3

# Estimate the CDF of population 2 and 4 at the observed data
# points.
cdf_est <- cdfDRM(k=c(2, 4), drmfit=drmfit)
# Show the output.
names(cdf_est)
cdf_est$F2
cdf_est$F4
}
\keyword{methods}
