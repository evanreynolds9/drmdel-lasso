\name{quantileCompWald}
\alias{quantileCompWald}
\title{Wald test for pairwise comparison and linear
  hypotheses about quantiles under the DRM}
\description{
  Suppose we have m+1 samples, labeled as
  \eqn{0, \, 1, \, \ldots, \, m}{0, 1, ..., m},
  whose population distributions satisfy the density ratio
  model (DRM) (see \code{\link{drmdel}} for the definition
  of DRM). We now want to test the linear hypothesis about a
  vector of quantiles
  \eqn{q = (q_1, q_2, \ldots, q_s)^T}{q = (q_1, q_2,
    ..., q_s)^T}
  of probably different populations:
  \deqn{H_0: \, Aq = b \ \  \mbox{against} \  \ H_1: \,
    Aq\neq b,}{H_0: Aq = b  against H_1:  Aq != b,}
  where \eqn{A}{A} is a \eqn{t \times s}{t by s},
  \eqn{t \le s}{t <= s}, non-singular matrix and b is a
  vector. The \code{\link{quantileCompWald}} function
  performs a Wald-test for the above hypothesis and also
  pairwise comparisons of the population quantiles.
}
\usage{
quantileCompWald(quantileDRMObject, n_total, pairwise=TRUE,
                 p_adj_method="none", A=NULL, b=NULL)
}
\arguments{
  \item{quantileDRMObject}{an output from the
    \code{\link{quantileDRM}} function. It must be a list
    containing a vector of quantile estimates
    (quantileDRMobject$est) and a estimated covariance
    matrix of the quantile estimates
    (quantileDRMobject$cov). That is, the argument 'cov'
    must be set to 'TRUE' when running quantileDRM.}
  \item{n_total}{total sample size.}
  \item{pairwise}{a logical variable specifying whether to
    perform pairwise comparisons of the quantiles. The
    default is TRUE.}
  \item{p_adj_method}{when pairwise=TRUE, how should the
    p-values be adjusted for multiple comparisons. The
    available methods are: "holm", "hochberg", "hommel",
    "bonferroni", "BH", "BY", "fdr" and "none". See
    \code{help(p.adjust)} for details. The default is
    "none", i.e., no adjustment.}
  \item{A}{the left-hand side t by s, t \eqn{\le}{<=} s,
    matrix in the linear hypothesis.}
  \item{b}{the right-hand side t-dimensional vector in the
    linear hypothesis.}
}
\value{
  \item{p_val_pair}{p-values of pairwise comparisons, in the
    form of a lower triangular matrix. Available only if
    argument pairwise=TRUE}
  \item{p_val}{p-value of the linear hypothesis. Available
    only if argumen 'A' and 'b' are not NULL.}
}
\details{
  Denote the EL quantile estimate of the \eqn{q}{q} vector
  as \eqn{\hat{q}}{qHat}, and the estimate of the
  corresponding covariance matrix as
  \eqn{\hat{\Sigma}}{SigmaHat}. \eqn{\hat{q}}{qHat} and
  \eqn{\hat{\Sigma}}{SigmaHat} can be calculated using
  function \code{\link{quantileDRM}} with 'cov=TRUE'.
  
  It is known that, \eqn{\sqrt{n_{total}}(\hat{q} -
    q)}{sqrt(n_total) * (qHat - q)} converges in
  distribution to a normal distribution with 0 mean and
  covariance matrix \eqn{\Sigma}{Sigma}. Also,
  \eqn{\hat{\Sigma}}{SigmaHat} is a consistent estimator of
  \eqn{\Sigma}{Sigma}. Hence, under the null of the linear
  hypothesis,
  \deqn{n_total (A \hat{q} - b)^T {(A \hat\Sigma A^T)}^{-1}
    (A \hat{q} - b)}{n_total * (A * qHat - b)^T (A *
    \hat\Sigma * A^T)^{-1} (A * qHat - b)}
  has a chi-square limiting distribution with t (=ncol(A))
  degrees of freedom.
}
\references{
  J. Chen and Y. Liu (2013),
  Quantile and quantile-function estimations under density
  ratio model.
  \emph{The Annals of Statistics}, 41(3):1669-1692.
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

# Fit a DRM with the basis function q(x) = (x, log(abs(x))),
# which is the basis function for gamma family. This basis
# function is the built-in basis function 6.
drmfit <- drmdel(x=x, n_samples=n_samples, basis_func=6)

# Quantile comparisons
# Compare the 5^th percentile of population 0, 1, 2 and 3.

# Estimate these quantiles first
qe <- quantileDRM(k=c(0, 1, 2, 3), p=0.05, drmfit=drmfit)

# Create a matrix A and a vector b for testing the equality
# of all these 5^th percentiles. Note that, for this test,
# the contrast matrix A is not unique.
A <- matrix(rep(0, 12), 3, 4)
A[1,] <- c(1, -1, 0, 0)
A[2,] <- c(0, 1, -1, 0)
A[3,] <- c(0, 0, 1, -1)
b <- rep(0, 3)

# Quantile comparisons
# No p-value adjustment for pairwise comparisons
(qComp <- quantileCompWald(qe, n_total=sum(n_samples), A=A,
                           b=b))

# Adjust the p-values for pairwise comparisons using the
# "holm" method.
(qComp1 <- quantileCompWald(qe, n_total=sum(n_samples),
                            p_adj_method="holm", A=A, b=b))
}
\keyword{methods}
