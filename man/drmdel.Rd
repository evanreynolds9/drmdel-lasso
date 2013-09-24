\name{drmdel}
\alias{drmdel}
\title{Fit a density ratio model}
\description{
  Fit a semiparametric density ratio model (DRM) to m+1 (m>=1)
  samples using maximum dual empirical likelihood method.

  Denote the population cumulative distribution functions of the
  m+1 samples as \eqn{F_k(x)}{F_k(x)}'s,
  \eqn{k = 0, \, 1, \, \ldots, \, m}{k = 0, 1, ..., m}.
  We pick \eqn{F_0(x)}{F_0(x)} as a baseline distribution. The DRM
  assumes that the ratio of the density of each non-baseline
  distribution to the density of the baseline distribution
  satisfies
  \deqn{dF_k(x)/dF_0(x) = \exp(\alpha + \beta^T q(x)), \  k=1, \,
    \ldots, \, m}{dF_k(x)/dF_0(x) = exp(alpha + beta^T*q(x)), k =
    1, ..., m}
  where \eqn{q(x)}{q(x)} is a pre-specified d-dimensional basis
  function of data, and \eqn{\alpha}{alpha}, \eqn{\beta}{beta} are
  model parameters. No parametric form for baseline distribution
  \eqn{F_0}{F_0} is assumed.
}
\usage{
drmdel(x, n_samples, basis_func, par_init=NULL, g_null=NULL,
       g_null_gr=NULL, par_dim_null=NULL, par_init_null=NULL, ...)
}
\arguments{
  \item{x}{a vector formed by concatenating multiple samples,
    \eqn{x_0}{x_0}, \eqn{x_1}{x_1}, ..., \eqn{x_m}{x_m}, in the
    order of baseline sample (\eqn{x_0}{x_0}), non-baseline sample
    1 (\eqn{x_1}{x_1}), ..., non-baseline sample m
    (\eqn{x_m}{x_m}).}
  \item{n_samples}{a vector of length m+1 specifying the sizes of
    the multiple samples, in the order of 0, 1, ..., m.}
  \item{basis_func}{basis function q(x) of the DRM; must either be
    an integer between 1 and 11 or a function of the data, x. The
    integers represents built-in basis-functions:

           1 -- \eqn{q(x) = x}{q(x) = x}.

           2 -- \eqn{q(x) = \log(|x|)}{q(x) = log(|x|)}.

           3 -- \eqn{q(x) = \sqrt{|x|}}{q(x) = sqrt(|x|)}.

           4 -- \eqn{q(x) = x^2}{q(x) = x^2}.

           5 -- \eqn{q(x) = (x, \, x^2)}{q(x) = (x, x^2)};
           Normal model.

           6 -- \eqn{q(x) = (x, \, \log(|x|))}{
             q(x) = (x, log(|x|))}; Gamma model.

           7 -- \eqn{q(x) = (\log(|x|), \, \sqrt{|x|}, \, x)}{
             q(x) = (log(|x|), sqrt(|x|), x)}.

           8 -- \eqn{q(x) = (\log(|x|), \, \sqrt{|x|}, \, x^2)}{
             q(x) = (log(|x|), sqrt(|x|), x^2)}.

           9 -- \eqn{q(x) = (\log(|x|), \, x, \, x^2)}{
             q(x) = (log(|x|), x, x^2)}.

           10 -- \eqn{q(x) = (\sqrt{|x|}, \, x, \, x^2)}{
             q(x) = (sqrt(|x|), x, x^2)}.

           11 -- \eqn{q(x) = (\log(|x|), \, \sqrt{|x|}, \, x, \,
             x^2)}{ q(x) = (log(|x|), sqrt(|x|), x, x^2)}.

    If the basis function one wants to use is in the above list,
    one should use the built-in function to maximize the speed of
    model fitting.}
  \item{par_init}{a vector of length m*(d+1) specifying the
    initial value of the parameter vector; if not specified, all
    set to zeors. In fact, it is better to always set initial
    values to zeros, i.e. default, for all parameters to ensure
    that, at initial values, the dual log-likelihood is finite!}
  \item{g_null}{the function specifying the null hypothesis about
    DRM parameter \eqn{\beta}{beta} if there is one; The default
    is NULL.}
  \item{g_null_gr}{a funciton specifying the gradient of g_null,
    which must return a matrix of dimension m*(d+1) by
    dim(par_null), if available. The default is NULL.}
  \item{par_dim_null}{dimension of the parameter vector in null
    hypothesis if there is one. The default is NULL.}
  \item{par_init_null}{a vector of length par_dim_null specifying
    the initial value of the parameter vector under null; if not
    specified (default), all set to zeros (recommended).}
  \item{...}{further arguments to be passed to the R function
    \code{optim} for maximizing the dual empirical likelihood. See
    \code{help(optim)} for details. In the \code{\link{drmdel}}
    function, by default, the "control$method" and "control$maxit"
    arguments of \code{optim} are set to "BFGS" and 10000,
    respectively.
  }
}
\value{
  \item{drm_info}{a list of basic information about the DRM:

            m -- number of samples - 1.

            d -- dimension of the basis function.

            n_smaples -- the input vector of length m+1 specifying
            the size of each sample.

            n_total -- total sample size.

            basis_func -- the input basis function of the DRM.

            rho -- sample proportion: n_samples/n_total.
  }
  \item{mdele}{maximum dual empirical likelihood estimator (MDELE)
    of the model parameters. The output is a vector organized in
    the following form:
    \deqn{
    (\alpha_1, \, \beta_{1,1}, \, \beta_{1,2}, \, ..., \,
     \beta_{1,d}, \,
     \alpha_2, \, \beta_{2,1}, \, \beta_{2,2}, \, ..., \,
     \beta_{2,d}, \,
     ..., \,
     \alpha_m, \, \beta_{m,1}, \, \beta_{m,2}, \, ..., \,
     \beta_{m,d}).
    }{
    (alpha[1], beta[1,1], beta[1,2], ..., beta[1,d],
     alpha[2], beta[2,1], beta[2,2], ..., beta[2,d],
     ...,
     alpha[m], beta[m,1], beta[m,2], ..., beta[m,d]).
    }
  }
  \item{info_mat}{estimated information matrix.}
  \item{negldl}{negative log dual-likelihood evaluated at mdele.}
  \item{mdele_null}{mdele of the parameters under the null
    hypothesis, if available.}
  \item{negldl_null}{negative log dual-likelihood evaluated at
    mdele under the null hypothesis, if available.}
  \item{delr}{dual-empirical-likelihood ratio statistic evaluated
    under the null hypothesis. If no null hypotheis (g_null) is
    given, this is simply -2*negldl.}
  \item{df}{degrees of freedom of the chi-square limiting
    distribution for DELR statistic under the null.}
  \item{p_val}{p-vale of the DELR test.}
  \item{p_est}{estimated \eqn{dF_k(x)}{dF_k(x)}'s at the observed
    data points, under the DRM. This is a data frame with the
    following three columns:

            k -- label for the populations, k = 0, 1, ..., m.

            x -- data points; at which 'x' value
            \eqn{dF_k(x)}{dF_k(x)} is estimated.

            p_est -- estimated \eqn{dF_k(x)}{dF_k(x)}.

    NOTE: To estimate the density of \eqn{F_k(x)}{F_k(x)}, it is
    recommended to use \code{\link{densityDRM}} function.
  }
  \item{cdf_est}{estimated CDFs, \eqn{F_k(x)}{F_k(x)}'s, at the
    observed data points, under the DRM. This is a data frame with
    the following three columns:

            k -- label for the populations, k = 0, 1, ..., m.

            x -- data points; at which 'x' value
            \eqn{F_k(x)}{F_k(x)} is estimated.

            cdf_est -- estimated \eqn{F_k(x)}{F_k(x)}.

    NOTE: To estimate CDF \eqn{dF_k(x)}{dF_k(x)}, it is
    recommended to use \code{\link{cdfDRM}} function instead of
    looking at this output.
  }
}
\references{
  S. Cai, J. Chen and J. V. Zidek (2013),
  Dual-empirical-likelihood ratio test under density ratio models
  for multiple samples.
  Manuscript.

  A. Keziou and S. Leoni-Aubin (2008),
  On empirical likelihood for semiparametric two-sample density
  ratio models.
  \emph{Journal of Statistical Planning and Inference},
  138:915-928.
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
# is the basis function for gamma family.

# There are 11 built-in basis function in drmdel(). And q(x) = (x,
# log(abs(x))) is the 6th basis function, so we can fit the model
# by specifying basis_func=6 in drmdel() as follows:
drmfit <- drmdel(x=x, n_samples=n_samples, basis_func=6)
names(drmfit)

# A brief summary of the DRM fit
summaryDRM(drmfit)

# Another way of specifying basis function for drmdel() is to pass
# a user-specified R function to the basis_func argument of the
# drmdel() function.
# NOTE: If the basis function one wants to use is included in the
# built-in function list, one should use the built-in functions by
# passing an integer between 1 to 11 to the drmdel() function,
# because the computation will be faster with a built-in function
# than with a user-specified function.
basis_gamma <- function(x) return(c(x, log(abs(x))))
drmfit1 <- drmdel(x=x, n_samples=n_samples,
                  basis_func=basis_gamma)

# One can see the summary of this DRM fit is exactly the same as
# that of the previous fit with basis_func=6
summaryDRM(drmfit1)
}
\keyword{models}
\keyword{methods}
