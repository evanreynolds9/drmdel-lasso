/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3, April 08, 2014.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "basisFuncs.h"
#include "utilities.h"

void lp_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    double *restrict lp /*outputs*/)
/* lp (linear predictor) is a vector of length m for the m non-baseline
 *   samples. */
{

  /* loop indices */
  unsigned long i, j;

  /* calculating lp[i] */
  for (i = 0; i < m; ++i) {
    lp[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      lp[i] = lp[i] + par_mat[i][j] * h[j-1];
    }
  }

}

void r_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    double *restrict r /*outputs*/)
/* r is a vector of length m for the m non-baseline samples. */
{

  /* loop indices */
  unsigned long i, j;

  /* calculating r_i */
  for (i = 0; i < m; ++i) {
    r[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      r[i] = r[i] + par_mat[i][j] * h[j-1];
    }
    r[i] = exp(r[i]);
  }

}

void R_val(unsigned long m, unsigned long d, double * restrict h, /*inputs*/
    double *restrict* restrict par_mat, double * restrict n_samples, /*inputs*/
    double *restrict R /*outputs*/)
/* R is a vector of length m for the m non-baseline samples.
 * R[i] = n_samples[i] * exp(alpha_i + beta_i^T h(x)), i = 1, 2, ..., m*/
{

  /* loop indices */
  unsigned long i, j;

  /* calculating r_i */
  for (i = 0; i < m; ++i) {
    R[i] = par_mat[i][0];
    for (j = 1; j < d + 1; ++j) {
      R[i] = R[i] + par_mat[i][j] * h[j-1];
    }
    R[i] = n_samples[i+1] * exp(R[i]);
  }

}

double logDualL(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double *restrict* restrict x_mat /*inputs*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k;

  double * restrict lp;
  lp = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (lp == NULL) errMsg("malloc() allocation failure for lp!");
  for (i = 0; i < m; ++i) {
    lp[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;

  /*define output*/
  double ldl_val;
  ldl_val = 0.0;

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      (*h_func)(x_mat[i][j], h); /*update h*/

      lp_val(m, d, h, par_mat, lp); /*update lp*/

      /* calculating q_i */
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += rho[k+1] * exp(lp[k]);
      }

      if (i == 0) {
        ldl_val = ldl_val - log(S);  
      } else {
        ldl_val = ldl_val + lp[i-1] - log(S);
      }

    }

  }

  /* free arrays */
  free((void *) lp);
  free((void *) h);
  free((void *) rho);

  return ldl_val;

}

void logDualLWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict ldl_val /*output*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;

  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1x, x_mat);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1logx, x_mat);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1sqrtx, x_mat);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 1,
          par_mat, &h1xSquare, x_mat);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 2,
          par_mat, &h2Normal, x_mat);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 2,
          par_mat, &h2Gamma, x_mat);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3a, x_mat);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3b, x_mat);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3c, x_mat);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 3,
          par_mat, &h3d, x_mat);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
          n_samples_use, (unsigned long)*m, 4,
          par_mat, &h4a, x_mat);
      break;
    
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      *ldl_val = logDualL((unsigned long)*n_total,
                        /*n_samples_use, (unsigned long)*m, (unsigned long)*d,*/
                        n_samples_use, (unsigned long)*m, 5,
                        par_mat, &h5a, x_mat);
    break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);

}

void logDualLGr(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double *restrict* restrict x_mat, /*inputs*/
    double *restrict* restrict ldl_gr_mat /*outputs*/)
/* Calculating the gradient of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr_mat -- a 2-D pointer array of dimension m by (d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k, l;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;
  double tmp_double;

  /*initializing ldl_gr_mat as safeguard*/
  for (i = 0; i < m; ++i) {
    for (j = 0; j < (d+1); ++j) {
      ldl_gr_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update H = (1, h^T)^T*/
      (*h_func)(x_mat[i][j], H+1); /*update H*/

      R_val(m, d, H+1, par_mat, rho, R); /*update R*/

      /*calculating S*/
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += R[k];
      }

      /* calculating the gradient of ldl */
      for (k = 0; k < m; ++k) {

        tmp_double = -R[k]/S;

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[k][l] += tmp_double * H[l];
        }

      }

      if (i > 0) {

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[i-1][l] = H[l] + ldl_gr_mat[i-1][l];
        }

      }

    }

  }

  /* free arrays */
  free((void *) R);
  free((void *) H);
  free((void *) rho);

}

void logDualLGrWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict ldl_gr /*output*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr -- a vector of length m*(d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;
  double *restrict* restrict ldl_gr_mat;

  /* converting x, par and ldl_gr to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_gr_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (ldl_gr_mat == NULL){
    errMsg("malloc() allocation failure for ldl_gr_mat!");
  }
  ldl_gr_mat[0] = ldl_gr;
  for (i = 1; i < (unsigned long)*m; ++i){
    ldl_gr_mat[i] = ldl_gr_mat[i-1] + ((unsigned long)*d + 1);
  }


  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1x, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1logx, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1sqrtx, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1xSquare, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Normal, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Gamma, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3a, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3b, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3c, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3d, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
          /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
          n_samples_use, (unsigned long)*m, 4, /*inputs*/
          par_mat, &h4a, x_mat, /*inputs*/
          ldl_gr_mat /*outputs*/);
      break;
      
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      logDualLGr((unsigned long)*n_total, /*inputs*/
        /*n_samples_use, (unsigned long)*m, (unsigned long)*d, [>inputs<]*/
        n_samples_use, (unsigned long)*m, 5, /*inputs*/
        par_mat, &h5a, x_mat, /*inputs*/
        ldl_gr_mat /*outputs*/);
    break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);
  free((void *) ldl_gr_mat);

}

void logDualLHessian(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double * restrict x, /*inputs*/
    double *restrict* restrict ldl_hessian_mat /*outputs*/)
/* Calculating the Hessian of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_hessian_mat -- the m(d+1) by m(d+1) Hessian matrix evaluated at the
 *     parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i, k, l, kh, lh;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }

  /*qaa matrix*/
  double *restrict* restrict qaa;
  qaa = (double *restrict* restrict) malloc((size_t) (m*sizeof(double*)));
  if (qaa == NULL) errMsg("malloc() allocation failure for qaa!");

  qaa[0] = (double * restrict) malloc((size_t) ((m*m)*sizeof(double)));
  if (qaa[0] == NULL) errMsg("malloc() allocation failure for qaa[0]!");
  for(i = 1; i < m; ++i) {
    qaa[i] = qaa[i-1] + m;
  }

  for (i = 0; i < m; ++i) {
    for (k = 0; k < m; ++k) {
      qaa[i][k] = 0.0;
    }
  }

  /*other variables*/
  double S;

  /*initializing ldl_hessian_mat as safeguard*/
  for (i = 0; i < m*(d+1); ++i) {
    for (k = 0; k < m*(d+1); ++k) {
      ldl_hessian_mat[i][k] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update H = (1, h^T)^T*/
    (*h_func)(x[i], H+1); /*update H*/

    R_val(m, d, H+1, par_mat, n_samples, R); /*update R*/

    /*calculating S*/
    S = n_samples[0];
    for (k = 0; k < m; ++k) {
      S += R[k];
    }

    for (k = 0; k < m; ++k) {

      for (l = 0; l < m; ++l) {
        qaa[k][l] = R[k]*R[l]/(S*S);
      }

    }
    for (k = 0; k < m; ++k) {
        qaa[k][k] -= R[k]/S;
    }

    /*calculating the Hessian matrix m(d+1) by m(d+1)*/
    /*set very entry of ldl_hessian_mat to be 0;*/
    /*k and l are for vertical (row) indices; kh and lh are for horizontal
     * (column) indices.*/
    for (k = 0; k < m; ++k) {
      for (l = 0; l < (d+1); ++l) {

        for (kh = 0; kh < m; ++kh) {
          for (lh = 0; lh < (d+1); ++lh) {

            ldl_hessian_mat[(k*(d+1) + l)][(kh*(d+1) + lh)] +=
              qaa[k][kh]*H[l]*H[lh];

          }
        }

      }
    }

  }

  /* free arrays */
  free((void *) R);
  free((void *) H);

  free((void *) qaa[0]);
  free((void *) qaa);

}

void logDualLHessianWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d, /*inputs*/
    double * restrict par, double * restrict model, /*inputs*/
    double * restrict x, /*inputs*/
    double * restrict ldl_hessian /*output*/)
/* Calculating the Hessian of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1,
 *     \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m;
 * Outputs:
 *   ldl_hessian -- the vector form (row by row) of the m(d+1) by m(d+1)
 *     Hessian matrix evaluated at the parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  double *restrict* restrict ldl_hessian_mat;

  /* converting par and ldl_hessian to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  ldl_hessian_mat = (double *restrict* restrict) malloc((size_t)
      ((((unsigned long)*m) * ((unsigned long)*d + 1))*sizeof(double*)));
  if (ldl_hessian_mat == NULL){
    errMsg("malloc() allocation failure for ldl_hessian_mat!");
  }
  ldl_hessian_mat[0] = ldl_hessian;
  for (i = 1; i < (((unsigned long)*m) * ((unsigned long)*d + 1)); ++i){
    ldl_hessian_mat[i] = ldl_hessian_mat[i-1] + 
      (((unsigned long)*m) * ((unsigned long)*d + 1));
  }

  /* calculating log dual empirical likelihood value at 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1x, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1logx, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1sqrtx, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 1, /*inputs*/
          par_mat, &h1xSquare, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Normal, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 2, /*inputs*/
          par_mat, &h2Gamma, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3a, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3b, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3c, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 3, /*inputs*/
          par_mat, &h3d, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long)*m, 4, /*inputs*/
          par_mat, &h4a, x, /*inputs*/
          ldl_hessian_mat /*outputs*/);
      break;
      
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      logDualLHessian((unsigned long)*n_total, n_samples, /*inputs*/
        (unsigned long)*m, 5, /*inputs*/
        par_mat, &h5a, x, /*inputs*/
        ldl_hessian_mat /*outputs*/);
      break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) par_mat);
  free((void *) ldl_hessian_mat);

}

void probEst(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    void (*h_func)(double, double * restrict), /*input*/
    double * restrict x, /*inputs*/
    unsigned short normalize,
    double *restrict* restrict pEst_mat /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst_mat -- a (m+1) by n_total matrix; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i, j;

  double * restrict r;
  r = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (r == NULL) errMsg("malloc() allocation failure for r!");
  for (i = 0; i < m; ++i) {
    r[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }

  double p_bl_tmp;

  /* auxiliary variable, used only when normalize == 1 */
  double * restrict pEst_sum;
  pEst_sum = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (pEst_sum == NULL) errMsg("malloc() allocation failure for pEst_sum!");
  for (i = 0; i < m+1; ++i) {
    pEst_sum[i] = 0.0;
  }

  /* initializing output variable pEst_mat to 0 as safeguard */
  for (i = 0; i < (m + 1); ++i) {
    for (j = 0; j < n_total; ++j) {
        pEst_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    (*h_func)(x[i], h); /*update h*/

    r_val(m, d, h, par_mat, r); /*update r*/

    /* calculate baseline probabilities */
    p_bl_tmp = n_samples[0];
    for (j = 0; j < m; ++j) {
      p_bl_tmp += n_samples[j+1] * r[j];
    }
    p_bl_tmp = 1.0/p_bl_tmp;

    /* calculate probabilities for all m+1 sampels*/
    pEst_mat[0][i] = p_bl_tmp; /*baseline probablities*/
    if (normalize==1) pEst_sum[0] += pEst_mat[0][i];
    for (j = 1; j < m+1; ++j) {
      pEst_mat[j][i] = r[j-1]*p_bl_tmp;
      if (normalize==1) pEst_sum[j] += pEst_mat[j][i];
    }

  }

  if (normalize==1) {

    for (i = 0; i < m+1; ++i) {

      for (j = 0; j < n_total; ++j) {
          pEst_mat[i][j] /= pEst_sum[i];
      }

    }

  }

  /* free arrays */
  free((void *) r);
  free((void *) h);
  free((void *) pEst_sum);

}

void probEstWrapper(double * restrict n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    double * restrict m, double * restrict d,
    double * restrict par, /*inputs*/
    double * restrict model, double * restrict x, /*inputs*/
    double * restrict normalize,
    double * restrict pEst /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst -- a vector of length (m+1)*n_total; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par;
  for (i = 1; i < (unsigned long)*m; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d + 1);
  }

  double *restrict* restrict pEst_mat;
  pEst_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m+1)*sizeof(double*)));
  if (pEst_mat == NULL) errMsg("malloc() allocation failure for pEst_mat!");
  pEst_mat[0] = pEst;
  for (i = 1; i < ((unsigned long)*m + 1); ++i){
    pEst_mat[i] = pEst_mat[i-1] + (unsigned long)*n_total;
  }

  /* calculating estimated probability for baseline at a given 'par' */
  switch ((unsigned long)*model)
  {
    case 1 :
      /*printf("h(x) = (x)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 1, h(x) = x, d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1x, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 2 :
      /*printf("h(x) = (log(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 2, h(x) = log(x), d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1logx, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 3 :
      /*printf("h(x) = (sqrt(x))");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 3, h(x) = sqrt(x), d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1sqrtx, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 4 :
      /*printf("h(x) = (x^2)");*/
      if ((unsigned long)*d != 1) {
        errMsg("For model 4, h(x) = x^2, d must be 1!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 1, /*inputs*/
          par_mat, &h1xSquare, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 5 :
      /*printf("h(x) = (x, x^2) -- Normal model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 5 (Normal model), h(x) = (x, x^2), d must be 2!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 2, /*inputs*/
          par_mat, &h2Normal, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 6 :
      /*printf("h(x) = (x, log(x)) -- Gamma model");*/
      if ((unsigned long)*d != 2) {
        errMsg("For model 6 (Gamma model), h(x) = (x, log(x)), d must be 2!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 2, /*inputs*/
          par_mat, &h2Gamma, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 7 :
      /*printf("h(x) = (log(x), sqrt(x), x)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 7, h(x) = (log(x), sqrt(x), x), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3a, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 8 :
      /*printf("h(x) = (log(x), sqrt(x), x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 8, h(x) = (log(x), sqrt(x), x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3b, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 9 :
      /*printf("h(x) = (log(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 9, h(x) = (log(x), x, x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3c, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 10 :
      /*printf("h(x) = (sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 3) {
        errMsg("For model 10, h(x) = (sqrt(x), x, x^2), d must be 3!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 3, /*inputs*/
          par_mat, &h3d, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;

    case 11 :
      /*printf("h(x) = (log(x), sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 4) {
        errMsg("For model 11, h(x) = (log(x), sqrt(x), x, x^2), d must be 4!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
          (unsigned long) *m, 4, /*inputs*/
          par_mat, &h4a, x, (unsigned short) *normalize, /*inputs*/
          pEst_mat /*output*/);
      break;
    
    case 12 :
      /*printf("h(x) = (log(x), log(x)^2, sqrt(x), x, x^2)");*/
      if ((unsigned long)*d != 5) {
        errMsg("For model 11, h(x) = (log(x), log(x)^2, sqrt(x), x, x^2), d must be 5!");
      }
      probEst((unsigned long)*n_total, n_samples, /*inputs*/
        (unsigned long) *m, 5, /*inputs*/
        par_mat, &h5a, x, (unsigned short) *normalize, /*inputs*/
        pEst_mat /*output*/);
      break;

    default :
      errMsg("'Model' must be an integer between 1 and 11 or a function of a single data point");
      break;

  }

  /* free arrays */
  free((void *) par_mat);
  free((void *) pEst_mat);

}

/* ***** User specified basis function (Uf) version ***** */

double logDualLUf(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double *restrict* restrict x_mat /*inputs*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k;

  double * restrict lp;
  lp = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (lp == NULL) errMsg("malloc() allocation failure for lp!");
  for (i = 0; i < m; ++i) {
    lp[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;

  /*define output*/
  double ldl_val;
  ldl_val = 0.0;

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update h*/
      REAL(h_args)[0] = x_mat[i][j];
      PROTECT(h_r_call = lang2(h_func, h_args));
      PROTECT(h_r_call_result = eval(h_r_call, env));
      for (k = 0; k < d; ++k) {
        h[k] = REAL(h_r_call_result)[k];
      }
      UNPROTECT(2);

      lp_val(m, d, h, par_mat, lp); /*update lp*/

      /* calculating q_i */
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += rho[k+1] * exp(lp[k]);
      }

      if (i == 0) {
        ldl_val = ldl_val - log(S);  
      } else {
        ldl_val = ldl_val + lp[i-1] - log(S);
      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) lp);
  free((void *) h);
  free((void *) rho);

  return ldl_val;

}

SEXP logDualLUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
    SEXP h_func, SEXP env, SEXP x)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_val -- value of ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* loop indices */
  unsigned long i;

  double ldl_val_c;
  SEXP ldl_val;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m_c + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m_c + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples_c[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;

  /* converting x and par to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m_c+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_val_c = logDualLUf((unsigned long)*n_total_c, n_samples_use,
      (unsigned long)*m_c, (unsigned long)*d_c, par_mat,
      h_func, env, x_mat);

  PROTECT(ldl_val = allocVector(REALSXP, 1));
  REAL(ldl_val)[0] = ldl_val_c;

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);

  UNPROTECT(1);

  return(ldl_val);

}

void logDualLGrUf(unsigned long n_total, /*inputs*/
    unsigned long * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double *restrict* restrict x_mat, /*inputs*/
    double *restrict* restrict ldl_gr_mat /*outputs*/)
/* Calculating the gradient of log dual empirical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_gr_mat -- a 2-D pointer array of dimension m by (d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
 */
{
  /* loop indices */
  unsigned long i, j, k, l;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double * restrict rho;
  rho = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (rho == NULL) errMsg("malloc() allocation failure for rho!");
  for (i = 0; i < m+1; ++i) {
    rho[i] = (double)n_samples[i]/(double)n_total;
  }

  /*other variables*/
  double S;
  double tmp_double;

  /*initializing ldl_gr_mat as safeguard*/
  for (i = 0; i < m; ++i) {
    for (j = 0; j < (d+1); ++j) {
      ldl_gr_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < m+1; ++i) {

    for (j = 0; j < n_samples[i]; ++j) {

      /*update H = (1, h^T)^T*/
      REAL(h_args)[0] = x_mat[i][j];
      PROTECT(h_r_call = lang2(h_func, h_args));
      PROTECT(h_r_call_result = eval(h_r_call, env));
      for (k = 0; k < d; ++k) {
        H[k+1] = REAL(h_r_call_result)[k];
      }
      UNPROTECT(2);

      R_val(m, d, H+1, par_mat, rho, R); /*update R*/

      /*calculating S*/
      S = rho[0];
      for (k = 0; k < m; ++k) {
        S += R[k];
      }

      /* calculating the gradient of ldl */
      for (k = 0; k < m; ++k) {

        tmp_double = -R[k]/S;

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[k][l] += tmp_double * H[l];
        }

      }

      if (i > 0) {

        for (l = 0; l < (d+1); ++l) {
          ldl_gr_mat[i-1][l] = H[l] + ldl_gr_mat[i-1][l];
        }

      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) R);
  free((void *) H);
  free((void *) rho);

}

SEXP logDualLGrUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* Outputs:
   *   ldl_gr -- a vector of length m*(d+1); value of the gradient of the ldl (log dual empirical likelihood) at a given "par" value.
   */
  SEXP ldl_gr;
  PROTECT(ldl_gr = allocVector(REALSXP, ((unsigned long)*m_c) * ((unsigned
            long)*d_c + 1)));
  double * restrict ldl_gr_c;
  ldl_gr_c = REAL(ldl_gr);

  /* loop indices */
  unsigned long i;

  unsigned long * restrict n_samples_use;
  n_samples_use = (unsigned long * restrict) malloc((size_t) (((unsigned
            long)*m_c + 1)*sizeof(unsigned long)));
  if (n_samples_use == NULL) {
    errMsg("malloc() allocation failure for m_samples_use!");
  }
  for (i = 0; i < ((unsigned long)*m_c + 1); ++i) {
    n_samples_use[i] = (unsigned long)n_samples_c[i];
  }

  double *restrict* restrict par_mat;
  double *restrict* restrict x_mat;
  double *restrict* restrict ldl_gr_mat;

  /* converting x, par and ldl_gr to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  x_mat = (double *restrict* restrict) malloc((size_t) (((unsigned
            long)*m_c+1)*sizeof(double*)));
  if (x_mat == NULL) errMsg("malloc() allocation failure for x_mat!");
  x_mat[0] = x_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    x_mat[i] = x_mat[i-1] + n_samples_use[i - 1];
  }

  ldl_gr_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (ldl_gr_mat == NULL){
    errMsg("malloc() allocation failure for ldl_gr_mat!");
  }
  ldl_gr_mat[0] = ldl_gr_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    ldl_gr_mat[i] = ldl_gr_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  logDualLGrUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_use, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_mat, /*inputs*/
      ldl_gr_mat /*outputs*/);

  /* free arrays */
  free((void *) n_samples_use);
  free((void *) x_mat);
  free((void *) par_mat);
  free((void *) ldl_gr_mat);

  UNPROTECT(1);

  return(ldl_gr);

}

void logDualLHessianUf(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*input*/
    double * restrict x, /*inputs*/
    double *restrict* restrict ldl_hessian_mat /*outputs*/)
/* Calculating the Hessian of log dual emprical likelihood (+ n \log n) at a given
 * parameter value.
 * Inputs:
 *   n_total -- length of data;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 * Outputs:
 *   ldl_hessian_mat -- the m(d+1) by m(d+1) Hessian matrix evaluated at the
 *     parameter value par_mat.
 */
{
  /* loop indices */
  unsigned long i, k, l, kh, lh;

  double * restrict R;
  R = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (R == NULL) errMsg("malloc() allocation failure for R!");
  for (i = 0; i < m; ++i) {
    R[i] = 0.0;
  }

  double * restrict H;
  H = (double * restrict) malloc((size_t) ((d+1)*sizeof(double)));
  if (H == NULL) errMsg("malloc() allocation failure for H!");
  H[0] = 1.0;
  for (i = 1; i < (d+1); ++i) {
    H[i] = 0.0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  /*qaa matrix*/
  double *restrict* restrict qaa;
  qaa = (double *restrict* restrict) malloc((size_t) (m*sizeof(double*)));
  if (qaa == NULL) errMsg("malloc() allocation failure for qaa!");

  qaa[0] = (double * restrict) malloc((size_t) ((m*m)*sizeof(double)));
  if (qaa[0] == NULL) errMsg("malloc() allocation failure for qaa[0]!");
  for(i = 1; i < m; ++i) {
    qaa[i] = qaa[i-1] + m;
  }

  for (i = 0; i < m; ++i) {
    for (k = 0; k < m; ++k) {
      qaa[i][k] = 0.0;
    }
  }

  /*other variables*/
  double S;

  /*initializing ldl_hessian_mat as safeguard*/
  for (i = 0; i < m*(d+1); ++i) {
    for (k = 0; k < m*(d+1); ++k) {
      ldl_hessian_mat[i][k] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update H = (1, h^T)^T*/
    REAL(h_args)[0] = x[i];
    PROTECT(h_r_call = lang2(h_func, h_args));
    PROTECT(h_r_call_result = eval(h_r_call, env));
    for (k = 0; k < d; ++k) {
      H[k+1] = REAL(h_r_call_result)[k];
    }
    UNPROTECT(2);

    R_val(m, d, H+1, par_mat, n_samples, R); /*update R*/

    /*calculating S*/
    S = n_samples[0];
    for (k = 0; k < m; ++k) {
      S += R[k];
    }

    for (k = 0; k < m; ++k) {

      for (l = 0; l < m; ++l) {
        qaa[k][l] = R[k]*R[l]/(S*S);
      }

    }
    for (k = 0; k < m; ++k) {
        qaa[k][k] -= R[k]/S;
    }

    /*calculating the Hessian matrix m(d+1) by m(d+1)*/
    /*set very entry of ldl_hessian_mat to be 0;*/
    /*k and l are for vertical (row) indices; kh and lh are for horizontal
     * (column) indices.*/
    for (k = 0; k < m; ++k) {
      for (l = 0; l < (d+1); ++l) {

        for (kh = 0; kh < m; ++kh) {
          for (lh = 0; lh < (d+1); ++lh) {

            ldl_hessian_mat[(k*(d+1) + l)][(kh*(d+1) + lh)] += 
              qaa[k][kh]*H[l]*H[lh];

          }
        }

      }
    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) R);
  free((void *) H);

  free((void *) qaa[0]);
  free((void *) qaa);

}


SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)
/* Calculating log dual empirical likelihood (+ n \log n) at a given parameter value.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par -- values of parameters (length of m(d+1)), organized as \theta_1, \cdots, \theta_m;
 *   h_func -- the basis function of the DRM;
 *   x -- data, organized as x_0, x_1, ..., x_m.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);

  /* Outputs:
   *   ldl_hessian -- a m(d+1) x m(d+1) vector for hessian matrix of ldl at the parameter value par.
   */
  SEXP ldl_hessian;
  PROTECT(ldl_hessian = allocMatrix(REALSXP,
        ((unsigned long)*m_c) * ((unsigned long)*d_c + 1),
        ((unsigned long)*m_c) * ((unsigned long)*d_c + 1)
        ));
  double * restrict ldl_hessian_c;
  ldl_hessian_c = REAL(ldl_hessian);

  /* loop indices */
  unsigned long i;

  double *restrict* restrict par_mat;
  double *restrict* restrict ldl_hessian_mat;

  /* converting x, par and ldl_hessian to matrices */
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  ldl_hessian_mat = (double *restrict* restrict) malloc((size_t)
      ((((unsigned long)*m_c) * ((unsigned long)*d_c + 1))*sizeof(double*)));
  if (ldl_hessian_mat == NULL){
    errMsg("malloc() allocation failure for ldl_hessian_mat!");
  }
  ldl_hessian_mat[0] = ldl_hessian_c;
  for (i = 1; i < (((unsigned long)*m_c) * ((unsigned long)*d_c + 1)); ++i){
    ldl_hessian_mat[i] = ldl_hessian_mat[i-1] + 
      (((unsigned long)*m_c) * ((unsigned long)*d_c + 1));
  }

  logDualLHessianUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_c, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_c, /*inputs*/
      ldl_hessian_mat /*outputs*/);

  /* free arrays */
  free((void *) par_mat);
  free((void *) ldl_hessian_mat);

  UNPROTECT(1);

  return(ldl_hessian);

}

void probEstUf(unsigned long n_total, /*inputs*/
    double * restrict n_samples, /*inputs*/
    unsigned long m, unsigned long d, /*inputs*/
    double *restrict* restrict par_mat, /*inputs*/
    SEXP h_func, SEXP env, /*inputs*/
    double * restrict x, /*inputs*/
    unsigned short normalize, /*inputs*/
    double *restrict* restrict pEst_mat /*output*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a (double, not unsigned long!) vector of length m+1
 *     specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x -- data, a long vector in order of as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 * Outputs:
 *   pEst_mat -- a (m+1) by n_total matrix; estimated probability for x's at
 *     parameter value 'par_mat' for each sample.
 */
{
  /* loop indices */
  unsigned long i, j;

  double * restrict r;
  r = (double * restrict) malloc((size_t) (m*sizeof(double)));
  if (r == NULL) errMsg("malloc() allocation failure for r!");
  for (i = 0; i < m; ++i) {
    r[i] = 0;
  }

  double * restrict h;
  h = (double * restrict) malloc((size_t) (d*sizeof(double)));
  if (h == NULL) errMsg("malloc() allocation failure for h!");
  for (i = 0; i < d; ++i) {
    h[i] = 0;
  }
  SEXP h_args, h_r_call, h_r_call_result;
  PROTECT(h_args = allocVector(REALSXP, 1));

  double p_bl_tmp;

  /* auxiliary variable, used only when normalize == 1 */
  double * restrict pEst_sum;
  pEst_sum = (double * restrict) malloc((size_t) ((m+1)*sizeof(double)));
  if (pEst_sum == NULL) errMsg("malloc() allocation failure for pEst_sum!");
  for (i = 0; i < m+1; ++i) {
    pEst_sum[i] = 0.0;
  }

  /* initializing output variable pEst_mat to 0 as safeguard */
  for (i = 0; i < (m + 1); ++i) {
    for (j = 0; j < n_total; ++j) {
        pEst_mat[i][j] = 0.0;
    }
  }

  for (i = 0; i < n_total; ++i) {

    /*update h*/
    REAL(h_args)[0] = x[i];
    PROTECT(h_r_call = lang2(h_func, h_args));
    PROTECT(h_r_call_result = eval(h_r_call, env));
    for (j = 0; j < d; ++j) {
      h[j] = REAL(h_r_call_result)[j];
    }
    UNPROTECT(2);

    r_val(m, d, h, par_mat, r); /*update r*/

    /* calculate baseline probabilities */
    p_bl_tmp = n_samples[0];
    for (j = 0; j < m; ++j) {
      p_bl_tmp += n_samples[j+1] * r[j];
    }
    p_bl_tmp = 1.0/p_bl_tmp;

    /* calculate probabilities for all m+1 sampels*/
    pEst_mat[0][i] = p_bl_tmp; /*baseline probablities*/
    if (normalize==1) pEst_sum[0] += pEst_mat[0][i];
    for (j = 1; j < m+1; ++j) {
      pEst_mat[j][i] = r[j-1]*p_bl_tmp;
      if (normalize==1) pEst_sum[j] += pEst_mat[j][i];
    }

  }

  if (normalize==1) {

    for (i = 0; i < m+1; ++i) {

      for (j = 0; j < n_total; ++j) {
          pEst_mat[i][j] /= pEst_sum[i];
      }

    }

  }

  UNPROTECT(1);

  /* free arrays */
  free((void *) r);
  free((void *) h);
  free((void *) pEst_sum);

}

SEXP probEstUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
    SEXP par, SEXP h_func, SEXP env, SEXP x, SEXP normalize /*input*/)
/* Estimating \hat{p} for all samples in density ratio models. \hat{p}'s are
 * useful for estimating quantiles, distribution function and density of each
 * population.
 * Inputs:
 *   n_total -- total sample size;
 *   n_samples -- a vector of length m+1 specifying the size of each sample;
 *   m -- number of samples - 1;
 *   d -- dimension of h(x); 
 *   par_mat -- values of parameters (pointer matrix of dimension m by (d+1));
 *   h_func -- the basis function of the DRM;
 *   x_mat -- 2-D pointer array of data, organized as x_0, x_1, ..., x_m.
 *   normalize -- indicator; wether normalize the probabilies so that the sum
 *     is strictly one: 1 for yes, 0 for no.
 */
{
  double * restrict n_total_c;
  n_total_c = REAL(n_total);
  double * restrict n_samples_c;
  n_samples_c = REAL(n_samples);
  double * restrict m_c;
  m_c = REAL(m);
  double * restrict d_c;
  d_c = REAL(d);
  double * restrict par_c;
  par_c = REAL(par);
  double * restrict x_c;
  x_c = REAL(x);
  double * restrict normalize_c;
  normalize_c = REAL(normalize);

  /* Outputs:
   *   pEst -- a vector of length (m+1)*n_total; estimated probability for x's at parameter value 'par_mat' for each sample.
   */
  SEXP pEst;
  PROTECT(pEst = allocVector(REALSXP, ((unsigned long)*m_c + 1) * (unsigned
          long)*n_total_c));
  double * restrict pEst_c;
  pEst_c = REAL(pEst);

  /* loop indices */
  unsigned long i;

  /* converting x and par to matrices */
  double *restrict* restrict par_mat;
  par_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c)*sizeof(double*)));
  if (par_mat == NULL) errMsg("malloc() allocation failure for par_mat!");
  par_mat[0] = par_c;
  for (i = 1; i < (unsigned long)*m_c; ++i){
    par_mat[i] = par_mat[i-1] + ((unsigned long)*d_c + 1);
  }

  double *restrict* restrict pEst_mat;
  pEst_mat = (double *restrict* restrict) malloc((size_t)
      (((unsigned long)*m_c+1)*sizeof(double*)));
  if (pEst_mat == NULL) errMsg("malloc() allocation failure for pEst_mat!");
  pEst_mat[0] = pEst_c;
  for (i = 1; i < ((unsigned long)*m_c + 1); ++i){
    pEst_mat[i] = pEst_mat[i-1] + (unsigned long)*n_total_c;
  }

  probEstUf((unsigned long)*n_total_c, /*inputs*/
      n_samples_c, (unsigned long)*m_c, (unsigned long)*d_c, /*inputs*/
      par_mat, h_func, env, x_c, (unsigned short)*normalize_c, /*inputs*/
      pEst_mat /*outputs*/);

  /* free arrays */
  free((void *) par_mat);
  free((void *) pEst_mat);

  UNPROTECT(1);

  return(pEst);

}
