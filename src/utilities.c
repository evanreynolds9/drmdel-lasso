/* This software is written by Song Cai and published under GPLv3.
 *
 * Version 1.3.2, October 23, 2021.
 */

#include <stdlib.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include "utilities.h"

/* **********Error messages handling********** */
/*void errMsg(char err_text[])*/
/*[> Standard error handler <]*/
/*{*/
    /*fprintf(stderr, "Program run-time error:\n");*/
    /*fprintf(stderr, "    %s\n", err_text);*/
    /*fprintf(stderr, "...now exiting to system...\n");*/
    /*exit(1);*/
/*}*/

void errMsg(char err_text[])
/* Standard error handler; R version */
{
    error(err_text);
}
/* ******************** */

/* **********Matrix operations********** */
void kroneckerProd(double *restrict* restrict A, /*inputs*/
    unsigned long m, unsigned long n, /*inputs*/
    double *restrict* restrict B, /*inputs*/
    unsigned long p, unsigned long q, /*inputs*/
    double *restrict* restrict C /*output*/)
/* Kronecker product of A and B
 * Inputs:
 *   A -- a m by n matrix;
 *   B -- a p by q matrix.
 * Output:
 *   C -- a mp by nq matrix.
 */
{
  /* loop indices */
  unsigned long i, j, k, l;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < p; ++j) {

      for (k = 0; k < n; ++k) {
        for (l = 0; l < q; ++l) {

          C[(i*p + j)][(k*q + l)] = A[i][k]*B[j][l];

        }
      }

    }
  }

}
/* ******************** */

/* **********C wrappers to LAPACK functions********** */
void dgetrfCWrapper(double * restrict m, /*input*/
    double * restrict n, /*input*/
    double * restrict a, /*input and output*/
    double * restrict lda, /*input*/
    double * restrict ipiv, /*output*/
    double * restrict info /*output*/)
{
  /* LAPACK subroutine DGETRF (double precision) computes an
   * LU factorization # of a general M-by-N matrix A using
   * partial pivoting with row interchanges. It also
   * provides test for exact singularity.
   *
   * Fortran routine:
   * subroutine dgetrf(integer M, 
   *     integer N,
   *     double precision, dimension( lda, * ) A,
   *     integer LDA,
   *     integer, dimension( * ) IPIV,
   *     integer INFO)
   *
   * Input: M, N, A, LDA
   * Output: A, IPIV, INFO
   *
   * R LAPACK header (included in R_ext/Lapack.h)
   * F77_NAME(dgetrf)(const int* m, const int* n, double* a,
   *     const int* lda, int* ipiv, int* info);
   */

  /* loop indices */
  unsigned long i;

  int m_c;
  m_c = (int)*m;
  int n_c;
  n_c = (int)*n;
  int lda_c;
  lda_c = (int)*lda;
  int info_c;
  info_c = (int)*info;

  int ipiv_dim;
  if (m_c > n_c) {
    ipiv_dim = n_c;
  } else {
    ipiv_dim = m_c;
  }

  int * restrict ipiv_c;
  ipiv_c = (int * restrict) malloc((size_t) (ipiv_dim*sizeof(int)));
  if (ipiv_c == NULL) errMsg("malloc() allocation failure for ipiv_c!");
  for (i = 0; i < ipiv_dim; ++i) {
    ipiv_c[i] = 0;
  }

  F77_CALL(dgetrf)(&m_c, &n_c, a, &lda_c, ipiv_c, &info_c);

  *info = (double)info_c;

  for (i = 0; i < ipiv_dim; ++i) {
    ipiv[i] = (double)ipiv_c[i];
  }
}
/* ******************** */
