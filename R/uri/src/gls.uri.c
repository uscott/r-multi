
#include "uri.h"


void gls_sumsq1(int T, int n, int p, double *y, double **X,
                double *beta, double *phi, double *sumsq, double *res)

  /********************************************************************
  *
  *   Description: Computes the residuals and sum of squared
  *     residuals in the regression equation
  *     y = X * beta + e where the residuals e[t] follow an AR(p)
  *     process with coefficients phi[0], ..., phi[p - 1].
  *
  *   Inputs:
  *     T =
  *     number of observations y[0], ..., y[T - 1] and also the
  *     number of rows of the matrix (represented by) X.
  *
  *     n =
  *     number of columns of X.
  *
  *     p =
  *     order of AR process for residuals.
  *
  ********************************************************************/

{
  const double *phi0  = phi; /* For keeping track of memory starting location */
  const double *beta0 = beta;

  double *x, *res2;
  int t, j, lag;

  res2 = res;

  for (t = 0; t < T; t++, res2++) {

    beta  = (double *) beta0;
    x     = *X++;
    *res2 = *y++;
    
    for (j = 0; j < n; j++)
      *res2 -= (*beta++) * (*x++);

  }

  res    = res + T - 1;
  *sumsq = 0.0;


  /* Following loop replaces res[t] from regression equation y = X * beta 
     with residual from AR(p) process for t = p, ..., T - 1. */

  for (t = T - 1; t >= p; t--, res--) {

    res2 = res - 1;
    phi  = (double *) phi0;

    for (lag = 1; lag <= p; lag++)
      *res -= (*phi++) * (*res2--);

    *sumsq += SQR(*res);
  }

  for (t = p - 1; t > 0; t--, res--) {

    res2 = res - 1;
    phi  = (double *) phi0;

    for (lag = 1; lag < t; lag++)
      *res -= (*phi++) * (*res2--);


    *sumsq += SQR(*res);
  }

  *sumsq += SQR(res[0]);

}






SEXP gls_sumsq2(SEXP y, SEXP X, SEXP beta, SEXP phi)
{

  SEXP sumsq, ans, ans_names, res;
  carray X_array;
  int numprot = 0;
  const int T = length(y);
  const int n = ncols(X);
  const int p = length(beta);

  if (!isMatrix(X)) 
    error("Non-matrix where matrix expected in gls_sumsq2");

  if (length(beta) != ncols(X) || length(y) != nrows(X))
    error("Dimension mismatch in gls_sumsq2");

  if (length(y) <= length(phi))
    error("Too many AR(p) parameters passed to gls_sumsq2");


  PROTECT(y         = AS_NUMERIC(y));
  PROTECT(sumsq     = NEW_NUMERIC(1));
  PROTECT(res       = NEW_NUMERIC(T));
  PROTECT(ans       = NEW_LIST(2));
  PROTECT(ans_names = NEW_STRING(2));
  numprot += 5;

  X_array = sexp_to_carray(X, /* dup = */ 1);

  gls_sumsq1(T, n, p, 
             REAL(y), 
             ARRAY2(X_array), 
             REAL(beta),
	     REAL(phi), 
             REAL(sumsq), 
             REAL(res));

  SET_ELT(ans, 0, sumsq);
  SET_ELT(ans, 1, res);

  CHARACTER_POINTER(ans_names)[0] = mkChar("sum.of.squares");
  CHARACTER_POINTER(ans_names)[1] = mkChar("residuals");

  SET_NAMES(ans, ans_names);

  UNPROTECT(numprot);

  return ans;
}






