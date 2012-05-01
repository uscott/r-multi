
#include "uri.h"

/*
** General functions related to copulas.
*/


void empcdf(const double *x,    const double *y, 
            const int    *xlen, const int    *ylen, double *ans)

  /********************************************************************
   *
   *   Description: Evaluates the empirical cdf of x[0],..., x[xlen - 1]
   *    on y[0],..., y[ylen - 1] storing the values in addresses
   *    ans,..., ans + ylen - 1.  For this function to work properly
   *    none of the x[i] should be NA or NaN.
   *
   ********************************************************************/

{
  double const *xB = x;
  double const *xE = x + *xlen;
  double const *yE = y + *ylen;

  for ( ; y < yE; y++, ans++) {

    if (ISNAN(*y))
      *ans = NA_REAL;
    else {
      for (*ans = 0.0, x = xB; x < xE; x++)
        if (*x < *y)
          *ans += 1.0;

      *ans += 0.5;
      *ans /= (double) *xlen;
    }

  }
}




SEXP empcdf2(SEXP x, const SEXP y)

  /********************************************************************
   *
   *   Description: Evaluates the empirical cdf of 
   *    REAL(x)[0],..., REAL(x)[length(x) - 1]
   *    on REAL(y)[0],..., REAL(y)[length(y) - 1].
   *    This function removes NA's and NaN's from x.
   *
   ********************************************************************/

{
  int xlen, ylen;
  SEXP ans;

  ylen = length(y);

  PROTECT(ans = NEW_NUMERIC(ylen));
  PROTECT(x   = na_omit(x));

  xlen = length(x);

  empcdf(REAL(x), REAL(y), &xlen, &ylen, REAL(ans));

  UNPROTECT(2);

  return ans;

}




void empunif(const double *x, const int *xlen, double *ans)

  /********************************************************************
   *
   *   Description: Stores the values of the empirical uniform
   *    random variates corresponding to x[0],...,x[xlen - 1] in
   *    ans,...,ans + xlen - 1.  None of the x[i] should be
   *    NA or NaN.
   *
   ********************************************************************/

{
  empcdf(x, x, xlen, xlen, ans);
}




void empnorm(const double *x, const int *xlen, double *ans)

  /*******************************************************************
   *
   *  Description: Stores the values of the empirical normal
   *     random variates corresponding to x[0],...,x[xlen - 1] in
   *     ans,...,ans + xlen - 1.
   *
   *******************************************************************/
{
  double const *fin = ans + *xlen;

  empcdf(x, x, xlen, xlen, ans);

  for ( ; ans < fin; ans++)
    *ans = qnorm(*ans, 0.0, 1.0, 1, 0);
}





SEXP empunif2(SEXP X)
{
  int  j, n, numprot = 0;
  SEXP ans = R_NilValue, tmp = R_NilValue;

  if (isMatrix(X)) {

    PROT2(ans = NEW_LIST(n = ncols(X)), numprot);

    for (j = 0; j < n; j++) {

      PROTECT(tmp = matcol2(X, j));
      SET_ELT(ans, j, empcdf2(tmp, tmp));
      UNPROTECT(1);

    }

  }
  else if (isNewList(X)) {

    PROT2(ans = NEW_LIST(n = length(X)), numprot);

    for (j = 0; j < n; j++)
      SET_ELT(ans, j, empcdf2(GET_ELT(X, j), GET_ELT(X, j)));

  }
  else if (isNumeric(X)) {
    
    PROT2(ans = empcdf2(X, X), numprot);

  }

  UNPROTECT(numprot);

  return ans;
}




SEXP empnorm2(SEXP X)
{
  int numprot = 0;
  long i, j, m;
  register double *p;
  SEXP ans = R_NilValue;

  PROT2(ans = empunif2(X), numprot);

  if (isNewList(ans))
    for (j = 0; j < length(ans); j++) {
      
      p = REAL(GET_ELT(ans, j));
      m = length(GET_ELT(ans, j));
      
      for (i = 0; i < m; i++, p++)
        *p = qnorm(*p, 0.0, 1.0, 1, 0);

    }
  else if (isNumeric(ans)) {
    
    p = REAL(ans);
    m = length(ans);

    for (i = 0; i < m; i++, p++)
      *p = qnorm(*p, 0.0, 1.0, 1, 0);
  }
  else
    ans = R_NilValue;

  UNPROTECT(numprot);

  return ans;
}


/*
** Functions related to Plackett copula follow.
*/



void C_plackett(const double *u,   const double *v, 
                const double *phi, const int    *n, double *ans)
{
  const double    p      = *phi - 1.0;
  register double t1, t2;
  double const   *ansEND = ans + *n;

  if (ABS(p) < DOUBLE_EPS)
    for ( ; ans < ansEND; ans++, u++, v++)
      *ans = (*u) * (*v);
  else
    for ( ; ans < ansEND; ans++, u++, v++) {

      t1 = 1.0 + p * (*u + *v);
      t2 = sqrt(SQR(1.0 + p * (*u + *v)) - 4.0 * (*u) * (*v) * (*phi) * p);
      *ans = 0.5 / p * (t1 - t2);

    }

}





void ll_plackett(const double *u,   const double *v,
                 const double *phi, const int    *n, double *ll)
{
  const double    p  = *phi - 1.0;
  const double  lphi = log(*phi);
  double const   *uE = u + *n;
  register double a, b;

  *ll = 0.0;

  for ( ; u < uE && R_FINITE(*ll); u++, v++) {
    a = lphi + log(1.0 + (*u - 2.0 * (*u) * (*v) + *v) * p);
    b = 1.5 * log(SQR(1.0 + p * (*u + *v)) - 4.0 * (*u) * (*v) * (*phi) * p);
    *ll += a - b;
  }

  if (!R_FINITE(*ll))
    *ll = R_NaN;
}




void sim_plackett_1(double       *u,   double    *v, 
                    const double *phi, const int *n)

     /***************************************************************
      *
      *   Description: Simulates n draws from the bivariate Plackett
      *     copula and puts the result in u[0], ..., u[*n - 1],
      *    v[0], ..., v[*n - 1].
      *
      ***************************************************************/

{
  const    double a0 = -(*phi + 1.0) / (*phi - 1.0);
  const    double p  = *phi - 1.0;
  register double d, e, a1, b1, b2, c0, c1, c2, t;
  double const   *uE = u + *n;

  if (*n < 1)
    error ("Negative number of draws requested in sim_plackett_1");
  if (*phi < 0)
    error ("Negative phi passed to sim_plackett_1");

  GetRNGstate();

  for ( ; u < uE; u++, v++) {
    
    *u = unif_rand();
     t = unif_rand();

    d = SQR(1.0 - 2.0 * t);
    e = t < 0.5 ? -1.0 : 1.0;

    a1 =  2.0 * (*phi) / p + 2.0 * (*phi) * (*u);
    b1 = -4.0 * (*u) * (*phi);
    b2 =  4.0 * (*u) * (*phi) + 4 * SQR(*u) * (*phi) * p;
  
    c0 = SQR(a0) - d;
    c1 = 2.0 * a0 * a1 - d * b1;
    c2 = SQR(a1) - d * b2;


    *v = (-c1 + e * sqrt(c1 * c1 - 4.0 * c0 * c2)) / (2 * c0);
    *v = (*v - 1.0) / p - *u;
  }
  
  PutRNGstate();

}


/*
** Gaussian copula stuff.
*/



SEXP ll_gaussian_mvt(const SEXP x, const SEXP rho, const SEXP transform)
{

  const long T = nrows(x);
  const int  n = ncols(x);
  int  numprot = 0;
  register int      i, j;
  register long     t;
  register double   ll;
  register double **xp = NULL;
  carray             r, a;
  SEXP              xn = R_NilValue;

  if (!isMatrix(x) || !isMatrix(rho) || n != nrows(rho) || n != ncols(rho))
    error ("C error: ll_gaussian_mvt.  Bad arguments.");

  ll = R_NaN;

  if (2 <= n) {

    r = sexp_to_carray(rho, 1);
    a = make_zero_matrix(n, n);

    if (matrix_inverse(r, a)) {

      array_op(a, make_identity_matrix(n), '-', a);
    
      if (asInteger(transform)) {

        PROTECT(xn = empnorm2(x));
        numprot   += 1;
        xp         = ARRAY2(sexp_to_carray(xn, 1));

      }
      else
        xp = ARRAY2(sexp_to_carray(x, 1));

      ll = T * log(det(r));

      for (t = 0; t < T && R_FINITE(ll); t++) {
        for (i = 0; i < n; i++)
          for (j = 0; j < n; j++)
            ll += xp[t][i] * ARRAY2(a)[i][j] * xp[t][j];
      }

      ll *= -0.5;

    } 

  }

  UNPROTECT(numprot);
  
  return ScalarReal(R_FINITE(ll) ? ll : R_NaN);
}




SEXP ll_gaussian_bvt(const SEXP x,   const SEXP y, 
                     const SEXP rho, const SEXP transform)
{

  const    long    T       = length(x);
  register long    t;
  int              numprot = 0;
  register double  ll      = 0.0,         r = 0.0;
  register double *xp      = NULL,      *yp = NULL;
  SEXP             xn      = R_NilValue, yn = R_NilValue;

  r = asReal(rho);

  if (ABS(r) <= 1.0) {
    
    if (asInteger(transform)) {

      PROTECT(xn = empnorm2(AS_NUMERIC(x)));
      PROTECT(yn = empnorm2(AS_NUMERIC(y)));
      numprot   += 2;

    }
    else {

      PROTECT(xn = AS_NUMERIC(x));
      PROTECT(yn = AS_NUMERIC(y));
      numprot   += 2;

    }

    if (T != length(xn) || T != length(yn))
      error ("Bad arguments.");

    xp = REAL(xn);
    yp = REAL(yn);

    for (t = 0, ll = 0.0; t < T && R_FINITE(ll); t++, xp++, yp++)
      ll += r * SQR(*xp) - 2.0 * (*xp) * (*yp) + r * SQR(*yp);

    ll = -0.5 * T * log(1.0 - SQR(r)) - 0.5 * r / (1.0 - SQR(r)) * ll;

  } 
  else 
    ll = R_NaN;

  UNPROTECT(numprot);
  
  return ScalarReal(ll);
}




/*
** Gumbel copula stuff below.
*/




void gumbel_ml(int *T, double *u, double *v, double *alpha, double *loglik)
{
  int t;
  double uu, vv, z;

  const double a = *alpha;

  *loglik = 0;

  for (t = 0; t < (*T); t++) {

    uu = -log(*u);
    vv = -log(*v);
    z = R_pow(uu, a) + R_pow(vv, a);

    *loglik += (a - 1) * log(uu * vv) - log((*u++) * (*v++))
      + (1.0/a - 2.0) * log(z) - R_pow(z, 1.0/a)
      + log(R_pow(z, 1.0/a) + a - 1.0);
  }
}











