#include "uri.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Digital option stuff.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


void digital0(const double T, const double S, const double K, 
              const double vol, const double r, const double q,   
              char ret, double *ans)

  /*******************************************************************
   *
   *  Description: Computes the digital option price, delta or vega and sets *ans to it.
   *
   *  Caution: Assumes T is given in *years* (or whatever time
   *    unit vol, r and q are given in).
   *
   *  Last modified: July 28 2003
   *
   *******************************************************************/

{
  int ok;
  double d1, d2, tmp;

  d1 = (log(S/K)+(r-q+0.5*DSQR(vol))*T) / (vol*sqrt(T));
  d2 = d1 - vol*sqrt(T);
  ret = tolower(ret);

  if (!ans)
    error ("null pointer passed");
  
  *ans = NA_REAL;

  ok = S >= 0 && K >= 0 && R_FINITE(q) &&
    vol >= 0 && R_FINITE(T) && R_FINITE(S) && 
    R_FINITE(K) && R_FINITE(vol) && R_FINITE(r) && T >= 0.0;
    
  if (!ok)
    return;

  switch(ret) {

  case 'p':
    *ans = exp(-r*T) * NORM_DIST(d2);
    break;

  case 'd':
    tmp = 1.0 / (S*vol*sqrt(T));
    *ans = tmp * exp(-r*T) * NORM_DENS(d2);
    break;

  case 'v':
    tmp = -(log(S/K)+(r-q)*T)/(vol*vol*sqrt(T)) - 0.5*sqrt(T);
    *ans = tmp * exp(-r*T)*NORM_DENS(d2);
    break;

  default:
    break;

  }

}




#define N_NUMARGS 6
static void digital_chkargs(SEXP *tau, SEXP *S, SEXP *K, 
                            SEXP *vol, SEXP *r, SEXP *q,
                            SEXP *ret, int  *numprot)
{
  int i, N;

  SEXP *args[N_NUMARGS]  = {tau, S, K, vol, r, q};
  char *names[N_NUMARGS] = {"tau", "S", "K", "vol", "r", "q"};

  /* Coerce types  */
  for (i = 0; i < N_NUMARGS; i++)
    ENSURE_NUMERIC(*args[i], *numprot);

  if (length(*ret) && !isNull(*ret))
    ENSURE_CHAR(*ret, *numprot);

  /* Find maximal length */
  for (i = N = 0; i < N_NUMARGS; i++)
    N = MAX(N, length(*args[i]));

  if (N > 1) /* Repeat arguments of length == 1. */
    for (i = 0; i < N_NUMARGS; i++)
      if (1 == length(*args[i]))
        PROT2(*args[i] = numrep(*args[i], N), *numprot);

  /* Check lengths */
  for (i = 0; i < N_NUMARGS; i++)
    if (N != length(*args[i]))
      error ("Argument %s of wrong length", names[i]);

}
#undef N_NUMARGS



SEXP digital(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP r, SEXP q, SEXP ret)

  /*******************************************************************
   *
   *  Description: Returns SEXP of type matrix with columns
   *    corresponding to the BS price, delta & vega.
   *
   *  Caution: Assumes tau is in *calender days*.
   *
   *******************************************************************/

{
  int numprot = 0;
  long i, j, m, retlen = 0;
  SEXP ans, dimnames;
  carray ans_matrix;

  digital_chkargs(&tau, &S, &K, &vol, &r, &q, &ret, &numprot);

  m = length(tau); /* By now all arguments should have same length */
  retlen = length(ret);

  ans_matrix = make_zero_matrix(m, retlen);
  
  for (j = 0; j < retlen; j++)
    for (i = 0; i < m; i++)
      digital0(REAL(tau)[i]/365.0, REAL(S)[i], REAL(K)[i],
               REAL(vol)[i], REAL(r)[i], REAL(q)[i],
               *CHAR(GET_ELT(ret, j)),
               ARRAY2(ans_matrix)[i] + j);

  PROT2(ans = carray_to_sexp(ans_matrix), numprot);
  PROT2(dimnames = NEW_LIST(2), numprot);
  
  SET_ELT(dimnames, 0, GET_NAMES(tau));
  SET_ELT(dimnames, 1, ret);

  setAttrib(ans, R_DimNamesSymbol, dimnames);

  UNPROTECT(numprot);

  return ans;
}


