#include "uri.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This file contains functions for calculating the BS formula, greeks,
**  implied vol and other stuff.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/




void bs0(const double T, const double S, const double K, 
         const double vol, const double r, const double q,   
         char opt, char ret, double *ans)

  /*******************************************************************
   *
   *  Description: Computes the BS price, delta or vega and sets *ans to it.
   *
   *  Caution: Assumes T is given in *years* (or whatever time
   *    unit vol, r and q are given in).
   *
   *  Last modified: July 28 2003
   *
   *******************************************************************/

{
  int ok;
  double callval, a, b;
  double d1, d2;

  a = 1.0 + ('s' == opt);
  d1 = (log(S/K)+(r-q+0.5*DSQR(vol))*T) / (vol*sqrt(T));
  d2 = d1 - vol*sqrt(T);

  opt = tolower(opt);
  ret = tolower(ret);

  if (!ans)
    error ("Illegal null pointer passed");
  
  *ans = NA_REAL;

  ok = S >= 0 && K >= 0 && R_FINITE(q) &&
    vol >= 0 && R_FINITE(T) && R_FINITE(S) && 
    R_FINITE(K) && R_FINITE(vol) && R_FINITE(r) && T >= 0.0;
    
  if (!ok)
    return;

  switch(ret) {

  case 'p':
    callval = (0.0 == T) ? DMAX(S-K,0) : exp(-q*T)*S*NORM_DIST(d1)-exp(-r*T)*K*NORM_DIST(d2);
    b = ('p' == opt) ? -exp(-q*T)*S + exp(-r*T)*K : 0.0;
    break;

  case 'd':
    callval = (0.0 == T) ? (S > K ? 1.0 : 0.0) : exp(-q*T)*NORM_DIST(d1);
    b = ('p' == opt) ? -exp(-q*T) : 0.0;
    break;

  case 'v':
    callval = (0.0 == T) ? 0.0 : exp(-q*T)*S*sqrt(T)*NORM_DENS(d1);
    b = 0.0;
    break;

  case 'g':
    callval = (0.0==T)?(exp(-q*T)*S==exp(-r*T)*K?R_PosInf:0.0):exp(-q*T)*NORM_DENS(d1)/(S*vol*sqrt(T));
    b = 0.0;
    break;

  default:
    callval = b = NA_REAL;
    break;

  }

  *ans = a*callval + b;

}



#define N_NUMARGS 6
static void bs_chkargs(SEXP *tau, SEXP *S, SEXP *K, 
                       SEXP *vol, SEXP *r, SEXP *q,
                       SEXP *opt, SEXP *ret, int  *numprot)
{
  int i, N;

  SEXP *args[N_NUMARGS+1]  = {tau, S, K, vol, r, q, opt};
  char *names[N_NUMARGS+1] = {"tau", "S", "K", "vol", "r", "q", "opt"};

  /* Coerce types  */
  for (i = 0; i < N_NUMARGS; i++)
    ENSURE_NUMERIC(*args[i], *numprot);

  ENSURE_CHAR(*opt, *numprot);

  if (length(*ret) && !isNull(*ret))
    ENSURE_CHAR(*ret, *numprot);

  /* Find maximal length */
  for (i = N = 0; i < N_NUMARGS + 1; i++)
    N = MAX(N, length(*args[i]));

  if (N > 1) {/* Repeat arguments of length == 1. */

    for (i = 0; i < N_NUMARGS; i++)
      if (1 == length(*args[i]))
        PROT2(*args[i] = numrep(*args[i], N), *numprot);

    if (1 == length(*opt))
      PROT2(*opt = charrep(*opt, N), *numprot);

  }

  /* Check lengths */
  for (i = 0; i < N_NUMARGS + 1; i++)
    if (N != length(*args[i]))
      error ("Argument %s of wrong length", names[i]);

}
#undef N_NUMARGS





SEXP bs(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP r, SEXP q, SEXP opt, SEXP ret)

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

  bs_chkargs(&tau, &S, &K, &vol, &r, &q, &opt, &ret, &numprot);

  m = length(tau); /* By now all arguments should have same length */
  retlen = length(ret);

  ans_matrix = make_zero_matrix(m, retlen);
  
  for (j = 0; j < retlen; j++)
    for (i = 0; i < m; i++)
      bs0(REAL(tau)[i]/365.0, REAL(S)[i], REAL(K)[i],
          REAL(vol)[i], REAL(r)[i], REAL(q)[i],
          *CHAR(GET_ELT(opt, i)), *CHAR(GET_ELT(ret, j)),
          ARRAY2(ans_matrix)[i] + j);

  PROT2(ans = carray_to_sexp(ans_matrix), numprot);
  PROT2(dimnames = NEW_LIST(2), numprot);
  
  SET_ELT(dimnames, 0, GET_NAMES(tau));
  SET_ELT(dimnames, 1, ret);

  setAttrib(ans, R_DimNamesSymbol, dimnames);

  UNPROTECT(numprot);

  return ans;
}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Functions for getting the strike corresponding to a certain delta.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/




static void dtoK(const double T, const double S, 
                 const double vol, const double del, 
                 const char  opt, const double r,   
                 const double q, double *ans)

  /*******************************************************************
   *
   *  Description: 
   *
   *  Caution: Assumes T is in years.
   *
   *******************************************************************/

{

  register double t1, t2;

  if (T < 0 || S < 0 || vol < 0) {

    *ans = NA_REAL; 
    return;

  }

  switch(opt) {

    case 'c':

      t1   = -vol * sqrt(T) * qnorm(exp(q * T) * del,0, 1, 1, 0);
      t2   = (r - q + 0.5 * SQR(vol)) * T;
      *ans = S * exp(t1 + t2);
      break;

    case 'p':

      t1   = -vol * sqrt(T) * qnorm(1.0 + exp(q * T) * del, 0, 1, 1, 0);
      t2   = (r - q + 0.5 * SQR(vol)) * T;
      *ans = S * exp(t1 + t2);
      break;

    case 's':

      t1   = -vol * sqrt(T) * qnorm(0.5 + 0.5 * exp(q * T) * del, 0, 1, 1, 0);
      t2   = (r - q + 0.5 * SQR(vol)) * T;
      *ans = S * exp(t1 + t2);
      break;

    default:

      *ans = NA_REAL;
      break;

  }

}




#define NARGS 7
static void deltaToStrike_chk(SEXP *T, SEXP *S, SEXP *vol, 
                              SEXP *del, SEXP *optionType,
                              SEXP *r, SEXP *q, int *numprot)
{
  SEXP *args[NARGS] = {T, S, vol, del, optionType, r, q};
  const SEXPTYPE type[NARGS] = {REALSXP, REALSXP, REALSXP,
                                REALSXP, STRSXP,
                                REALSXP, REALSXP};
  long lens[NARGS] = {0};
  long maxLen;
  int i;
  
  for (i = 0; i < NARGS; i++) {

    PROT2(*args[i] = coerceVector(*args[i], type[i]), *numprot);
    lens[i] = length(*args[i]);

  }

  maxLen = lmax(lens, NARGS);

  if (maxLen > 1)
    for (i = 0; i < NARGS; i++) {
      if (1 == lens[i])
        PROT2(*args[i] = rep(*args[i], type[i], maxLen), *numprot);
      else if (lens[i] != maxLen)
        error ("mismatched arg lengths");
    }
  
}



SEXP deltaToStrike(SEXP T, SEXP S, SEXP vol, 
                   SEXP del, SEXP optionType,
                   SEXP r, SEXP q)
{
  int numprot = 0;
  long len, i;
  SEXP K;

  deltaToStrike_chk(&T, &S, &vol, &del, &optionType, &r, &q, &numprot);

  len = length(T);
  PROT2(K = NEW_NUMERIC(len), numprot);

  for (i = 0; i < len; i++)    
    dtoK(REAL(T)[i] / 365.0, REAL(S)[i], REAL(vol)[i],
         REAL(del)[i], *CHAR(GET_ELT(optionType, i)),
         REAL(r)[i], REAL(q)[i],
         REAL(K) + i);

  UNPROTECT(numprot);

  return K;

}









