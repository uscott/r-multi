
#include "uri.h"


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  The following functions are for computing implied volatility.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/






static double bserr(const double T, const double S,
                    const double K, const double vol, const double r, 
                    const double q, const char   opt, const double mprice)

  /*******************************************************************
   *
   *  Description: Computes the BS options price x for the 
   *    given parameters T, S, K, vol, r, q, opt and returns
   *    log(x / mprice).
   *
   *  Caution: T is in years (or whatever units vol, r and q
   *    are in).
   *
   *******************************************************************/

{
  double bsval;

  bs0(T, S, K, vol, r, q, opt, 'p', &bsval);

  return log(bsval / mprice);

}





static void impvol0(const double T, const double S, const double K, 
                    const double mprice, const double r, const double q, 
                    const char   optype, const double tol, int maxit, 
                    double *iv)

  /*******************************************************************
   *
   *  Description: Computes the implied volatility and assigns it
   *    to *iv.
   *
   *  Caution: T should be in the same time units as r and q
   *    (usually years).
   *
   *******************************************************************/

{

  int it = 0;
  int ok;
  double lower, upper;
  double err, err_l, err_u;
  double errA, errB;
  double price_bound_lower = 0.0;
  double price_bound_upper = 0.0;

  iv ? 0.0 : error ("Null pointer passed where not allowed");

  /* Get lower and upper arbitrage bounds for the price */
  switch(optype) {

  case 'c':

    price_bound_lower = MAX(exp(-q * T) * S - exp(-r * T) * K, 0);
    price_bound_upper = S;
    break;

  case 'p':

    price_bound_lower = MAX(exp(-r * T) * K - exp(-q * T) * S, 0);
    price_bound_upper = K;
    break;

  case 's':
    
    price_bound_lower = 
      MAX(exp(-q * T) * S - exp(-r * T) * K, 0) +
      MAX(exp(-r * T) * K - exp(-q * T) * S, 0);
    price_bound_upper = S + K;
    break;

  default:
    break;

  }

  ok = price_bound_lower < mprice && mprice < price_bound_upper;

  if (!ok) {

    *iv = NA_REAL;
    return;

  }

  maxit = MAX(maxit, 20);
  lower = 1.0;
  it    = 0;
  
  do {

    lower *= 0.5;
    it    += 1;
    err_l  = bserr(T, S, K, lower, r, q, optype, mprice);

  } while (it < maxit && err_l > 0);

  it    = 0;
  upper = 0.5;

  do {

    upper *= 2.0;
    it    += 1;
    err_u  = bserr(T, S, K, upper, r, q, optype, mprice);

  } while (it < maxit && err_u < 0);

  
  if (0 == err_l)

    *iv = lower;

  else if (0 == err_u)

    *iv = upper;

  else if (err_l > 0 || err_u < 0)

    *iv = NA_REAL;

  else {

    it = 0;

    do {

      *iv  = 0.5 * (lower + upper);
      err  = bserr(T, S, K, *iv, r, q, optype, mprice);
      it  += 1;

      errA = bserr(T, S, K, upper, r, q, optype, mprice);
      errB = bserr(T, S, K, lower, r, q, optype, mprice);

      if (SIGN_OF(err) == SIGN_OF(errA))

        upper = *iv;

      else if (SIGN_OF(err) == SIGN_OF(errB))

        lower = *iv;

      else
        break;

    } while (ABS(err) > tol && it < maxit);

  }    
  
}







#define N_NUMERIC_ARGS 6
static void imp_vol_chkargs(SEXP *tau, SEXP *S, SEXP *K, SEXP *mprice, 
                            SEXP *r, SEXP *q, SEXP *op_type,
                            int  *numprot)
{

  int i, N = 0;

  SEXP *args[]  = {tau, S, K, mprice, r, q, op_type};
  char *names[] = {"tau", "S", "K", "mprice", "r", "q", "op_type"};
  
  /* Coerce types  */
  for (i = 0; i < N_NUMERIC_ARGS; i++)
    PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

  /* Find maximal length */
  for (i = 0, N = 0; i < N_NUMERIC_ARGS + 1; i++)
    N = MAX(N, length(*args[i]));

  if (N > 1) { /* Repeat arguments of length == 1. */

    for (i = 0; i < N_NUMERIC_ARGS; i++)
      if (1 == length(*args[i]))
        PROT2(*args[i] = numrep(*args[i], N), *numprot);

    if (1 == length(*op_type))
      PROT2(*op_type = charrep(*op_type, N), *numprot);

  }

  /* Check lengths */

  for (i = 0; i < N_NUMERIC_ARGS + 1; i++)
    if (N != length(*args[i]))
      error ("Argument %s of wrong length", names[i]);
  

}
#undef N_NUMERIC_ARGS



SEXP imp_vol(SEXP tau, SEXP S, SEXP K, 
             SEXP mprice, SEXP r, SEXP q,
             SEXP op_type, SEXP tol, SEXP maxit)

  /*******************************************************************
   *
   *  Description: Computes the implied BS volatility of an option.
   *    
   *  Caution: Assumes tau is in days.
   *
   *******************************************************************/

{  
  SEXP ans;
  
  int N, i;
  int numprot = 0;

  imp_vol_chkargs(&tau, &S, &K, &mprice, &r, &q, &op_type, &numprot);
  
  /* By now tau, S, K, mprice, r, q, op_type should all have same length */

  N = length(tau); 

  PROT2(ans = NEW_NUMERIC(N), numprot);

  for (i = 0; i < N; i++) {

    impvol0(REAL(tau)[i] / 365.0, 
            REAL(S)[i], REAL(K)[i], REAL(mprice)[i],
            REAL(r)[i], REAL(q)[i], *CHAR(GET_ELT(op_type, i)),
            asReal(tol), asInteger(maxit), REAL(ans) + i);

  }
   
  UNPROTECT(numprot);

  return ans;

}


