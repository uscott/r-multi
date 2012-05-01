
#include "uri.h"

#include <stdlib.h>

#define MAX_NUM_OPTNS 100
#define THREE 3




/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  The following section is for computing the P&L and transaction
**  costs for dynamic discrete delta hedging using the BS deltas.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


static void bsd1_chk(const int n, const double *T, const double *S, 
                     const double  K, const double *vol, const double r, 
                     const double q, const int hedgePeriod, 
                     const char optionType, const double *del)

  /*******************************************************************
   *
   *  Description: Calls R error or warning function if arguments are 
   *    not as expected and does nothing otherwise.
   *
   *******************************************************************/

{
  int ok;

  ok = T && S && K && vol && del;

  if (!ok)
    error ("null pointer passed to bsd1");

  ok = R_FINITE((double)n) && R_FINITE(K) && R_FINITE((double)hedgePeriod);
    
  if (!ok)  
    warning("Improper n, K or hedgePeriod");

  ok = R_FINITE(r) && R_FINITE(q);

  if (!ok)
    warning ("Improper r or q");

  ok = n > 0  || *T >= 0 || K >= 0 || *vol >= 0;

  if (!ok)
    warning("n <= 0 or *T < 0 or K < 0 or *vol < 0");

  ok = ('c' == optionType || 'p' == optionType || 's' == optionType);
  
  if (!ok)
    warning ("optionType not one of c, p, s");
}








static void bsd1(const int Tlen, const double *T, 
                 const double *S, const double  K, 
                 const double *vol, const double r, 
                 const double q, const int hedgePeriod, 
                 const char optionType, double *del)

  /*******************************************************************
   *
   *  Description: Evaluates the deltas from the basic BS
   *    dynamic hedging strategy applied to a single 
   *    path S[0], ..., S[Tlen] of the underlying measured at days
   *    T[0], ..., T[Tlen - 1] from expiration.  The deltas are rebalanced
   *    every (hedgePeriod) days.  Note that vol[i] is the vol for
   *    computing the BS delta on day T[i] from expiration if a 
   *    rebalancing occurs on this day.
   *    Deltas are put in del[0], ..., del[Tlen].  Note that del[Tlen]
   *    will be 0, 1 or -1 corresponding to whether the option expires ITM.
   *
   *  Parameters:
   *
   *    n - Gives the length of the contiguous memory block
   *    pointed to by T.
   *
   *    T - T[0], ..., T[n - 1] are the times to expiration
   *    at which the prices S[0], ..., S[n - 1] are measured.
   *
   *    S - S[0], ..., S[n] constitute the price path with
   *    S[n] being the price at expiration and S[0], ..., S[n - 1]
   *    being measured at T[0], ..., T[n - 1].
   *
   *    vol - vol[0], ..., vol[n - 1] are the vols plugged into
   *    the BS function to get the deltas.
   *
   *    r, q - Interest rate and dividend
   *    yield per annum.
   *
   *    optionType - Should be one of 'c', 'p' or 's'.
   *
   *  Return Value:
   *    Overwrite del[0], ..., del[n] to be the deltas T[0], ...,
   *    T[n - 1] days from expiration and at expiration
   *    respectively.
   *    
   *  
   *  Status: Tested & OK.  Avoid making changes.
   *
   *******************************************************************/

{
  register int i, ok;
  register double d1, d2, tau, x, Nd1, divDiscFac, intDiscFac, ctr;

  bsd1_chk(Tlen, T, S, K, vol, r, q, hedgePeriod, optionType, del);

  /* First compute del[0] & del[Tlen] */

  tau = *T / 365.0;
  x   = log(*S / K);

  d1  = (x + (r - q + 0.5 * SQR(*vol)) * tau) / (*vol * sqrt(tau));
  d2  = d1 - *vol * sqrt(tau);

  Nd1  = pnorm(d1, 0, 1, 1, 0);
  divDiscFac  = exp(-q * tau);
  intDiscFac  = exp(-r * tau);

  if ('c' == optionType) {

    del[0]    = divDiscFac * Nd1;
    del[Tlen] = S[Tlen] > K ? 1.0 : 0.0;

  }
  else if ('p' == optionType) {

    del[0]    = divDiscFac * (Nd1 - 1.0);
    del[Tlen] = K > S[Tlen] ? -1.0 : 0.0;

  }
  else if ('s' == optionType) {

    del[0]    = divDiscFac * (2.0 * Nd1 - 1.0);
    del[Tlen] = SIGN_OF(S[Tlen] - K);

  }
  else
    del[0] = del[Tlen] = NA_REAL;

  /* Get rest of deltas */
  for (i = 1, ctr = 0; i < Tlen; i++) {

    ok = R_FINITE(T[i]) && R_FINITE(S[i]) && 
         T[i] >= 0      && vol[i] >= 0    && T[i - 1] > T[i];
    
    if (!ok)
      break;

    ctr += T[i - 1] - T[i];

    if (ctr < hedgePeriod)

      del[i] = del[i - 1];

    else {

      ctr = 0;

      tau = T[i] / 365.0;
      x   = log(S[i] / K);

      d1  = x + (r - q + 0.5 * SQR(vol[i])) * tau;
      d1 /= vol[i] * sqrt(tau);

      d2  = d1 - vol[i] * sqrt(tau);

      Nd1 = pnorm(d1, 0, 1, 1, 0);

      divDiscFac  = exp(-q * tau);

      switch(optionType) {

      case 'c':
        del[i] = divDiscFac *  Nd1;
        break;
      case 'p':
        del[i] = divDiscFac * (Nd1 - 1.0);
        break;
      case 's':
        del[i] = divDiscFac * (2.0 * Nd1 - 1.0);
        break;
      default:
        del[i] = NA_REAL;
        break;

      }
      
    }


  }


}













static int bsdh1(const long Tlen, const int numStrikes,
                 const double *T, const double *S, 
                 const double *K, const double *vol,
                 const double *posn, const double r, 
                 const double q, const int hedgePeriod,   
                 char *optionTypes[], const double *transCosts, 
                 const int reltc, double *deltas,
                 double *PL, double *dly_PL, 
                 double *totalTransCosts)
     
  /*******************************************************************
   *
   *  Description: Computes the deltas and the corresponding P&L
   *    for delta hedgas done from time T[0] to T[Tlen - 1] _DAYS_
   *    away from expiration rebalanced every *hedgePeriod days.
   *    The position being hedged is an option position
   *    with strikes K[0], ..., K[numStrikes - 1], positions
   *    posn[0], ..., posn[numStrikes - 1], and option types
   *    *optionTypes[0], ..., *optionTypes[numStrikes - 1] 
   *    (each of which should
   *    be one of the 3 characters 'c', 'p', 's').
   *    Note there should be *Tlen + 1 values S[0], ...., S[Tlen]
   *    with S[Tlen] corresponding to the price at expiration.
   *    vol[i] is the vol used for computing the BS delta T[i]
   *    days from expiration.  
   *    The values transCosts[0], ..., transCosts[Tlen]
   *    provide transaction costs per contract 
   *    which are relative if *reltc is
   *    non-zero and absolute otherwise.  The P&L (not including
   *    transaction costs) is put into PL[0].
   *    Parameters r, q stand for interest rate and
   *    dividend yield respectively.  Total daily realized 
   *    discounted transaction costs are
   *    kept track of seperately in *totalTransCosts.
   *
   *    Status: Tested & OK.
   *
   *******************************************************************/

{
  int i, j, ok;
  register double dS, div, discFac, dt, transCostToday;
  double *tmpDeltas;

  ok = 
    Tlen > 0 && T && S && K && 
    vol && posn  && numStrikes > 0 && 
    optionTypes && deltas && 
    PL && transCosts && totalTransCosts && 
    dly_PL;

  if (!ok)
    error ("null pointer passed or negative lengths");

  tmpDeltas = (double *) R_alloc(Tlen + 1, sizeof(double));

  /* Compute all the deltas */

  dsetzero(tmpDeltas, Tlen + 1);
  dsetzero(deltas,    Tlen + 1);

  for (j = 0; j < numStrikes; j++) {

    bsd1(Tlen, T, S, *K++, 
         vol, r, q, hedgePeriod, 
         *optionTypes[j], tmpDeltas);

    for (i = 0; i <= Tlen; i++)
      deltas[i] += posn[j] * tmpDeltas[i];

  }

  discFac = 1.0; /* discFac = discount factor */

  *PL = *dly_PL = 0.0;

  *totalTransCosts = transCosts[0] * ABS(deltas[0]) * (reltc ? S[0] : 1.0);

  for (i = 1; i <= Tlen; i++) {

    dt  = (i < Tlen ? T[i - 1] - T[i]: T[Tlen - 1]) / 365.0;

    discFac *= exp(-r * dt);

    dS  = S[i] - S[i - 1];

    div = S[i - 1] * expm1((q - r) * dt);

    transCostToday  = transCosts[i] * ABS(deltas[i] - deltas[i - 1]);
    transCostToday *= discFac;

    reltc ? transCostToday *= S[i] : 0;

    *totalTransCosts += transCostToday;

    dly_PL[i] = -discFac * (deltas[i - 1] * (dS + div));
    
    *PL += dly_PL[i];

  }

  return 0;

}










#define NARGS 11
static int bsDeltaHedge1_chk(SEXP *T, SEXP *S, SEXP *K,
                             SEXP *vol, SEXP *posn, SEXP *r,
                             SEXP *q, SEXP *hedgePeriod,
                             SEXP *optionTypes, SEXP *transCosts,
                             SEXP *reltc, int *numprot)
{
  long T_len, S_len;
  int i, numStrikes;

  SEXP *args[NARGS] = {T, S, K, vol, posn, r, q, hedgePeriod,
                       optionTypes, transCosts, reltc};

  for (i = 0; i < 8; i++)
    if ( !IS_NUMERIC(*args[i]) )
      PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

  if (!IS_CHARACTER(*optionTypes))
    PROT2(*optionTypes = coerceVector(*optionTypes, STRSXP), *numprot);

  if (!IS_NUMERIC(*transCosts))
    PROT2(*transCosts = AS_NUMERIC(*transCosts), *numprot);

  T_len      = length(*T);
  S_len      = length(*S);
  numStrikes = length(*K);


  if ( S_len != T_len && S_len != T_len + 1)
    error ("length mismatch between S and T series");

  if (1 == length(*vol))
    PROT2(*vol = numrep(*vol, T_len), *numprot);
  
  if (1 < length(*r))
    warning ("only 1st value of r used");

  if (1 < length(*q))
    warning ("only 1st value of q used");

  if (1 == length(*transCosts))
    PROT2(*transCosts = numrep(*transCosts, S_len), *numprot);

  if (T_len != length(*vol))
    error ("vol wrong length");
  
  if (S_len != length(*transCosts))
    error ("transCosts wrong length");

  if (length(*posn) != numStrikes)
    error ("length mismatch between posn and K");

  if (length(*optionTypes) != numStrikes)
    error ("length mismatch between optionTypes and K");

  return 0;

}
#undef NARGS







#define ANS_LEN 4
SEXP bsDeltaHedge1(SEXP T, SEXP S, SEXP K, 
                   SEXP vol, SEXP posn, SEXP r,
                   SEXP q, SEXP hedgePeriod,
                   SEXP optionTypes, SEXP transCosts,
                   SEXP reltc)
{
  int numprot = 0;

  long j;
  long T_len, S_len;
  long numStrikes;
  const int hPeriod = asInteger(hedgePeriod);
  char *names[ANS_LEN] = {"deltas", "PL", "dailyPL", "totalTransCosts"};
  char *opt[MAX_NUM_OPTNS];

  double PL, totalTransCosts;
  
  SEXP ans, deltas, dly_PL;



  bsDeltaHedge1_chk(&T, &S, &K, 
                    &vol, &posn, &r, 
                    &q, &hedgePeriod,
                    &optionTypes, &transCosts, 
                    &reltc, &numprot);


  T_len      = length(T);
  S_len      = length(S);
  numStrikes = length(K);

  if (numStrikes > MAX_NUM_OPTNS)
    error ("too many strikes!");

  for (j = 0; j < numStrikes; j++) {

     opt[j] = (char *) R_alloc(1, sizeof(char));
    *opt[j] = *CHAR(STRING_ELT(optionTypes, j));

  }



  PROT2(deltas = NEW_NUMERIC(S_len), numprot);
  PROT2(dly_PL = NEW_NUMERIC(S_len), numprot);
  PROT2(ans    = NEW_LIST(ANS_LEN),  numprot);

  SET_NAMES(deltas, GET_NAMES(S));
  SET_NAMES(dly_PL, GET_NAMES(S));

  bsdh1(T_len, numStrikes, REAL(T), REAL(S),
        REAL(K), REAL(vol), REAL(posn), *REAL(r),
        *REAL(q), hPeriod, opt, REAL(transCosts),
        asInteger(reltc), REAL(deltas), &PL,
        REAL(dly_PL), &totalTransCosts);


  SET_ELT(ans, 0, deltas);
  SET_ELT(ans, 1, ScalarReal(PL));
  SET_ELT(ans, 2, dly_PL);
  SET_ELT(ans, 3, ScalarReal(totalTransCosts));

  set_names(ans, names);
  

  UNPROT2;

  return ans;
}
#undef ANS_LEN







#define NARGS 12
static void bsDeltaHedge2_chk(SEXP *tau, SEXP *S, SEXP *K, 
                              SEXP *vol, SEXP *posn, SEXP *r,
                              SEXP *q, SEXP *hedgePeriod, SEXP *optionTypes,
                              SEXP *transCosts, SEXP *reltc, 
                              SEXP *returnDaily,
                              int  *numprot)
{
  long S_len_m1;
  int numStrikes, i, j;
  SEXP tmp;
  SEXP *args[NARGS] = {tau, S, K, vol, 
                       posn, r, q, hedgePeriod,
                       optionTypes, transCosts, 
                       reltc, returnDaily};

  const SEXPTYPE type[NARGS] = {REALSXP, REALSXP, REALSXP, REALSXP,
                                REALSXP, REALSXP, REALSXP, REALSXP,
                                STRSXP,  REALSXP, INTSXP,  INTSXP};

  for (i = 0; i < NARGS; i++)
    PROT2(*args[i] = coerceVector(*args[i], type[i]), *numprot);

  S_len_m1 = length(*S) - 1;

  if ( S_len_m1 <= 0 )
    error ("need at least 2 price observations");

  if (!IS_MATRIX(*K) && IS_NUMERIC(*K)) {
    
    numStrikes = length(*K);

    tmp = PROT2(duplicate(*K), *numprot);
    *K  = PROT2(allocMatrix(REALSXP, S_len_m1, numStrikes), *numprot);
    
    for (i = 0; i < S_len_m1; i++)
      for (j = 0; j < numStrikes; j++)
        REAL(*K)[i + S_len_m1 * j] = REAL(tmp)[j];

    warning("Replacing K by matrix(K, length(S)-1, length(K), byrow = T)");

  }


  if (!IS_MATRIX(*K))
    error ("K is not a matrix");

  numStrikes = ncols(*K);

  if (numStrikes <= 0)
    error ("need positive number of strikes");

  if (1 == length(*r))
    PROT2(*r = numrep(*r, S_len_m1), *numprot);

  if (1 == length(*vol))
    PROT2(*vol = numrep(*vol, S_len_m1), *numprot);
  
  if (1 == length(*q))
    PROT2(*q = numrep(*q, S_len_m1), *numprot);

  if (1 == length(*transCosts))
    PROT2(*transCosts = numrep(*transCosts, S_len_m1 + 1), *numprot);


  if (numStrikes != length(*optionTypes))
    error ("numStrikes != length(optionTypes)");

  if (S_len_m1 != nrows(*K) && S_len_m1 + 1 != nrows(*K))
    error ("S_len_m1 != nrows(K) && S_len_m1 + 1 != nrows(K)");
  if (numStrikes != length(*posn))
    error ("numStrikes != length(posn)");
  if (S_len_m1 != length(*tau) || S_len_m1 != length(*r))
    error ("S_len_m1 != length(tau) || S_len_m1 != length(r)");
  if (S_len_m1 != length(*vol) && S_len_m1 + 1 != length(*vol))
    error ("S_len_m1 != length(vol) && S_len_m1 + 1 != length(vol)");
  if (S_len_m1 != length(*q) || S_len_m1 + 1 != length(*transCosts))
    error("S_len_m1!=length(q)||S_len_m1 + 1!=length(transCosts)");

  if (numStrikes > MAX_NUM_OPTNS)
    error ("numStrikes > MAX_NUM_OPTNS");


}
#undef NARGS



#define ANS_LEN 4
SEXP bsDeltaHedge2(SEXP tau, SEXP S, SEXP K, 
                   SEXP vol, SEXP posn, SEXP r, 
                   SEXP q, SEXP hedgePeriod, 
                   SEXP optionTypes, SEXP transCosts,
                   SEXP reltc, SEXP returnDaily)

  /*******************************************************************
   *
   *  Description: Returns SEXP 'ans' of type list including
   *    elements pl and tc.
   *    In R ans$pl[i] and ans$tc[i] are the the P&L and transaction costs
   *    respectively from the basic BS delting
   *    hedging strategy performed from REAL(tau)[i - 1] days from
   *    expiration to expiration.  Rebalancing occurs
   *    every asInteger(hedgePeriod) days.
   *    Here parameters K, posn and optionTypes are matrices each of whose
   *    rows taken together correpond
   *    to an option combination on the corresponding time to expiration. 
   *    
   *  Status: Tested & OK.
   *
   *******************************************************************/

{

  long S_len_m1, hPeriod, numStrikes, relTc, retDly;

  char *names[ANS_LEN]  = {"PL", "totalTransCosts", "dailyPL", "deltas"};

  int numprot = 0;

  int m = 0, i = 0, j = 0;
  double *tmpDeltas = NULL;
  double *tmpDlyPL  = NULL;
  double *tmpK      = NULL;
  char *opt[MAX_NUM_OPTNS];
 
  SEXP PL = R_NilValue, deltas = R_NilValue, ans = R_NilValue, 
    dly_PL = R_NilValue, totalTransCosts = R_NilValue;


  bsDeltaHedge2_chk(&tau, &S, &K, &vol, 
                    &posn, &r, &q, &hedgePeriod, 
                    &optionTypes, &transCosts, &reltc, &returnDaily, 
                    &numprot);



  S_len_m1   = length(S) - 1;
  hPeriod    = asInteger(hedgePeriod);
  numStrikes = ncols(K);
  relTc      = asInteger(reltc);
  retDly     = asInteger(returnDaily);


  PROT2(PL = NEW_NUMERIC(S_len_m1), numprot);
  PROT2(totalTransCosts = NEW_NUMERIC(S_len_m1), numprot);
  PROT2(ans = NEW_LIST(ANS_LEN), numprot);

  SET_NAMES(PL, GET_NAMES(tau));
  SET_NAMES(totalTransCosts, GET_NAMES(tau));


  if (retDly) {

    PROT2(deltas = NEW_LIST(S_len_m1), numprot);
    PROT2(dly_PL = NEW_LIST(S_len_m1), numprot);

    SET_NAMES(deltas, GET_NAMES(SET_LENGTH(tau, S_len_m1)));
    SET_NAMES(dly_PL, GET_NAMES(SET_LENGTH(tau, S_len_m1)));

  }


  dsetzero(REAL(PL), length(PL));
  dsetzero(REAL(totalTransCosts), length(totalTransCosts));

  tmpDeltas = (double *) R_alloc(S_len_m1 + 1, sizeof(double));
  tmpDlyPL  = (double *) R_alloc(S_len_m1 + 1, sizeof(double));
  tmpK      = (double *) R_alloc(numStrikes,   sizeof(double));


  if ( numStrikes > MAX_NUM_OPTNS )
    error ("too many strikes!");


  for (j = 0; j < numStrikes; j++) {
     opt[j] = (char *) R_alloc (1, sizeof(char));
    *opt[j] = *CHAR(GET_ELT(optionTypes, j));
  }



  for (i = 0; i < S_len_m1; i++) {
    
    m = S_len_m1 - i;

    dsetzero(tmpDeltas, S_len_m1 + 1); 
    dsetzero(tmpDlyPL,  S_len_m1 + 1);

    if (retDly) {
      SET_ELT(deltas, i, NEW_NUMERIC(m + 1));
      SET_ELT(dly_PL, i, NEW_NUMERIC(m + 1));
    }

    dmatrow1(K, i, tmpK);

    bsdh1(m, numStrikes, REAL(tau) + i, 
          REAL(S) + i, tmpK, REAL(vol) + i, 
          REAL(posn), REAL(r)[i], REAL(q)[i], 
          hPeriod, opt, REAL(transCosts) + i, 
          relTc, tmpDeltas, REAL(PL) + i, 
          tmpDlyPL, REAL(totalTransCosts) + i);


    if (retDly) {
      memcpy(REAL(GET_ELT(deltas, i)), tmpDeltas, (m + 1) * sizeof(double));
      memcpy(REAL(GET_ELT(dly_PL, i)), tmpDlyPL,  (m + 1) * sizeof(double));
    }

  }


  SET_ELT(ans, 0, PL);
  SET_ELT(ans, 1, totalTransCosts);

  if (retDly) {
    SET_ELT(ans, 2, dly_PL);
    SET_ELT(ans, 3, deltas);
  }

  set_names(ans, names);

  UNPROT2;

  return ans;

}
#undef ANS_LEN






#define NARGS 12
static void bsDeltaHedge3_chk(SEXP *tau, SEXP *S, SEXP *K, 
                              SEXP *vol, SEXP *posn, SEXP *r,
                              SEXP *q, SEXP *hedgePeriod, SEXP *optionTypes,
                              SEXP *transCosts, SEXP *reltc, 
                              SEXP *returnDaily,
                              int *numprot)

{
  int i, listLen;

  SEXP *args[NARGS] = {tau, S, K, 
                       vol, posn, r, 
                       q, hedgePeriod, optionTypes, 
                       transCosts, reltc, returnDaily};

  const char *names[NARGS] = {"tau", "S", "K", 
                              "vol", "posn", "r",
                              "q", "hedgePeriod", "optionTypes", 
                              "transCosts", "reltc", "returnDaily"};

  for (i = 0; i < NARGS; i++)
    if (1 == length(*args[i]) && !isNewList(*args[i])) {

      PROT2(*args[i] = tolist(*args[i]), *numprot);
      warning ("Replacing argument %s with list(%s)", names[i], names[i]);

    }

  for (i = 0; i < NARGS; i++)
    PROT2(*args[i] = AS_LIST(*args[i]), *numprot);

  listLen = length(*tau);

  if (1 != listLen)
    for (i = 1; i < NARGS; i++)
      if (1 == length(*args[i])) {

        PROT2(*args[i] = listrep(*args[i], listLen), *numprot);

        warning ("Replacing %s with rep(%s, %d)", names[i], names[i], listLen);

      }


  for (i = 1; i < NARGS; i++)
    if (listLen != length(*args[i]))
      error ("Argument %s not of length %d", names[i], listLen);

}
#undef NARGS







SEXP bsDeltaHedge3(SEXP tau, SEXP S, SEXP K, 
                   SEXP vol, SEXP posn, SEXP r,
                   SEXP q, SEXP hedgePeriod, SEXP optionTypes,
                   SEXP transCosts, SEXP reltc, 
                   SEXP returnDaily)

  /*******************************************************************
   *
   *  Description: Function for applying bsDeltaHedge2 
   *    when tau, S, K, ...
   *    are vectors or lists each of whose elements are the appropriate
   *    type for bsDeltaHedge2.  Returns SEXP of type list
   *    each of whose objects are the corresponding P&L and 
   *    transaction costs.
   *
   *  Status: ?
   *
   *******************************************************************/

{
  int numprot = 0;

  long listLen;
  int i;

  SEXP ans = R_NilValue, tmp = R_NilValue;

  bsDeltaHedge3_chk(&tau,  &S, &K, &vol,
                    &posn, &r, &q, &hedgePeriod, 
                    &optionTypes, &transCosts, &reltc, 
                    &returnDaily, &numprot);

  listLen = length(tau);

  PROT2(ans = NEW_LIST(listLen), numprot);

  for (i = 0; i < listLen; i++) {

    tmp = bsDeltaHedge2(GET_ELT(tau, i), GET_ELT(S, i),
                        GET_ELT(K, i), GET_ELT(vol, i),
                        GET_ELT(posn, i), GET_ELT(r, i),
                        GET_ELT(q, i), GET_ELT(hedgePeriod, i),
                        GET_ELT(optionTypes, i),
                        GET_ELT(transCosts, i),
                        GET_ELT(reltc, i),
                        GET_ELT(returnDaily, i));

    SET_ELT(ans, i, duplicate(tmp));

  }

  SET_NAMES(ans, GET_NAMES(S));

  UNPROT2;

  return ans;

}
#undef MAX_NUM_OPTNS
#undef THREE


#ifdef NOT_YET


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Functions below are for combination delta & vega-hedging.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/





static void dvHedge_chk(SEXP *posn, SEXP *posnStrikes, SEXP *optionTypes,
                        SEXP *T, SEXP *S, SEXP *vols, 
                        SEXP *volStrikes, SEXP *r, SEXP *q,
                        SEXP *hedgePeriod, SEXP *undTransCosts,
                        SEXP *optTransCosts, SEXP *isRelUndCosts,
                        SEXP *isRelOptCosts, int *numprot)
{
  int ok, i;
  int posnLen;
  long T_len, S_len;

  SEXP *args[] = {posn, posnStrikes, optionTypes,
                  T, S, vols,
                  volStrikes, r, q,
                  hedgePeriod, undTransCosts,
                  optTransCosts, isRelUndCosts,
                  isRelOptCosts};

  const SEXPTYPE type[] = {REALSXP, REALSXP, STRSXP,
                           REALSXP, REALSXP, REALSXP,
                           REALSXP, REALSXP, REALSXP,
                           REALSXP, INTSXP,
                           INTSXP};

  for (i = 0; i < 12; i++)
    PROT2(*args[i] = coerceVector(*args[i], type[i]), *numprot);

  posnLen = length(*posn);

  if (length(*posnStrikes) != posnLen || length(*optionTypes) != posnLen)
    error ("posn and optionType must equal posnStrikes in length");
  
  T_len = length(*T);
  S_len = length(*S);

  if (S_len != T_len + 1 )
    error ("T must be 1 less than S in length");

  if (S_len < 2)
    error ("Price series must have length >= 2");

  if (!isMatrix(*vols) || !isMatrix(*volStrikes))
    error ("vols and volStrikes must be matrices");

  ok = (nrows(*vols) == nrows(*volStrikes) && ncols(*vols) == ncols(*volStrikes));

  if (!ok)
    error ("Non-conforming vol and volStrike matrices");

  ok = (nrows(*vols) == T_len || nrows(*vols) == S_len);

  if (!ok)
    error ("Time-series matrices vols and volStrikes have wrong # of rows");

  if (length(*r) > 1 || length(*q) > 1)
    warning ("Only using 1st values of r and q");

  for (i = 0; i < T_len - 1; i++)
    if ( REAL(*T)[i] <= REAL(*T)[i+1] )
      error ("times to expiration must be strictly decreasing");

}













static void getPosnDeltasVegas(SEXP posn, SEXP posnStrikes, SEXP optionTypes,
                               SEXP T, SEXP S, SEXP vols, 
                               SEXP volStrikes, SEXP r, SEXP q,
                               SEXP atmVol, SEXP atmStrike, 
                               SEXP deltas, SEXP vegas, 
                               double *sigma, double *K)

  /*******************************************************************
   *
   *  Description: Utility function for dvHedge.  Computes deltas
   *    and vegas for given option position and also gets
   *    atm vols & strikes.
   *    
   *  Warning: Function does no argument checking as it assumes
   *    that this was done from dvHedge.
   *
   *******************************************************************/

{
  long i, j, ix, S_len = length(S), T_len = S_len - 1;
  int numVolCols = ncols(vols), numOptions = length(posn);
  double tau = 0, vol = 0, bsOut[3] = {0.0};

  /* Get deltas & vegas of base position */
  for (i = 0; i < S_len; i++) {

    if (i < T_len) {

      dmatrow1(volStrikes, i, K);
      dmatrow1(vols, i, sigma);

      ix = nearestIndex(REAL(S)[i], 0, K, numVolCols);

      REAL(atmVol)[i] = sigma[ix];
      REAL(atmStrike)[i] = K[ix];

    }

    tau = (i < T_len ? REAL(T)[i] / 365.0 : 0);
    vol = (i < T_len ? REAL(atmVol)[i] : 1e-8);

    REAL(deltas)[i] = REAL(vegas)[i] = 0;

    for (j = 0; j < numOptions; j++) { 

      bs0(tau, REAL(S)[i], REAL(posnStrikes)[j], 
          vol, *REAL(r), *REAL(q), 
          *CHAR(STRING_ELT(optionTypes, j)), bsOut);

      REAL(deltas)[i] += REAL(posn)[j] * bsOut[1];
      REAL(vegas)[i]  += REAL(posn)[j] * bsOut[2];
      
    }
  }

}





static void getHedges(SEXP T, SEXP S, SEXP vols, SEXP volStrikes,
                      SEXP r, SEXP q, SEXP atmVol, SEXP atmStrike, 
                      SEXP deltas, SEXP vegas, SEXP strdVol,
                      SEXP prevStrdVol, SEXP strdPrice, SEXP prevStrdPrice,
                      SEXP strdPosn, SEXP strdStrike, SEXP undPosn,
                      double *sig, double *K, int period)

  /*******************************************************************
   *
   *  Description: Utility function for dvHedge.  Computes
   *    underlying and straddle positions for delta & vega
   *    hedging strategy.
   *    
   *  Warning: Function does no argument checking as it assumes
   *    that this was done from dvHedge.
   *
   *******************************************************************/

{
  long i, ix, dayCtr = 0, S_len = length(S), T_len = S_len - 1;
  const int numVolCols = ncols(vols);
  double tau = 0, bsOut[3] = {0};

  REAL(prevStrdVol)[0] = REAL(prevStrdPrice)[0] = NA_REAL;

  for (i = 0; i < T_len; i++) {

    dmatrow1(volStrikes, i, K);
    dmatrow1(vols, i, sig);

    if (dayCtr >= period || !i) { 

      dayCtr = 0;

      REAL(strdStrike)[i] = REAL(atmStrike)[i];
      REAL(strdVol)[i] = REAL(atmVol)[i];
      
      bs0(REAL(T)[i] / 365.0, REAL(S)[i], REAL(strdStrike)[i], 
          REAL(strdVol)[i], *REAL(r), *REAL(q), 's', bsOut);
      
      REAL(strdPosn)[i] = -REAL(vegas)[i] / bsOut[2];
      REAL(undPosn)[i] = -(REAL(deltas)[i] + REAL(strdPosn)[i] * bsOut[1]);
      REAL(strdPrice)[i] = bsOut[0];

    }
    else {

      REAL(strdStrike)[i] = REAL(strdStrike)[i-1];
      REAL(strdPosn)[i] = REAL(strdPosn)[i-1];
      REAL(undPosn)[i] = REAL(undPosn)[i-1];

      ix = nearestIndex(REAL(strdStrike)[i], 0, K, numVolCols);
      
      REAL(strdVol)[i] = sig[ix];

      bs0(REAL(T)[i] / 365.0, REAL(S)[i], REAL(strdStrike)[i], 
          REAL(strdVol)[i], *REAL(r), *REAL(q), 's', bsOut);

      REAL(strdPrice)[i] = bsOut[0];

    }

    if (i) {
      ix = nearestIndex(REAL(strdStrike)[i-1], 0, K, numVolCols);
      REAL(prevStrdVol)[i] = sig[ix];
    
      tau = REAL(T)[i] / 365.0;

      bs0(tau, REAL(S)[i], REAL(strdStrike)[i-1],
          REAL(prevStrdVol)[i], *REAL(r), *REAL(q), 's', bsOut);
    
      REAL(prevStrdPrice)[i] = bsOut[0];
    }

    dayCtr += i < T_len - 2 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len - 1];

  }


  bs0(0, REAL(S)[S_len-1], REAL(strdStrike)[S_len-2], 0.01, 0, 0, 's', bsOut);

  REAL(prevStrdPrice)[S_len-1] = bsOut[0];
  REAL(strdVol)[S_len-1] = REAL(strdPrice)[S_len-1] = NA_REAL;
  REAL(strdStrike)[S_len-1] = NA_REAL;
  REAL(undPosn)[S_len-1] = 
    -REAL(deltas)[S_len-1] - bsOut[1] * REAL(strdPosn)[S_len-2];
  REAL(strdPosn)[S_len-1] = 0;

}





static void getHedgePLandTC(SEXP T, SEXP S, SEXP r, SEXP q, 
                            SEXP strdPrice, SEXP prevStrdPrice, SEXP strdPosn, 
                            SEXP strdStrike, SEXP undPosn, 
                            double tcUnd, double tcOpt,
                            int isRelUndTc, int isRelOptTc,
                            double *deltaHedgePL, double *vegaHedgePL,
                            double *totalUndTC, double *totalOptTC)

  /*******************************************************************
   *
   *  Description: Utility function for dvHedge.
   *    Computes hedging P&L and transaction costs.
   *
   *  Warning: Function does no argument checking as it assumes
   *    that this was done from dvHedge.
   *
   *******************************************************************/

{
  long i, S_len = length(S), T_len = S_len - 1;
  double undTCtoday, optTCtoday, strdPrice1, strdPrice2, dt;
  double ir = asReal(r), div = asReal(q), PLtoday = 0, discFac = 1.0;
  double tc1, tc2;

  *deltaHedgePL = *vegaHedgePL = 0;

  /* Get hedging P&L */
  for (i = 0; i < T_len; i++) {
    
    dt  = i < T_len - 1 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len-1];
    dt /= 365.0;

    discFac *= exp(-ir * dt);

    PLtoday  = REAL(undPosn)[i] * (REAL(S)[i+1] - REAL(S)[i]);
    PLtoday += REAL(undPosn)[i] * REAL(S)[i] * expm1((div-ir)*dt);

    *deltaHedgePL += discFac * PLtoday;

    strdPrice1 = exp(ir * dt) * REAL(strdPrice)[i];
    strdPrice2 = REAL(prevStrdPrice)[i+1];

    PLtoday  = REAL(strdPosn)[i] * (strdPrice2 - strdPrice1);

    *vegaHedgePL += discFac * PLtoday;

  }


  /* Get transaction costs */
  undTCtoday  = ABS(REAL(undPosn)[0]) * tcUnd;
  undTCtoday *= isRelUndTc ? REAL(S)[0] : 1;

  optTCtoday  = ABS(REAL(strdPosn)[0]) * tcOpt;
  optTCtoday *= isRelOptTc ? REAL(strdPrice)[0] : 1;
    
  *totalUndTC = undTCtoday;
  *totalOptTC = optTCtoday;

  for (i = 1, discFac = 1.0; i < T_len; i++) {

    dt  = (REAL(T)[i-1] - REAL(T)[i]) / 365.0;
    discFac *= exp(-ir * dt);

    undTCtoday  = ABS(REAL(undPosn)[i] - REAL(undPosn)[i-1]) * tcUnd;
    undTCtoday *= isRelUndTc ? REAL(S)[i] : 1;

    if ( ABS(REAL(strdStrike)[i-1] - REAL(strdStrike)[i]) < DOUBLE_EPS ) {

      optTCtoday = ABS(REAL(strdPosn)[i] - REAL(strdPosn)[i-1]) * tcOpt;
      optTCtoday *= isRelOptTc ? REAL(strdPrice)[i] : 1.0;

    }
    else {

      tc1 = ABS(REAL(strdPosn)[i]) * tcOpt;
      tc2 = ABS(REAL(strdPosn)[i-1]) * tcOpt;
      tc1 *= isRelOptTc ? REAL(strdPrice)[i] : 1;
      tc2 *= isRelOptTc ? REAL(prevStrdPrice)[i] : 1;

      optTCtoday = tc1 + tc2;

    }

    *totalUndTC += discFac * undTCtoday;
    *totalOptTC += discFac * optTCtoday;

  }

  discFac *= exp(-ir * REAL(T)[T_len-1] / 365);

  undTCtoday  = ABS(REAL(undPosn)[T_len] - REAL(undPosn)[T_len-1]) * tcUnd;
  undTCtoday *= isRelUndTc ? REAL(S)[T_len] : 1;

  *totalUndTC += discFac * undTCtoday;

}




#define ANS_LEN 13
SEXP dvHedge(SEXP posn, SEXP posnStrikes, SEXP optionTypes,
             SEXP hedgePeriod, SEXP T, SEXP S, SEXP vols, 
             SEXP volStrikes, SEXP r, SEXP q, 
             SEXP undTransCosts, 
             SEXP optTransCosts, SEXP isRelUndCosts, 
             SEXP isRelOptCosts)

  /*******************************************************************
   *
   *  Description: Computes the P&L and transaction costs for vega
   *    & delta-hedging an option position.  The vega hedging is
   *    done using ATM straddles and rebalanced using vertical 
   *    straddle swaps.
   *    
   *  Parameters:
   *
   *    posn - Vector giving positions in various options.
   *    posnStrikes - Vector of position strikes.
   *    optionTypes - Vector of type "character" giving the options
   *        in which the position is held.
   *    T - Time series of times-to-expiration at which the observations
   *        of underlying price and volatilities et. al. take place.
   *        Assumes each T[i] is in days.
   *    S - Time series of underlying prices.
   *    vols - Time series matrix of volatilities.
   *    volStrikes - Time series matrix of strikes corresponding to
   *        vols.
   *    r - Per annum risk-free rate.
   *    q - Per annum dividend yield.
   *    hedgePeriod - Number of days between
   *        hedge rebalances.
   *    undTransCosts - transaction costs per contract traded in
   *        the underlying, either relative or absolute.
   *    optTransCosts - transaction costs per contract for options.
   *    isRelUndCosts - Nonzero iff undTransCosts are relative.
   *    isRelOptCosts - Nonzero iff optTransCosts are relative.
   *
   *  Status: Not finished.
   *
   *******************************************************************/

{

  int numprot = 0;
  long numOptions, T_len, S_len, numVolCols;

  const int period = asInteger(hedgePeriod);
  const int isRelUndTc = asInteger(isRelUndCosts);
  const int isRelOptTc = asInteger(isRelOptCosts);
  const double tcUnd = asReal(undTransCosts);
  const double tcOpt = asReal(optTransCosts);

  double *K, *sig;
  double deltaHedgePL = 0, vegaHedgePL = 0;
  double totalUndTC = 0, totalOptTC = 0;

  char *names[ANS_LEN] = {"strdPosn", "undPosn", "strdVol", 
                          "prevStrdVol", "strdStrike", "strdPrice",
                          "prevStrdPrice", "deltas", "vegas",
                          "undPL", "strdPL", "undTC", "strdTC"};

  SEXP strdVol, prevStrdVol, strdStrike,
    strdPrice, prevStrdPrice, atmVol, atmStrike,
    deltas, vegas, strdPosn, undPosn, ans;


  dvHedge_chk(&posn, &posnStrikes, &optionTypes,
                &T, &S, &vols,
              &volStrikes, &r, &q,
              &hedgePeriod, &undTransCosts,
              &optTransCosts, &isRelUndCosts,
              &isRelOptCosts, &numprot);


  S_len = length(S);
  T_len = S_len - 1;
  numOptions = length(posn);
  numVolCols = ncols(vols);


  PROT2(strdVol = NEW_NUMERIC(S_len), numprot);
  PROT2(prevStrdVol = NEW_NUMERIC(S_len), numprot);
  PROT2(strdStrike = NEW_NUMERIC(S_len), numprot);
  PROT2(deltas = NEW_NUMERIC(S_len), numprot);
  PROT2(vegas  = NEW_NUMERIC(S_len), numprot);
  PROT2(undPosn = NEW_NUMERIC(S_len), numprot);
  PROT2(strdPosn = NEW_NUMERIC(S_len), numprot);
  PROT2(strdPrice = NEW_NUMERIC(S_len), numprot);
  PROT2(prevStrdPrice = NEW_NUMERIC(S_len), numprot);
  PROT2(atmVol = NEW_NUMERIC(S_len), numprot);
  PROT2(atmStrike = NEW_NUMERIC(S_len), numprot);
  PROT2(ans    = NEW_LIST(ANS_LEN), numprot);


  K = (double *) R_alloc(numVolCols, sizeof(double));
  sig = (double *) R_alloc(numVolCols, sizeof(double));

  getPosnDeltasVegas(posn, posnStrikes, optionTypes,
                     T, S, vols, 
                     volStrikes, r, q,
                     atmVol, atmStrike, deltas, vegas, sig, K);


  getHedges(T, S, vols, volStrikes,
            r, q, atmVol, atmStrike, 
            deltas, vegas, strdVol,
            prevStrdVol, strdPrice, prevStrdPrice,
            strdPosn, strdStrike, undPosn, sig, 
            K, period);


  getHedgePLandTC(T, S, r, q, 
                  strdPrice, prevStrdPrice, strdPosn, 
                  strdStrike, undPosn, tcUnd, tcOpt,
                  isRelUndTc, isRelOptTc,
                  &deltaHedgePL, &vegaHedgePL,
                  &totalUndTC, &totalOptTC);

  SET_ELT(ans, 0, strdPosn);
  SET_ELT(ans, 1, undPosn);
  SET_ELT(ans, 2, strdVol);
  SET_ELT(ans, 3, prevStrdVol);
  SET_ELT(ans, 4, strdStrike);
  SET_ELT(ans, 5, strdPrice);
  SET_ELT(ans, 6, prevStrdPrice);
  SET_ELT(ans, 7, deltas);
  SET_ELT(ans, 8, vegas);
  SET_ELT(ans, 9, ScalarReal(deltaHedgePL));
  SET_ELT(ans, 10, ScalarReal(vegaHedgePL));
  SET_ELT(ans, 11, ScalarReal(totalUndTC));
  SET_ELT(ans, 12, ScalarReal(totalOptTC));

  set_names(ans, names);

  UNPROT2;

  return ans;

}
#undef ANS_LEN








static void dvHedge2_chk(SEXP *posn, SEXP *posnStrikes, SEXP *optionTypes,
                         SEXP *hedgePeriod, SEXP *T, SEXP *S, SEXP *vols, 
                         SEXP *volStrikes, SEXP *r, SEXP *q, 
                         SEXP *undTransCosts, SEXP *optTransCosts, 
                         SEXP *isRelUndCosts, SEXP *isRelOptCosts,
                         int *numprot)
{
  long len, T_len, S_len;
  int ok;

  PROT2(*T = AS_NUMERIC(*T), *numprot);

  len = length(*posnStrikes);

  S_len = length(*S);
  T_len = length(*T);

  ok = (T_len == S_len - 1);

  if (!ok)
    error ("bleh");

  ok = isMatrix(*posnStrikes) && isMatrix(*vols) && 
    isMatrix(*volStrikes) && 
    nrows(*posnStrikes) == T_len &&
    nrows(*volStrikes) == T_len &&
    nrows(*vols) == T_len && 
    ncols(*posnStrikes) == length(*posn) &&
    length(*posn) == length(*optionTypes);

  if (!ok)
    error ("bleh bleh");

  if (1 == length(*r))
    PROT2(*r = numrep(*r, T_len), *numprot);
  if (1 == length(*q))
    PROT2(*q = numrep(*q, T_len), *numprot);

  ok = length(*r) == T_len && length(*q) == T_len;

  if (!ok)
    error ("bleh bleh bleh");

  if (S_len < 2)
    error ("Price series must have length >= 2");

}







#define ANS_LEN 4
SEXP dvHedge2(SEXP posn, SEXP posnStrikes, SEXP optionTypes,
              SEXP hedgePeriod, SEXP T, SEXP S, SEXP vols, 
              SEXP volStrikes, SEXP r, SEXP q, 
              SEXP undTransCosts, SEXP optTransCosts, 
              SEXP isRelUndCosts, SEXP isRelOptCosts)
{

  int numprot = 0;

  int numVolCols, posnLen;
  long T_len, S_len, i;
  SEXP strikesTmp, ans, tmp, undPL, strdPL, undTC, strdTC, r0, q0;

  char *names[ANS_LEN] = {"undPL", "strdPL", "undTC", "strdTC"};

  dvHedge2_chk(&posn, &posnStrikes, &optionTypes,
               &hedgePeriod, &T, &S, &vols,
               &volStrikes, &r, &r,
               &undTransCosts, &optTransCosts,
               &isRelUndCosts, &isRelOptCosts,
               &numprot);

  S_len = length(S);
  T_len = S_len - 1;
  posnLen = length(posn);
  numVolCols = ncols(vols);

  PROT2(ans = NEW_LIST(ANS_LEN), numprot);
  PROT2(strikesTmp = NEW_NUMERIC(posnLen), numprot);
  PROT2(undPL  = NEW_NUMERIC(T_len), numprot);
  PROT2(strdPL = NEW_NUMERIC(T_len), numprot);
  PROT2(undTC  = NEW_NUMERIC(T_len), numprot);
  PROT2(strdTC = NEW_NUMERIC(T_len), numprot);
  PROT2(r0 = NEW_NUMERIC(1), numprot);
  PROT2(q0 = NEW_NUMERIC(1), numprot);


  set_names(ans, names);

  dsetna(REAL(undPL), T_len);
  dsetna(REAL(strdPL), T_len);
  dsetna(REAL(undTC), T_len);
  dsetna(REAL(strdTC), T_len);

  for (i = 0; i < T_len; i++) {
    
    dmatrow1(posnStrikes, i, REAL(strikesTmp));

    REAL(r0)[0] = REAL(r)[i];
    REAL(q0)[0] = REAL(q)[i];
    
    tmp = dvHedge(posn, strikesTmp, optionTypes, hedgePeriod,
                  PROTECT(dSubVector(T, i, T_len - 1)),
                  PROTECT(dSubVector(S, i, S_len - 1)),
                  PROTECT(dSubMatrix(vols, i, T_len - 1, 0, numVolCols - 1)),
                  PROTECT(dSubMatrix(volStrikes, i, T_len - 1, 0, numVolCols - 1)),
                  r0, q0,
                  undTransCosts, optTransCosts,
                  isRelUndCosts, isRelOptCosts);

    PROTECT(tmp);

    REAL(undPL)[i]  = *REAL(getListElt(tmp, "undPL"));
    REAL(strdPL)[i] = *REAL(getListElt(tmp, "strdPL"));
    REAL(undTC)[i]  = *REAL(getListElt(tmp, "undTC"));
    REAL(strdTC)[i] = *REAL(getListElt(tmp, "strdTC"));


    UNPROTECT(5);

  }

  SET_ELT(ans, 0, undPL);
  SET_ELT(ans, 1, strdPL);
  SET_ELT(ans, 2, undTC);
  SET_ELT(ans, 3, strdTC);



  UNPROT2;

  return ans;

}
#undef ANS_LEN


#endif
