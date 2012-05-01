#include "uri.h"
#include <stdlib.h>

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following section contains code for delta-hedging of spread
**  options.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

static int sprdOptDh1_chk(SEXP *T, SEXP *S1, SEXP *S2,
                          SEXP *K, SEXP *vol1, SEXP *vol2,
                          SEXP *rho, SEXP *r, SEXP *q1,
                          SEXP *q2, SEXP *opt, SEXP *nGrdPts,
                          int *numprot)
{
#define NNARGS 10
    int i, ok, allFinite;
    long t, T_len, S1_len, S2_len;
    SEXP *args[NNARGS+1] = {T, S1, S2, K, vol1, vol2,
                            rho, r, q1, q2, opt};

    for (i = 0; i < NNARGS; i++)
        if (0 == args[i])
            error("Null pointer.");

    for (i = 0; i < NNARGS; i++)
        ENSURE_NUMERIC(*args[i], *numprot);

    ENSURE_CHAR(*opt, *numprot);


    T_len = length(*T);
    S1_len = length(*S1);
    S2_len = length(*S2);
    if (T_len != S1_len - 1 || S1_len != S2_len)
        error ("Length problems");

    for (i = 3; i < NNARGS; i++)
        if (1 == length(*args[i]))
            PROT2(*args[i] = numrep(*args[i], T_len), *numprot);

    if (1 == length(*opt))
        PROT2(*opt = charrep(*opt, T_len), *numprot);

    ok = T_len == length(*opt);
    for (i = 3; i < NNARGS; i++)
        ok *= T_len == length(*args[i]);

    if (!ok)
        error ("More length problems");

    allFinite = R_FINITE(asReal(*nGrdPts));

    for (i = 0; i < NNARGS && allFinite; i++)
        for (t = 0; t < length(*args[i]) && allFinite; t++)
            allFinite *= R_FINITE(REAL(*args[i])[t]);

    return allFinite;
#undef NNARGS
}





/*
**  Description: Computes the delta-hedging P&L for a single spread
**  option.
*/

SEXP sprdOptDh1(SEXP T, SEXP S1, SEXP S2, 
                SEXP K, SEXP vol1, SEXP vol2,
                SEXP rho, SEXP r, SEXP q1,
                SEXP q2, SEXP opt, SEXP nGrdPts)
{
    int numprot = 0;
    int allFinite;
    long i, T_len;
    double d1, d2, PL, dt, PLtoday, discFac, tmp[PEARSON0_PTR_LEN];

    allFinite = sprdOptDh1_chk(&T, &S1, &S2, &K, &vol1, &vol2,
                               &rho, &r, &q1, &q2, &opt, &nGrdPts,
                               &numprot);

    if (!allFinite) {
        UNPROT2; 
        return ScalarReal(NA_REAL);
    }

    T_len = length(T);
    PL = 0;
    discFac = 1;

    for (i = 0; i < T_len && R_FINITE(PL); i++) {

        pearson0(REAL(T)[i]/365, REAL(S1)[i], REAL(S2)[i],
                 REAL(K)[i], REAL(vol1)[i], REAL(vol2)[i], 
                 REAL(rho)[i], REAL(r)[i], REAL(q1)[i], 
                 REAL(q2)[i], *CHAR(GET_ELT(opt, i)),
                 asInteger(nGrdPts), /* calcDeltas = */ 1,
                 tmp);

        d1 = tmp[1];
        d2 = tmp[2];

        dt = i < T_len - 1 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len-1];
        dt /= 365;

        if (dt <= 0)
            error ("Can't have negative time increments");

        discFac *= exp(-REAL(r)[i]*dt);

        PLtoday  = -d1 * (REAL(S1)[i+1] - REAL(S1)[i]);
        PLtoday += -d1 * REAL(S1)[i] * expm1((REAL(q1)[i]-REAL(r)[i])*dt);
        PLtoday += -d2 * (REAL(S2)[i+1] - REAL(S2)[i]);
        PLtoday += -d2 * REAL(S2)[i] * expm1((REAL(q2)[i]-REAL(r)[i])*dt);

        PL += discFac * PLtoday;

    }
    UNPROT2;
    return ScalarReal(PL);
}




static void listSprdOptDh1_chk(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, 
                               SEXP *vol1, SEXP *vol2, SEXP *rho, SEXP *r, 
                               SEXP *q1, SEXP *q2, SEXP *opt, SEXP *nGrdPts, 
                               int *numprot)
{
#define NARGS 11
    long i, N = 0;
    SEXP *args[NARGS] = {tau, S1, S2, K, 
                         vol1, vol2, rho, r,
                         q1, q2, opt};


    for (i = 0; i < NARGS; i++) {
    
        if (!isNewList(*args[i]))
            PROT2(*args[i] = tolist(*args[i]), *numprot);

        N = MAX(N, length(*args[i]));
  
    }

    if (1 < N)
        for (i = 0; i < NARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = listrep(*args[i], N), *numprot);

  
#undef NARGS
}



SEXP listSprdOptDh1(SEXP tau, SEXP S1, SEXP S2, 
                    SEXP K, SEXP vol1, SEXP vol2,
                    SEXP rho, SEXP r, SEXP q1,
                    SEXP q2, SEXP opt, SEXP nGrdPts)
{
    int numprot = 0;
    long i, N;

    SEXP ans, tmp;

    listSprdOptDh1_chk(&tau, &S1, &S2, &K, &vol1, &vol2, &rho,
                       &r, &q1, &q2, &opt, &nGrdPts,
                       &numprot);

    N = length(tau);
    PROT2(ans = NEW_NUMERIC(N), numprot);
  
    for (i = 0; i < N; i++) {

        tmp = sprdOptDh1(GET_ELT(tau, i), GET_ELT(S1, i), GET_ELT(S2, i),
                         GET_ELT(K, i), GET_ELT(vol1, i), GET_ELT(vol2, i),
                         GET_ELT(rho, i), GET_ELT(r, i), GET_ELT(q1, i),
                         GET_ELT(q2, i), GET_ELT(opt, i), nGrdPts);

        REAL(ans)[i] = REAL(tmp)[0];

    }

    SET_NAMES(ans, GET_NAMES(tau));
    UNPROT2;
    return ans;
}







/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following section contains code for delta-hedging spread options
**  when there is a quoted implied volatility skew for each of the
**  legs.  The convention is that at time t the volatility used
**  for leg S2 corresponds to S1 + K and vice-versa.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


static void sprdOptDh2_chk(SEXP *T, SEXP *S1, SEXP *S2,
                           SEXP *K, SEXP *vols1, SEXP *vols2,
                           SEXP *volStrikes1, SEXP *volStrikes2,
                           SEXP *rho, SEXP *r, SEXP *q1,
                           SEXP *q2, SEXP *opt, SEXP *nGrdPts,
                           int *numprot)
{
#define N_VEC_ARGS 8
#define N_MAT_ARGS 4

    long t, T_len, S_len;
  
    int i, ok, allFinite;

    SEXP *vecArgs[N_VEC_ARGS] = {T, S1, S2, K, rho, r, q1, q2};
    SEXP *matArgs[N_MAT_ARGS] = {vols1, vols2, volStrikes1, volStrikes2};

    for (i = 0; i < N_VEC_ARGS; i++)
        ENSURE_NUMERIC(*vecArgs[i], *numprot);

    ENSURE_CHAR(*opt, *numprot);

    for (i = 0; i < N_MAT_ARGS; i++)
        PROT2(*matArgs[i] = asRealMatrix(*matArgs[i]), *numprot);

    if (length(*S1) != length(*S2))
        error ("Unequal price series lengths");

    T_len = length(*T);
    S_len = length(*S1);

    if (S_len != T_len + 1)
        error("price series need be 1 longer than time-to-exp. series");

    for (i = 4; i < N_VEC_ARGS; i++)
        if (1 == length(*vecArgs[i]))
            PROT2(*vecArgs[i] = numrep(*vecArgs[i], T_len), *numprot);

    for (i = 0; i < N_MAT_ARGS; i++)
        if (nrows(*matArgs[i]) != T_len)
            error ("not enough rows");

    if (ncols(*vols1) != ncols(*volStrikes1))
        error ("column mismatch");

    if (ncols(*vols2) != ncols(*volStrikes2))
        error ("column mismatch");

    if (1 == length(*opt))
        PROT2(*opt = charrep(*opt, T_len), *numprot);

    ok = (T_len == length(*opt));

    for (i = 4; i < N_VEC_ARGS && ok; i++)
        ok *= (T_len == length(*vecArgs[i]));

    for (i = 0; i < N_MAT_ARGS && ok; i++)
        ok *= (T_len == nrows(*matArgs[i]));

    if (!ok)
        error ("More length problems");

    allFinite = R_FINITE(asReal(*nGrdPts));

    for (i = 0; i < N_VEC_ARGS && allFinite; i++)
        for (t = 0; t < length(*vecArgs[i]) && allFinite; t++)
            allFinite *= R_FINITE(REAL(*vecArgs[i])[t]);

    for (i = 0; i < N_MAT_ARGS && allFinite; i++)
        for (t = 0; t < length(*matArgs[i]) && allFinite; t++)
            allFinite *= R_FINITE(REAL(*matArgs[i])[t]);

    if (!allFinite)
        error ("ooh");
#undef N_VEC_ARGS
#undef N_MAT_ARGS
  
}





static void getDeltas(SEXP T, SEXP S1, SEXP S2, double K,
                      SEXP vols1, SEXP vols2, SEXP rho,
                      SEXP volStrikes1, SEXP volStrikes2,
                      SEXP r, SEXP q1, SEXP q2,
                      char optionType, long nGrdPts,
                      SEXP d1, SEXP d2, 
                      double *sig1, double *sig2,
                      double *K1, double *K2)
{
    long i, T_len=length(T);
    int numVolCols1 = ncols(vols1), numVolCols2 = ncols(vols2);
    double vol1 = 0, vol2, X1, X2, pearsonOut[PEARSON0_PTR_LEN] = {0.0};

    /* Get deltas */
    for (i = 0; i < T_len; i++) {

        X1 = REAL(S1)[i] + K;
        X2 = DMAX(DOUBLE_EPS, REAL(S2)[i] - K);

        dmatrow1(volStrikes1, i, K1);
        dmatrow1(volStrikes2, i, K2);
        dmatrow1(vols1, i, sig1);
        dmatrow1(vols2, i, sig2);
      
        vol1 = near(X1, 0, K1, sig1, numVolCols1);
        vol2 = near(X2, 0, K2, sig2, numVolCols2);
    
        pearson0(REAL(T)[i]/365.0, REAL(S1)[i], REAL(S2)[i], K,
                 vol1, vol2, REAL(rho)[i],
                 REAL(r)[i], REAL(q1)[i], REAL(q2)[i],
                 optionType, nGrdPts,
                 /* calcDeltas = */ 1, pearsonOut);

        REAL(d1)[i] = pearsonOut[PEARSON0_D1_INDEX];
        REAL(d2)[i] = pearsonOut[PEARSON0_D2_INDEX];

    }

    REAL(d1)[T_len] = EXP_DELTA(REAL(S1)[T_len],REAL(S2)[T_len]-K,optionType);
    REAL(d2)[T_len] = EXP_DELTA(REAL(S2)[T_len],REAL(S1)[T_len]+K,optionType);
  
}





static void getDeltaHedges(SEXP T, SEXP d1, SEXP d2, SEXP undPosn1,
                           SEXP undPosn2, int period)
{
    long i=0, dayCtr = 0, T_len=length(T), S_len = T_len+1;

    for (i = 0; i < S_len; i++) {

        if (dayCtr >= period || !i) { 

            dayCtr = 0;

            REAL(undPosn1)[i] = -REAL(d1)[i];
            REAL(undPosn2)[i] = -REAL(d2)[i];
      
        }
        else {

            REAL(undPosn1)[i] = REAL(undPosn1)[i-1];
            REAL(undPosn2)[i] = REAL(undPosn2)[i-1];

        }

        dayCtr += i < T_len - 2 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len - 1];

    }

}



static void getHedgePLandTC(SEXP T, SEXP S1, SEXP S2, 
                            SEXP r, SEXP q1, SEXP q2,
                            SEXP undPosn1, SEXP undPosn2,
                            double tcUnd1, double tcUnd2,
                            int isRelUndTc1, int isRelUndTc2,
                            double *deltaHedgePL, double *totalTC)
{
    long i, T_len=length(T), S_len = T_len+1;
    double und1TCtoday, und2TCtoday, dt;
    double ir=asReal(r), div1=asReal(q1), div2=asReal(q2),PLtoday=0, discFac=1;


    if (!deltaHedgePL || !totalTC)
        error ("null pointers");

    *deltaHedgePL = *totalTC = 0;

    /* Get hedging P&L */
    for (i = 0; i < T_len; i++) {
    
        dt  = i < T_len - 1 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len-1];
        dt /= 365.0;

        discFac *= exp(-ir*dt);

        PLtoday  = REAL(undPosn1)[i] * (REAL(S1)[i+1] - REAL(S1)[i]);
        PLtoday += REAL(undPosn1)[i] * REAL(S1)[i] * expm1((div1-ir)*dt);
        PLtoday += REAL(undPosn2)[i] * (REAL(S2)[i+1] - REAL(S2)[i]);
        PLtoday += REAL(undPosn2)[i] * REAL(S2)[i] * expm1((div2-ir)*dt);    

        *deltaHedgePL += discFac * PLtoday;

    }

    /* Get transaction costs */
    und1TCtoday  = ABS(REAL(undPosn1)[0]) * tcUnd1;
    und1TCtoday *= isRelUndTc1 ? REAL(S1)[0] : 1;
    und2TCtoday  = ABS(REAL(undPosn2)[0]) * tcUnd2;
    und2TCtoday *= isRelUndTc2 ? REAL(S2)[0] : 1;
      
    *totalTC = und1TCtoday + und2TCtoday;

    for (i = 1, discFac = 1.0; i < S_len; i++) {

        dt  = i < T_len ? (REAL(T)[i-1] - REAL(T)[i])/365.0: REAL(T)[T_len-1]/365.0;
        discFac *= exp(-ir*dt);

        und1TCtoday  = ABS(REAL(undPosn1)[i] - REAL(undPosn1)[i-1]) * tcUnd1;
        und1TCtoday *= isRelUndTc1 ? REAL(S1)[i] : 1;
        und2TCtoday  = ABS(REAL(undPosn2)[i] - REAL(undPosn2)[i-1]) * tcUnd2;
        und2TCtoday *= isRelUndTc2 ? REAL(S2)[i] : 1;

        *totalTC += discFac * (und1TCtoday+und2TCtoday);

    }

}




SEXP sprdOptDh2(SEXP T, SEXP S1, SEXP S2,
                SEXP K, SEXP vols1, SEXP vols2,
                SEXP volStrikes1, SEXP volStrikes2,
                SEXP rho, SEXP r, SEXP q1,
                SEXP q2, SEXP opt, SEXP nGrdPts)
{
    int numprot = 0;
    long T_len, S_len, numVolCols1, numVolCols2;
    char optionType;
    double deltaHedgePL=0, totalTC=0;
    double *sig1, *sig2, *K1, *K2;
    SEXP d1, d2, undPosn1, undPosn2;

    sprdOptDh2_chk(&T, &S1, &S2, &K,
                   &vols1, &vols2, &volStrikes1, &volStrikes2,
                   &rho, &r, &q1, &q2,
                   &opt, &nGrdPts, &numprot);
  
    optionType = *CHAR(GET_ELT(opt,0));
    T_len = length(T);
    S_len = T_len+1;
    numVolCols1 = ncols(vols1);
    numVolCols2 = ncols(vols2);

    sig1 = (double *) R_alloc(numVolCols1, DSIZE);
    sig2 = (double *) R_alloc(numVolCols2, DSIZE);
    K1 = (double *) R_alloc(numVolCols1, DSIZE);
    K2 = (double *) R_alloc(numVolCols2, DSIZE);

    PROT2(d1 = NEW_NUMERIC(S_len),numprot);
    PROT2(d2 = NEW_NUMERIC(S_len),numprot);
    PROT2(undPosn1 = NEW_NUMERIC(S_len),numprot);
    PROT2(undPosn2 = NEW_NUMERIC(S_len),numprot);

    getDeltas(T, S1, S2, asReal(K),
              vols1, vols2, rho,
              volStrikes1, volStrikes2,
              r, q1, q2,
              optionType, asInteger(nGrdPts),
              d1, d2, 
              sig1, sig2, K1, K2);

    getDeltaHedges(T, d1, d2, undPosn1, undPosn2, /* period = */ 1);
  
    getHedgePLandTC(T, S1, S2, r, q1, q2,
                    undPosn1, undPosn2,
                    /* tcUnd1 */ 0, /* tcUnd2 */ 0,
                    /* isRelUndTc1 */ 0, /* isRelUndTc2 */ 0,
                    &deltaHedgePL, &totalTC);

    UNPROT2;
    return ScalarReal(deltaHedgePL);
}
                



static void listSprdOptDh2_chk(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, 
                               SEXP *vols1, SEXP *vols2, 
                               SEXP *volStrikes1, SEXP *volStrikes2,
                               SEXP *rho, SEXP *r, 
                               SEXP *q1, SEXP *q2, SEXP *opt, SEXP *nGrdPts, 
                               int *numprot)
{
#define NUM_ARGS 13
    long i, N = 0;
    SEXP *args[NUM_ARGS] = {tau, S1, S2, K, 
                            vols1, vols2, 
                            volStrikes1, volStrikes2,
                            rho, r,
                            q1, q2, opt};


    for (i = 0; i < NUM_ARGS; i++) {
    
        if (!isNewList(*args[i]))
            PROT2(*args[i] = tolist(*args[i]), *numprot);

        N = MAX(N, length(*args[i]));
  
    }

    if (1 < N)
        for (i = 0; i < NUM_ARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = listrep(*args[i], N), *numprot);

#undef NUM_ARGS
}





SEXP listSprdOptDh2(SEXP tau, SEXP S1, SEXP S2, 
                    SEXP K, SEXP vols1, SEXP vols2,
                    SEXP volStrikes1, SEXP volStrikes2,
                    SEXP rho, SEXP r, SEXP q1,
                    SEXP q2, SEXP opt, SEXP nGrdPts)
{

    int numprot = 0;
    long i, N;

    SEXP ans, tmp;

    listSprdOptDh2_chk(&tau, &S1, &S2, &K, 
                       &vols1, &vols2, 
                       &volStrikes1, &volStrikes2,
                       &rho,
                       &r, &q1, &q2, &opt, &nGrdPts,
                       &numprot);

    N = length(tau);
    PROT2(ans = NEW_NUMERIC(N), numprot);
    
    for (i = 0; i < N; i++) {

        tmp = sprdOptDh2(GET_ELT(tau, i), GET_ELT(S1, i), GET_ELT(S2, i),
                         GET_ELT(K, i), GET_ELT(vols1, i), GET_ELT(vols2, i),
                         GET_ELT(volStrikes1,i), GET_ELT(volStrikes2,i),
                         GET_ELT(rho, i), GET_ELT(r, i), GET_ELT(q1, i),
                         GET_ELT(q2, i), GET_ELT(opt, i), nGrdPts);

        REAL(ans)[i] = REAL(tmp)[0];

    }

    SET_NAMES(ans, GET_NAMES(tau));

    UNPROT2;

    return ans;

}




/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following section contains functions that compute the hedging
**  P&L as in the 1st function sprdOptDh1 however the daily P&L
**  is reported as are the daily option values, so that a daily net
**  P&L may be calculated.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


static int sprdOptDh3_chk(SEXP *T, SEXP *S1, SEXP *S2,
                          SEXP *K, SEXP *vol1, SEXP *vol2,
                          SEXP *rho, SEXP *r, SEXP *q1,
                          SEXP *q2, SEXP *opt, SEXP *nGrdPts,
                          int *numprot)
{
#define NNARGS 10
    int i, ok, allFinite;
    long t, T_len, S1_len, S2_len;
    SEXP *args[NNARGS+1] = {T, S1, S2, K, vol1, vol2,
                            rho, r, q1, q2, opt};

    for (i = 0; i < NNARGS; i++)
        ENSURE_NUMERIC(*args[i], *numprot);

    ENSURE_CHAR(*opt, *numprot);

    T_len = length(*T);
    S1_len = length(*S1);
    S2_len = length(*S2);
    if (T_len != S1_len - 1 || S1_len != S2_len)
        error ("Length problems");

    for (i = 4; i < NNARGS; i++)
        if (1 == length(*args[i]))
            PROT2(*args[i] = numrep(*args[i], T_len), *numprot);

    ok = 1;
    for (i = 4; i < NNARGS; i++)
        ok *= T_len == length(*args[i]);

    if (!ok)
        error ("More length problems");

    allFinite = R_FINITE(asReal(*nGrdPts));
    for (i = 0; i < NNARGS && allFinite; i++)
        for (t = 0; t < length(*args[i]) && allFinite; t++)
            allFinite *= R_FINITE(REAL(*args[i])[t]);
    return allFinite;
#undef NNARGS
}

/*
**  Description: Computes the daily delta-hedging P&L and option
**  values for a single spread option.
*/

SEXP sprdOptDh3(SEXP T, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho, 
                SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP nGrdPts, SEXP maxdt)
{
    enum {i_hedge, i_option, i_o1, i_o2, i_o3, i_h1, i_h2, i_h3, i_d1, i_d2, 
          i_tau, i_S1, i_S2, i_K, i_vol1, i_vol2, i_rho, ANS_LEN};

    const int max_dt = asInteger(maxdt);
    int numprot = 0;
    const long grdp = asInteger(nGrdPts);
    long i, T_len;
    double del1, del2, pltoday, discFac, den, *dt, tmp[PEARSON0_PTR_LEN];
    char option_type;
    char *rlist_names[ANS_LEN] = {"hedge", "option","o1","o2","o3","h1","h2","h3",
                                  "d1", "d2", "tau","S1","S2","K","vol1","vol2","rho"};
    SEXP dpl, dval, rlist, h1, o1, h2, o2, o3, h3, d1, d2, names;
    double X;

    PROT2(rlist = NEW_LIST(ANS_LEN), numprot);
    SET_ELT(rlist,i_rho,rho);

    sprdOptDh3_chk(&T, &S1, &S2, &K, &vol1, &vol2,
                   &rho, &r, &q1, &q2, &opt, &nGrdPts,
                   &numprot);

    T_len = length(T);
    discFac = 1;
    PROT2(dpl = NEW_NUMERIC(T_len+1), numprot);
    PROT2(dval = NEW_NUMERIC(T_len+1), numprot);
    PROT2(h1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(h2 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(h3 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(o1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(o2 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(o3 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(d1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(d2 = NEW_NUMERIC(T_len+1), numprot);

    dt = (double *) R_alloc(T_len, DSIZE);

    REAL(dpl)[0] = 0;
    REAL(h1)[0] = REAL(o1)[0] = REAL(h2)[0] = REAL(o2)[0] = NA_REAL;
    REAL(h3)[0] = REAL(o3)[0] = NA_REAL;

    option_type = *CHAR(GET_ELT(opt, 0));
    X = asReal(K);

    for (i = 0; i < T_len; i++) {

        pearson0(REAL(T)[i]/365, REAL(S1)[i], REAL(S2)[i],
                 X, REAL(vol1)[i], REAL(vol2)[i], 
                 REAL(rho)[i], REAL(r)[i], REAL(q1)[i], 
                 REAL(q2)[i], option_type,
                 grdp, /* calcDeltas = */ 1, tmp);

        REAL(dval)[i] = discFac * tmp[PEARSON0_PRICE_INDEX];
        del1 = REAL(d1)[i] = tmp[PEARSON0_D1_INDEX];
        del2 = REAL(d2)[i] = tmp[PEARSON0_D2_INDEX];

        dt[i] = i < T_len - 1 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len-1];
        dt[i] /= 365;

        if (dt[i] <= 0)
            error ("Unexpected non-positive time increment.");

        discFac *= exp(-REAL(r)[i] * dt[i]);

        pltoday  = -del1 * (REAL(S1)[i+1] - REAL(S1)[i]);
        pltoday -=  del1 * REAL(S1)[i] * expm1((REAL(q1)[i] - REAL(r)[i]) * dt[i]);
        pltoday -=  del2 * (REAL(S2)[i+1] - REAL(S2)[i]);
        pltoday -=  del2 * REAL(S2)[i] * expm1((REAL(q2)[i] - REAL(r)[i]) * dt[i]);

        REAL(dpl)[i+1] = discFac * pltoday;

        dt[i] *= 365;
    }

    pearson0(0, REAL(S1)[T_len], REAL(S2)[T_len], X, 
             /*vol1*/ 0, /*vol2*/ 0, /*rho*/ 0,
             /*r*/0, /*q1*/0, /*q2*/0, option_type, 
             /*grid_pts_pow*/0, /*calcDeltas*/ 0, tmp);
    REAL(dval)[T_len] = discFac * tmp[PEARSON0_PRICE_INDEX];
    REAL(d1)[T_len] = tmp[PEARSON0_D1_INDEX];
    REAL(d2)[T_len] = tmp[PEARSON0_D2_INDEX];

    for (i = 1; i < T_len + 1; i++) {
        if (dt[i-1] > max_dt)
            REAL(h1)[i] = REAL(h2)[i] = REAL(h3)[i] = REAL(o1)[i] = REAL(o2)[i] = REAL(o3)[i] = NA_REAL;
        else {

            den = sqrt(REAL(S1)[i-1] * REAL(S2)[i-1]);

            REAL(h1)[i] = den > 0 ? REAL(dpl)[i] / den : NA_REAL;
            REAL(h2)[i] = den > 0 ? REAL(h1)[i] / sqrt(dt[i-1]) : NA_REAL;
            REAL(o1)[i] = den > 0 ? (REAL(dval)[i] - REAL(dval)[i-1]) / den : NA_REAL;
            REAL(o2)[i] = den > 0 ? REAL(o1)[i] / sqrt(dt[i-1]) : NA_REAL;

            den = REAL(dval)[i-1];

            REAL(h3)[i] = REAL(dpl)[i] / den;
            REAL(o3)[i] = (REAL(dval)[i] - REAL(dval)[i-1]) / den;
        }
    }
    
    names = isNull(GET_NAMES(S1)) ? GET_NAMES(S2) : GET_NAMES(S1);

    SET_NAMES(dpl, names);
    SET_NAMES(dval, names);
    SET_NAMES(h1, names);
    SET_NAMES(o1, names);
    SET_NAMES(h2, names);
    SET_NAMES(o2, names);
    SET_NAMES(h3, names);
    SET_NAMES(o3, names);
    SET_NAMES(d1, names);
    SET_NAMES(d2, names);

    SET_ELT(rlist, i_hedge, dpl);
    SET_ELT(rlist, i_option, dval);
    SET_ELT(rlist, i_h1, h1);
    SET_ELT(rlist, i_o1, o1);
    SET_ELT(rlist, i_h2, h2);
    SET_ELT(rlist, i_o2, o2);
    SET_ELT(rlist, i_h3, h3);
    SET_ELT(rlist, i_o3, o3);
    SET_ELT(rlist, i_d1, d1);
    SET_ELT(rlist, i_d2, d2);
    SET_ELT(rlist, i_tau, T);
    SET_ELT(rlist, i_S1, S1);
    SET_ELT(rlist, i_S2, S2);
    SET_ELT(rlist, i_K, K);
    SET_ELT(rlist, i_vol1, vol1);
    SET_ELT(rlist, i_vol2, vol2);
    set_names(rlist, rlist_names);
    UNPROT2;
    return rlist;
}




static void listSprdOptDh3_chk(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, 
                               SEXP *vol1, SEXP *vol2, SEXP *rho, SEXP *r, 
                               SEXP *q1, SEXP *q2, SEXP *opt, SEXP *nGrdPts, 
                               int *numprot)
{
#define NARGS 11
    long i, N = 0;
    SEXP *args[NARGS] = {tau, S1, S2, K, vol1, vol2, rho, r, q1, q2, opt};

    for (i = 0; i < NARGS; i++) {
    
        if (!isNewList(*args[i]))
            PROT2(*args[i] = tolist(*args[i]), *numprot);

        N = MAX(N, length(*args[i]));
  
    }

    if (1 < N)
        for (i = 0; i < NARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = listrep(*args[i], N), *numprot);

  
#undef NARGS
}



SEXP listSprdOptDh3(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho, SEXP r, SEXP q1,
                    SEXP q2, SEXP opt, SEXP nGrdPts, SEXP maxdt)
{
    int numprot = 0;
    long i, N;

    SEXP ans, tmp;

    listSprdOptDh3_chk(&tau, &S1, &S2, &K, &vol1, &vol2, &rho,
                       &r, &q1, &q2, &opt, &nGrdPts,
                       &numprot);

    N = length(tau);
    PROT2(ans = NEW_LIST(N), numprot);
  
    for (i = 0; i < N; i++) {

        tmp = sprdOptDh3(GET_ELT(tau, i), GET_ELT(S1, i), GET_ELT(S2, i),
                         GET_ELT(K, i), GET_ELT(vol1, i), GET_ELT(vol2, i),
                         GET_ELT(rho, i), GET_ELT(r, i), GET_ELT(q1, i),
                         GET_ELT(q2, i), GET_ELT(opt, i), nGrdPts, maxdt);

        SET_ELT(ans, i, tmp);

    }

    SET_NAMES(ans, GET_NAMES(tau));
    UNPROT2;
    return ans;
}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following section contains functions similar to that of the previous
**  section.  The difference is they allow for the short rate and the
**  discount rate for the term of the option to be different.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/




static int sprdOptDh4_chk(SEXP *T, SEXP *S1, SEXP *S2,
                          SEXP *K, SEXP *vol1, SEXP *vol2,
                          SEXP *rho, SEXP *shortRate, SEXP *longRate, 
                          SEXP *q1, SEXP *q2, SEXP *opt,
                          int *numprot)
{
#define NNARGS 11
    int i, ok, allFinite;
    long t, T_len, S1_len, S2_len;
    SEXP *args[NNARGS+1] = {T, S1, S2, K, vol1, vol2,
                            rho, shortRate, longRate, q1, q2, opt};

    for (i = 0; i < NNARGS + 1; i++)
        if (0 == args[i])
            error ("Null pointer.");

    for (i = 0; i < NNARGS; i++)
        ENSURE_NUMERIC(*args[i], *numprot);

    ENSURE_CHAR(*opt, *numprot);

    T_len = length(*T);
    S1_len = length(*S1);
    S2_len = length(*S2);

    if (T_len != S1_len - 1 || S1_len != S2_len)
        error ("Length problems");

    for (i = 4; i < NNARGS; i++)
        if (1 == length(*args[i]))
            PROT2(*args[i] = numrep(*args[i], T_len), *numprot);

    ok = 1;
    for (i = 4; i < NNARGS; i++)
        ok *= (T_len == length(*args[i]));

    if (!ok)
        error ("More length problems");

    allFinite = 1;

    for (i = 0; i < NNARGS && allFinite; i++)
        for (t = 0; t < length(*args[i]) && allFinite; t++)
            allFinite *= R_FINITE(REAL(*args[i])[t]);

    return allFinite;
#undef NNARGS
}



/*
**  Description: Computes the daily delta-hedging P&L and option
**  values for a single spread option.
*/

SEXP sprdOptDh4(SEXP T, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho, 
                SEXP shortRate, SEXP longRate, SEXP q1, SEXP q2, SEXP opt, SEXP nGrdPts, SEXP maxdt)
{
    enum {i_hedge, i_option, i_o1, i_o2, i_o3, i_h1, i_h2, i_h3, i_d1, i_d2, 
          i_tau, i_S1, i_S2, i_K, i_vol1, i_vol2, i_rho, i_totalHedgesProfit, i_totalOptionProfit, ANS_LEN};

    const int max_dt = asInteger(maxdt);
    int numprot = 0;
    long i, T_len;
    double optnProfit, del1, del2, pltoday, discFac, den, *dt, tmp[PEARSON0_PTR_LEN], X;
    double totalHedgesProfit = 0, totalOptionProfit = 0;
    char option_type;
    char *rlist_names[ANS_LEN] = {"hedgeProfits", "optionValues", 
                                  "o1","o2","o3","h1","h2","h3",
                                  "d1", "d2", 
                                  "tau","S1","S2","K","vol1","vol2","rho", 
                                  "totalHedgesProfit", "totalOptionProfit"};

    SEXP dailyProfitHedges, sprdOptionVal;
    SEXP rlist;
    SEXP hedgesProfit1, optionProfit1, hedgesProfit2, optionProfit2, optionProfit3, hedgesProfit3;
    SEXP d1, d2, names;

    const long grdp = asInteger(nGrdPts);
    const SEXP *r = &shortRate;

    PROT2(rlist = NEW_LIST(ANS_LEN), numprot);

    SET_ELT(rlist, i_rho, rho);

    sprdOptDh4_chk(&T, &S1, &S2, &K, &vol1, &vol2,
                   &rho, &shortRate, &longRate, &q1, &q2, &opt,
                   &numprot);

    T_len = length(T);
    discFac = 1;

    PROT2(dailyProfitHedges = NEW_NUMERIC(T_len+1), numprot);
    PROT2(sprdOptionVal = NEW_NUMERIC(T_len+1), numprot);
    PROT2(hedgesProfit1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(hedgesProfit2 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(hedgesProfit3 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(optionProfit1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(optionProfit2 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(optionProfit3 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(d1 = NEW_NUMERIC(T_len+1), numprot);
    PROT2(d2 = NEW_NUMERIC(T_len+1), numprot);

    dt = (double *) R_alloc(T_len, DSIZE);

    REAL(dailyProfitHedges)[0] = 0;
    REAL(hedgesProfit1)[0] = REAL(optionProfit1)[0] = NA_REAL;
    REAL(hedgesProfit2)[0] = REAL(optionProfit2)[0] = NA_REAL;
    REAL(hedgesProfit3)[0] = REAL(optionProfit3)[0] = NA_REAL;

    option_type = *CHAR(GET_ELT(opt, 0));
    X = asReal(K);


    for (i = 0; i < T_len; i++) {

        pearson0(REAL(T)[i]/365, REAL(S1)[i], REAL(S2)[i],
                 X, REAL(vol1)[i], REAL(vol2)[i], 
                 REAL(rho)[i], /*r*/ REAL(longRate)[i], REAL(q1)[i], 
                 REAL(q2)[i], option_type,
                 grdp, /* calcDeltas */ 1, tmp);

        REAL(sprdOptionVal)[i] = discFac * tmp[PEARSON0_PRICE_INDEX];

        if (0 == i)
            totalOptionProfit = -tmp[PEARSON0_PRICE_INDEX];

        del1 = REAL(d1)[i] = tmp[PEARSON0_D1_INDEX];
        del2 = REAL(d2)[i] = tmp[PEARSON0_D2_INDEX];

        dt[i] = i < T_len - 1 ? REAL(T)[i] - REAL(T)[i+1] : REAL(T)[T_len-1];
        dt[i] /= 365;

        if (dt[i] <= 0)
            error ("Unexpected non-positive time increment.");

        discFac *= exp(-REAL(*r)[i] * dt[i]);

        pltoday  = -del1 * (REAL(S1)[i+1] - REAL(S1)[i]);
        pltoday -=  del1 * REAL(S1)[i] * expm1((REAL(q1)[i] - REAL(*r)[i]) * dt[i]);
        pltoday -=  del2 * (REAL(S2)[i+1] - REAL(S2)[i]);
        pltoday -=  del2 * REAL(S2)[i] * expm1((REAL(q2)[i] - REAL(*r)[i]) * dt[i]);

        REAL(dailyProfitHedges)[i+1] = discFac * pltoday;
        totalHedgesProfit += REAL(dailyProfitHedges)[i+1];

        dt[i] *= 365;
    }

    pearson0(0, REAL(S1)[T_len], REAL(S2)[T_len], X,
             /*vol1*/ 0, /*vol2*/ 0, /*rho*/ 0,
             /*r*/0, /*q1*/0, /*q2*/0, option_type, 
             /*grid_pts_pow*/0, /*calcDeltas*/ 0, tmp);

    totalOptionProfit += exp(-REAL(longRate)[0] * REAL(T)[0] / 365.0) * tmp[PEARSON0_PRICE_INDEX];

    REAL(sprdOptionVal)[T_len] = discFac * tmp[PEARSON0_PRICE_INDEX];
    REAL(d1)[T_len] = tmp[PEARSON0_D1_INDEX];
    REAL(d2)[T_len] = tmp[PEARSON0_D2_INDEX];

    for (i = 1; i < T_len + 1; i++) {

        if (dt[i-1] > max_dt) {

            REAL(hedgesProfit1)[i] = REAL(hedgesProfit2)[i] = REAL(hedgesProfit3)[i] = NA_REAL;
            REAL(optionProfit1)[i] = REAL(optionProfit2)[i] = REAL(optionProfit3)[i] = NA_REAL;

        } else {

            optnProfit = REAL(sprdOptionVal)[i] - REAL(sprdOptionVal)[i-1];

            den = sqrt(REAL(S1)[i-1] * REAL(S2)[i-1]);

            REAL(hedgesProfit1)[i] = den > 0 ? REAL(dailyProfitHedges)[i] / den : NA_REAL;
            REAL(hedgesProfit2)[i] = den > 0 ? REAL(hedgesProfit1)[i] / sqrt(dt[i-1]) : NA_REAL;
            REAL(optionProfit1)[i] = den > 0 ? optnProfit / den : NA_REAL;
            REAL(optionProfit2)[i] = den > 0 ? REAL(optionProfit1)[i] / sqrt(dt[i-1]) : NA_REAL;

            den = REAL(sprdOptionVal)[i-1];

            REAL(hedgesProfit3)[i] = den > 0 ? REAL(dailyProfitHedges)[i] / den : NA_REAL;
            REAL(optionProfit3)[i] = den > 0 ? optnProfit / den : NA_REAL;

        }
    }
    
    names = isNull(GET_NAMES(S1)) ? GET_NAMES(S2) : GET_NAMES(S1);

    SET_NAMES(dailyProfitHedges, names);
    SET_NAMES(sprdOptionVal, names);
    SET_NAMES(hedgesProfit1, names);
    SET_NAMES(optionProfit1, names);
    SET_NAMES(hedgesProfit2, names);
    SET_NAMES(optionProfit2, names);
    SET_NAMES(hedgesProfit3, names);
    SET_NAMES(optionProfit3, names);
    SET_NAMES(d1, names);
    SET_NAMES(d2, names);

    SET_ELT(rlist, i_hedge, dailyProfitHedges);
    SET_ELT(rlist, i_option, sprdOptionVal);
    SET_ELT(rlist, i_h1, hedgesProfit1);
    SET_ELT(rlist, i_o1, optionProfit1);
    SET_ELT(rlist, i_h2, hedgesProfit2);
    SET_ELT(rlist, i_o2, optionProfit2);
    SET_ELT(rlist, i_h3, hedgesProfit3);
    SET_ELT(rlist, i_o3, optionProfit3);
    SET_ELT(rlist, i_d1, d1);
    SET_ELT(rlist, i_d2, d2);
    SET_ELT(rlist, i_tau, T);
    SET_ELT(rlist, i_S1, S1);
    SET_ELT(rlist, i_S2, S2);
    SET_ELT(rlist, i_K, K);
    SET_ELT(rlist, i_vol1, vol1);
    SET_ELT(rlist, i_vol2, vol2);
    SET_ELT(rlist, i_totalHedgesProfit, ScalarReal(totalHedgesProfit));
    SET_ELT(rlist, i_totalOptionProfit, ScalarReal(totalOptionProfit));

    set_names(rlist, rlist_names);

    UNPROT2;

    return rlist;
}




static void listSprdOptDh4_chk(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, 
                               SEXP *vol1, SEXP *vol2, SEXP *rho, 
                               SEXP *shortRate, SEXP *longRate, 
                               SEXP *q1, SEXP *q2, SEXP *opt,
                               int *numprot)
{
#define NARGS 12
    long i, N = 0;
    SEXP *args[NARGS] = {tau, S1, S2, K, vol1, vol2, rho, shortRate, longRate, q1, q2, opt};

    for (i = 0; i < NARGS; i++)
        if (0 == args[i])
            error ("Null pointer.");

    for (i = 0; i < NARGS; i++) {
    
        if (!isNewList(*args[i]))
            PROT2(*args[i] = tolist(*args[i]), *numprot);

        N = MAX(N, length(*args[i]));
  
    }

    if (1 < N)
        for (i = 0; i < NARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = listrep(*args[i], N), *numprot);

  
#undef NARGS
}



SEXP listSprdOptDh4(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho, 
                    SEXP shortRate, SEXP longRate, 
                    SEXP q1, SEXP q2, SEXP opt, SEXP nGrdPts, SEXP maxdt)
{
    int numprot = 0;
    long i, N;

    SEXP ans, tmp;

    listSprdOptDh4_chk(&tau, &S1, &S2, &K, &vol1, &vol2, &rho,
                       &shortRate, &longRate, &q1, &q2, &opt,
                       &numprot);

    N = length(tau);
    PROT2(ans = NEW_LIST(N), numprot);
  
    for (i = 0; i < N; i++) {

        tmp = sprdOptDh4(GET_ELT(tau, i), GET_ELT(S1, i), GET_ELT(S2, i),
                         GET_ELT(K, i), GET_ELT(vol1, i), GET_ELT(vol2, i),
                         GET_ELT(rho, i), GET_ELT(shortRate, i), GET_ELT(longRate, i),
                         GET_ELT(q1, i), GET_ELT(q2, i), GET_ELT(opt, i), nGrdPts, maxdt);

        SET_ELT(ans, i, tmp);

    }

    SET_NAMES(ans, GET_NAMES(tau));
    UNPROT2;
    return ans;
}

