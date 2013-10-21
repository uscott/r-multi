#include "uri.h"
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains functions for calculating the Margrabe
**  spread option formula and for evaluating a spread option with
**  nonzero strike using simulation with stratified stampling
**  and the control variate technique.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



/*
**  NB: ans must point to contiguous memory block of length
**  sizeof(double) * MARGRABE0_PTR_LEN.
*/
void margrabe0(const double T, const double S1, 
               const double S2, const double vol1, const double vol2,
               const double rho, const double r, const double q1,   
               const double q2, const char opt, double *ans)
{
    register int ok, i;
    register double d1, d2, Nd1, Nd2, div1, div2, Npd1, Npd2, vol, v;
    double c[MARGRABE0_PTR_LEN] = {0.0};
    double p[MARGRABE0_PTR_LEN] = {0.0};
    const int EXPIRATION = (0.0 == T);

    if (!ans)
        error ("Illegal null pointer passed");
	
    for (i = 0; i < MARGRABE0_PTR_LEN; i++)
        ans[i] = NA_REAL;

    ok = T >= 0 && S1 >= 0 && S2 >= 0 && R_FINITE(q1) &&
        vol1 >= 0 && vol2 >= 0 && R_FINITE(T) && R_FINITE(S1) && 
        R_FINITE(S2) && R_FINITE(vol1) && R_FINITE(r) && R_FINITE(vol2) &&
        R_FINITE(q2) && R_FINITE(rho);

    if (!ok)
        return;

    switch(opt) {

    case 'c':

        if (!EXPIRATION) {

            vol = sqrt(SQR(vol1) + SQR(vol2) - 2.0 * rho * vol1 * vol2);
            d1  = (log(S2/S1) + (q1 - q2 + 0.5 * SQR(vol)) * T) / (vol * sqrt(T));
            d2  = d1 - vol * sqrt(T);
            Nd1  = NORM_DIST(d1);
            Nd2  = NORM_DIST(d2);
            Npd1 = NORM_DENS(d1);
            Npd2 = NORM_DENS(d2);
            div1 = exp(-q1 * T);
            div2 = exp(-q2 * T);
        
            v = S2 * sqrt(T) * Npd1 * div2 / vol; /* Used in vega1, vega2 and eta */

        }

        ans[MARGRABE0_PRICE_INDEX] =  EXPIRATION ? DMAX(S2 - S1, 0) : div2 * S2 * Nd1 - div1 * S1 * Nd2;
        ans[MARGRABE0_D1_INDEX]    =  EXPIRATION ?  STEP_FUN(S2 - S1) : -div1 * Nd2;
        ans[MARGRABE0_D2_INDEX]    =  EXPIRATION ? -STEP_FUN(S1 - S2) :  div2 * Nd1;     
        ans[MARGRABE0_G1_INDEX]    =  EXPIRATION ? (S1 == S2 ? NA_REAL : 0) : div1 * Npd2 / (S1 * vol * sqrt(T));            
        ans[MARGRABE0_G2_INDEX]    =  EXPIRATION ? (S1 == S2 ? NA_REAL : 0) : div2 * Npd1 / (S2 * vol * sqrt(T));            
        ans[MARGRABE0_V1_INDEX]    =  EXPIRATION ? 0 : v * (vol1 - rho * vol2);
        ans[MARGRABE0_V2_INDEX]    =  EXPIRATION ? 0 : v * (vol2 - rho * vol1);          
        ans[MARGRABE0_ETA_INDEX]   =  EXPIRATION ? 0 : -v * vol2 * vol1;

        break;

    case 'p':

        margrabe0(T, S2, S1, vol2, vol1, rho, r, q2, q1, 'c', ans);
        dswap(ans + MARGRABE0_D1_INDEX, ans + MARGRABE0_D2_INDEX);
        dswap(ans + MARGRABE0_G1_INDEX, ans + MARGRABE0_G2_INDEX);
        dswap(ans + MARGRABE0_V1_INDEX, ans + MARGRABE0_V2_INDEX);

        break;

    case 's':

        margrabe0(T, S1, S2, vol1, vol2, rho, r, q1, q2, 'c', c);
        margrabe0(T, S1, S2, vol1, vol2, rho, r, q1, q2, 'p', p);

        for (i = 0; i < MARGRABE0_PTR_LEN; i++)
            ans[i] = c[i] + p[i];

        break;

    default:
        break;

    }

}


static void margrabe_chkargs(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *vol1,
                             SEXP *vol2, SEXP *rho, SEXP *r, SEXP *q1,
                             SEXP *q2, SEXP *opt, int  *numprot)
{
#define N_NUMARGS 9
    int i, N;
    SEXP *args [N_NUMARGS + 1] = {tau, S1, S2, vol1, vol2, rho, r, q1, q2, opt};
    char *names[N_NUMARGS + 1] = {"tau", "S1", "S2", "vol1", "vol2", "rho", "r", "q1", "q2", "opt"};
    /* Coerce types  */
    for (i = 0; i < N_NUMARGS; i++)
        PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);
    /* Find maximal length */
    N = 0;
    for (i = 0; i < N_NUMARGS + 1; i++)
        N = MAX(N, length(*args[i]));
    /* Repeat arguments of length == 1. */
    if (N > 1) {
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
#undef N_NUMARGS
}




/*
**  Description: Returns SEXP of type matrix with columns
**  corresponding to Margrabe formula price, deltas,
**  gammas & vegas.
**  Caution: Assumes tau is in *calender days*.
*/

SEXP margrabe(SEXP tau, SEXP S1, SEXP S2, SEXP vol1, SEXP vol2, SEXP rho,
              SEXP r, SEXP q1, SEXP q2, SEXP opt)
{
    int i, N, numprot = 0;
    SEXP ans, dimnames, colnames;
    carray ans_matrix;
    char *names[MARGRABE0_PTR_LEN] = MARGRABE0_RET_NAMES;

    margrabe_chkargs(&tau, &S1, &S2, &vol1, &vol2, &rho, &r, &q1, &q2, &opt, &numprot);

    N = length(tau); /* By now all arguments should have same length */

    ans_matrix = make_zero_matrix(N, MARGRABE0_PTR_LEN);
  
    for (i = 0; i < N; i++)
        margrabe0(REAL(tau)[i] / 365.0, REAL(S1)[i], REAL(S2)[i], 
                  REAL(vol1)[i], REAL(vol2)[i], REAL(rho)[i],
                  REAL(r)[i], REAL(q1)[i], REAL(q2)[i], 
                  *CHAR(GET_ELT(opt, i)), ARRAY2(ans_matrix)[i]);

    PROT2(ans = carray_to_sexp(ans_matrix), numprot);
    PROT2(colnames = NEW_STRING(MARGRABE0_PTR_LEN), numprot);
    PROT2(dimnames = NEW_LIST(2), numprot);

    for (i = 0; i < MARGRABE0_PTR_LEN; i++)
        SET_ELT(colnames, i, mkChar(names[i]));

    SET_ELT(dimnames, 0, R_NilValue);
    SET_ELT(dimnames, 1, colnames);

    setAttrib(ans, R_DimNamesSymbol, dimnames);

    UNPROTECT(numprot);
    return ans;
}








static void margrabe2_0(const double T, const double S1, const double S2,
                        const double K, const double vol1, const double vol2,
                        const double rho, const double r, const double q1,   
                        const double q2, const char opt, const int npaths,
                        double *e1, double *e2, double *ans)
{
    register int ok, i;
    register double discFac;
    const double u = 1.0 + 1e-4;
    register double price, price1, price2;
    register double x1, x2, mu1, mu2, v1, v2, z;
    double tmp[MAX(MARGRABE0_PTR_LEN, MARGRABE2_0_PTR_LEN)] = {0.0};

    if (!ans || !e1 || !e2)
        error ("Illegal null pointer passed");
	
    ans[0] = ans[1] = ans[2] = NA_REAL;

    ok = T >= 0 && S1 >= 0 && S2 >= 0 && R_FINITE(q1) &&
        vol1 >= 0 && vol2 >= 0 && R_FINITE(T) && R_FINITE(S1) && 
        R_FINITE(S2) && R_FINITE(vol1) && R_FINITE(r) && R_FINITE(vol2) &&
        R_FINITE(q2) && R_FINITE(K) && npaths > 0;

    if (!ok)
        return;

    if (0.0 == K) {
        margrabe0(T, S1, S2, vol1, vol2, rho, r, q1, q2, opt, tmp);
        ans[0] = tmp[MARGRABE0_PRICE_INDEX];
        ans[1] = tmp[MARGRABE0_D1_INDEX];
        ans[2] = tmp[MARGRABE0_D2_INDEX];
        return;
    }

    v1  = vol1 * sqrt(T); v2  = vol2 * sqrt(T);
    mu1 = exp((r - q1) * T); mu2 = exp((r - q2) * T);
    discFac  = exp(- r * T);
    price = price1 = price2 = 0.0;

    switch(opt) {

    case 'c':   
        for (i = 0; i < npaths; i++, e1++, e2++) {

            x1 = S1 * mu1 * exp(*e1 * v1 - 0.5 * SQR(v1));
            z  = rho * (*e1) + sqrt(1.0 - SQR(rho)) * (*e2);
            x2 = S2 * mu2 * exp( z  * v2 - 0.5 * SQR(v2));

            price  += MAX(x2 - x1 - K, 0) - MAX(x2 - x1, 0);
            price1 += MAX(x2 - u * x1 - K, 0) - MAX(x2 - u * x1, 0);
            price2 += MAX(u * x2 - x1 - K, 0) - MAX(u * x2 - x1, 0);
        }
        break;
    case 'p':
        for (i = 0; i < npaths; i++, e1++, e2++) {
            x1 = S1 * mu1 * exp(*e1 * v1 - 0.5 * SQR(v1));
            z  = rho * (*e1) + sqrt(1.0 - SQR(rho)) * (*e2);
            x2 = S2 * mu2 * exp(z * v2 - 0.5 * SQR(v2));

            price  += MAX(K - (x2 - x1), 0) - MAX(x1 - x2, 0);
            price1 += MAX(K - (x2 - u * x1), 0) - MAX(u * x1 - x2, 0);
            price2 += MAX(K - (u * x2 - x1), 0) - MAX(x1 - u * x2, 0);
        }
        break;
    case 's':
        for (i = 0; i < npaths; i++, e1++, e2++) {
            x1 = S1 * mu1 * exp(*e1 * v1 - 0.5 * SQR(v1));
            z  = rho * (*e1) + sqrt(1.0 - SQR(rho)) * (*e2);
            x2 = S2 * mu2 * exp( z  * v2 - 0.5 * SQR(v2));

            price  += ABS(x2 - x1 - K) - ABS(x2 - x1);
            price1 += ABS(x2 - u * x1 - K) - ABS(x2 - u * x1);
            price2 += ABS(u * x2 - x1 - K) - ABS(u * x2 - x1);
        }
        break;
    default:
        break;
    }

    price  = discFac * price  / (double)npaths;
    price1 = discFac * price1 / (double)npaths;
    price2 = discFac * price2 / (double)npaths;
 
    margrabe0(T, S1, S2, vol1, vol2, rho, r, q1, q2, opt, tmp);
    price = price + tmp[0];
    margrabe0(T, u * S1, S2, vol1, vol2, rho, r, q1, q2, opt, tmp);
    price1 = price1 + tmp[0];
    margrabe0(T, S1, u * S2, vol1, vol2, rho, r, q1, q2, opt, tmp);
    price2 = price2 + tmp[0];

    ans[0] = price;                                               
    ans[1] = (price1 - price) / ((u - 1.0) * S1);                   
    ans[2] = (price2 - price) / ((u - 1.0) * S2);                 

}

    

static void margrabe2_chkargs(SEXP *tau, SEXP *S1, 
                              SEXP *S2, SEXP *K, SEXP *vol1,
                              SEXP *vol2, SEXP *rho, SEXP *r,
                              SEXP *q1, SEXP *q2, SEXP *opt,
                              SEXP *npaths, SEXP *eps, int *numprot)
{
#define N_NUMARGS 10
    int i, N = 0, n;
    SEXP *args [N_NUMARGS + 1] = {tau, S1, S2, K, vol1, vol2, rho, r, q1, q2, opt};
    char *names[N_NUMARGS + 1] = {"tau", "S1",   "S2", "K",
                                  "vol1", "vol2", "rho",
                                  "r", "q1", "q2", "opt"};

    /* Coerce types  */
    for (i = 0; i < N_NUMARGS; i++)
        PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (isNull(*eps) || !IS_MATRIX(*eps) || 2 != ncols(*eps) || 0 == nrows(*eps)) {
        n = (int) ceil(sqrt(asInteger(*npaths))) + 1;
        *eps = PROTECT(bvt_normuri(ScalarInteger(n), ScalarReal(0.0), ScalarInteger(1)));
        *numprot += 1;
    }

    /* Find maximal length */
    for (i = 0, N = 0; i < N_NUMARGS + 1; i++)
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

    if (!IS_MATRIX(*eps) || 2 != ncols(*eps) || !nrows(*eps))
        error ("eps is messed up");
#undef N_NUMARGS
}



SEXP margrabe2(SEXP tau, SEXP S1, SEXP S2, SEXP K,
               SEXP vol1, SEXP vol2, SEXP rho,
               SEXP r, SEXP q1, SEXP q2,
               SEXP opt, SEXP npaths, SEXP eps)
{
    int i, N, numprot = 0;
    SEXP ans, dimnames, colnames;
    carray ans_matrix;
    double *e1, *e2;
    char *names[MARGRABE2_0_PTR_LEN] = {"p", "d1", "d2"};

    margrabe2_chkargs(&tau, &S1, &S2, &K, &vol1, &vol2, 
                      &rho, &r, &q1, &q2, &opt, &npaths, &eps, &numprot);

  
    N = length(tau); /* By now all arguments should have same length */

    ans_matrix = make_zero_matrix(N, MARGRABE2_0_PTR_LEN);
  
    e1  = matcol1(eps, 0);
    e2  = matcol1(eps, 1);

    for (i = 0; i < N; i++)
        margrabe2_0(REAL(tau)[i] / 365.0, REAL(S1)[i], REAL(S2)[i],
                    REAL(K)[i], REAL(vol1)[i], REAL(vol2)[i],
                    REAL(rho)[i], REAL(r)[i], REAL(q1)[i],
                    REAL(q2)[i], *CHAR(GET_ELT(opt, i)), nrows(eps),
                    e1, e2, ARRAY2(ans_matrix)[i]);

    PROT2(ans = carray_to_sexp(ans_matrix),  numprot);
    PROT2(colnames = NEW_STRING(MARGRABE2_0_PTR_LEN), numprot);
    PROT2(dimnames = NEW_LIST(2), numprot);

    for (i = 0; i < MARGRABE2_0_PTR_LEN; i++)
        SET_ELT(colnames, i, mkChar(names[i]));

    SET_ELT(dimnames, 0, R_NilValue);
    SET_ELT(dimnames, 1, colnames);

    setAttrib(ans, R_DimNamesSymbol, dimnames);

    UNPROTECT(numprot);
    return ans;
}

