#include "multi.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains functions for computing implied correlations.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

static double pearsonErr(const double T, const double S1, const double S2,
                         const double K, const double vol1, const double vol2,
                         const double rho, const double r, const double q1,
                         const double q2, const char opt, const int nGrdPts,
                         const double mprice)

/*******************************************************************
 *
 *  Description: Computes the Pearson spread options price x for the
 *    given parameters and returns log(x / mprice).
 *
 *  Caution: T is in years (or whatever units vol, r and q
 *    are in).
 *
 *******************************************************************/

{
    double val[PEARSON0_PTR_LEN] = {0.0};

    pearson0(T, S1, S2, K, vol1, vol2, rho, r, q1, q2, opt, nGrdPts,
             /* calcDeltas = */ 0, val);

    return log(val[PEARSON0_PRICE_INDEX] / mprice);
}

#define MAXIT (100)
static void getBounds(const double T, const double S1, const double S2,
                      const double K, const double vol1, const double vol2,
                      const double mprice, const double r, const double q1,
                      const double q2, const char optype, const int nGrdPts,
                      const double tol, int maxit, double impcorr,
                      double *lower, double *upper)
{
    double err = 0, errUpper = 0, errLower = 0, upperOld, width = 0.1;
    long   it = 0;

    if (!upper || !lower) error("null pointer where not allowed in getBounds");

    if (!R_FINITE(impcorr)) return;

    *upper = *lower = impcorr;

    err = pearsonErr(T, S1, S2, K, vol1, vol2, impcorr, r, q1, q2, optype,
                     nGrdPts, mprice);

    GetRNGstate();

    if (err > 0)
    { /* initial correlation too low */

        do
        {
            it       += 1;
            upperOld  = *upper;
            *upper   += width * unif_rand();

            errUpper = pearsonErr(T, S1, S2, K, vol1, vol2, *upper, r, q1, q2,
                                  optype, nGrdPts, mprice);

            if (!R_FINITE(errUpper))
            {
                errUpper  = 1;
                *upper    = upperOld;
                width    /= 1 + unif_rand();
            } else
                width *= 1 + unif_rand();

        } while (errUpper > 0 && it < MAXIT);

    } else if (err < 0)
    { /* initial correlation too high */

        do
        {
            it     += 1;
            *lower += lower[0] < -1 ? runif(0, lower[0])
                                    : -(width *= runif(1, 2)) * unif_rand();

            errLower = pearsonErr(T, S1, S2, K, vol1, vol2, *lower, r, q1, q2,
                                  optype, nGrdPts, mprice);

        } while (errLower < 0 && it < MAXIT && R_FINITE(errLower));

    } else
    {
        if (0 != err) *lower = *upper = NA_REAL;

        return;
    }

    if (it >= MAXIT || !R_FINITE(errLower)) *lower = *upper = NA_REAL;

    PutRNGstate();
}

static void impCorr0(const double T, const double S1, const double S2,
                     const double K, const double vol1, const double vol2,
                     const double mprice, const double r, const double q1,
                     const double q2, const char optype, const long nGrdPts,
                     const double tol, int maxit, double *impcorr)

/*******************************************************************
 *
 *  Description: Computes the implied correlation and assigns it
 *    to *impcorr.
 *
 *  Caution: T should be in the same time units as r and q
 *    (usually years).
 *
 *******************************************************************/

{
    int    it    = 0;
    double lower = 0, upper = 0, err, errUpper, errLower;
    double lowerPrice            = 0.0;
    double upperPrice            = 0.0;
    double tmp[PEARSON0_PTR_LEN] = {0.0};

    /* Do some argument checking */
    if (!impcorr) error("Null pointer passed where not allowed");

    if (!R_FINITE(*impcorr)) *impcorr = 0.9;

    if (ABS(*impcorr) > 1) *impcorr = 0.9 * SIGN_OF(*impcorr);

    QUIT_IF(!R_FINITE(mprice), *impcorr, NA_REAL);
    QUIT_IF(mprice <= 0, *impcorr, NA_REAL);

    pearson0(T, S1, S2, K, vol1, vol2,
             /* rho = */ 0, r, q1, q2, optype, nGrdPts,
             /* calcDeltas = */ 0, tmp);

    QUIT_IF(!R_FINITE(tmp[PEARSON0_PRICE_INDEX]), *impcorr, NA_REAL);

    maxit = MAX(maxit, 20);

    /* Get upper and lower bounds for bisection */
    getBounds(T, S1, S2, K, vol1, vol2, mprice, r, q1, q2, optype, nGrdPts, tol,
              maxit, *impcorr, &lower, &upper);

    QUIT_IF(!R_FINITE(upper) || !R_FINITE(lower), *impcorr, NA_REAL);

    pearson0(T, S1, S2, K, vol1, vol2,
             /* rho = */ lower, r, q1, q2, optype, nGrdPts,
             /* calcDeltas = */ 0, tmp);

    upperPrice = tmp[PEARSON0_PRICE_INDEX];

    pearson0(T, S1, S2, K, vol1, vol2,
             /* rho = */ upper, r, q1, q2, optype, nGrdPts,
             /* calcDeltas = */ 0, tmp);

    lowerPrice = tmp[PEARSON0_PRICE_INDEX];

    QUIT_IF(lowerPrice > mprice || mprice > upperPrice, *impcorr, NA_REAL);
    QUIT_IF(!R_FINITE(lowerPrice) || !R_FINITE(upperPrice), *impcorr, NA_REAL);

    /* bisection */
    it = 0;

    do
    {
        *impcorr = 0.5 * (lower + upper);
        err = pearsonErr(T, S1, S2, K, vol1, vol2, *impcorr, r, q1, q2, optype,
                         nGrdPts, mprice);
        errLower = pearsonErr(T, S1, S2, K, vol1, vol2, lower, r, q1, q2,
                              optype, nGrdPts, mprice);
        if (SIGN_OF(err) == SIGN_OF(errLower))
            lower = *impcorr;
        else
        {
            errUpper = pearsonErr(T, S1, S2, K, vol1, vol2, upper, r, q1, q2,
                                  optype, nGrdPts, mprice);
            if (SIGN_OF(err) == SIGN_OF(errUpper))
                upper = *impcorr;
            else
                break;
        }
    } while (ABS(err) > tol && ++it < maxit);

    if (ABS(err) > tol) *impcorr = NA_REAL;
}
#undef MAXIT

static void impCorr_chkargs(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, SEXP *vol1,
                            SEXP *vol2, SEXP *mprice, SEXP *r, SEXP *q1,
                            SEXP *q2, SEXP *op_type, const long grdPts,
                            SEXP *initGuess, int *numprot)
{
#define N_NUMARGS 10
    int i, N = 0;

    SEXP *args[N_NUMARGS + 2] = {tau,    S1, S2, K,  vol1,    vol2,
                                 mprice, r,  q1, q2, op_type, initGuess};

    char *names[N_NUMARGS + 2] = {"tau",  "S1",   "S2",      "K",
                                  "vol1", "vol2", "mprice",  "r",
                                  "q1",   "q2",   "op_type", "initGuess"};

    if (grdPts < 0) error("Negative number of grid pts.");

    /* Coerce types  */
    for (i = 0; i < N_NUMARGS; i++)
        PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (!length(*initGuess) || R_NilValue == *initGuess)
        PROT2(*initGuess = NEW_NUMERIC(1), *numprot);

    if (!IS_NUMERIC(*initGuess))
        PROT2(*initGuess = AS_NUMERIC(*initGuess), *numprot);

    /* Find maximal length */
    for (i = 0, N = 0; i < N_NUMARGS + 2; i++) N = MAX(N, length(*args[i]));

    if (N > 1)
    { /* Repeat arguments of length == 1. */

        for (i = 0; i < N_NUMARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = numrep(*args[i], N), *numprot);

        if (1 == length(*op_type))
            PROT2(*op_type = charrep(*op_type, N), *numprot);

        if (1 == length(*initGuess))
            PROT2(*initGuess = numrep(*initGuess, N), *numprot);
    }

    /* Check lengths */

    for (i = 0; i < N_NUMARGS + 2; i++)
        if (N != length(*args[i]))
            error("Argument %s of wrong length", names[i]);

#undef N_NUMARGS
}

/*******************************************************************
 *
 *  Description: Computes the implied correlation of a option.
 *
 *  Caution: Assumes tau is in days.
 *
 *******************************************************************/

SEXP impCorr(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2,
             SEXP mprice, SEXP r, SEXP q1, SEXP q2, SEXP op_type, SEXP nGrdPts,
             SEXP initGuess, SEXP tol, SEXP maxit)
{
    int        numprot = 0;
    long       N, i;
    const long grdPts = asInteger(nGrdPts);
    SEXP       ans;

    impCorr_chkargs(&tau, &S1, &S2, &K, &vol1, &vol2, &mprice, &r, &q1, &q2,
                    &op_type, grdPts, &initGuess, &numprot);

    /* By now tau, S, K, mprice, r, q, op_type should all have same length */

    N = length(tau);

    PROT2(ans = NEW_NUMERIC(N), numprot);
    memcpy(REAL(ans), REAL(initGuess), N * sizeof(double));

    for (i = 0; i < N; i++)
        impCorr0(REAL(tau)[i] / 365.0, REAL(S1)[i], REAL(S2)[i], REAL(K)[i],
                 REAL(vol1)[i], REAL(vol2)[i], REAL(mprice)[i], REAL(r)[i],
                 REAL(q1)[i], REAL(q2)[i], *CHAR(GET_ELT(op_type, i)),
                 /* nGrdPts */ 0, asReal(tol), asInteger(maxit), REAL(ans) + i);

    if (grdPts > 0)
        for (i = 0; i < N; i++)
            impCorr0(REAL(tau)[i] / 365.0, REAL(S1)[i], REAL(S2)[i], REAL(K)[i],
                     REAL(vol1)[i], REAL(vol2)[i], REAL(mprice)[i], REAL(r)[i],
                     REAL(q1)[i], REAL(q2)[i], *CHAR(GET_ELT(op_type, i)),
                     grdPts, asReal(tol), asInteger(maxit), REAL(ans) + i);

    UNPROTECT(numprot);

    return ans;
}
