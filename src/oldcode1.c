#ifdef YOU_SUCK
#include "uri.h"

static void oneStepGrdSprd_old(long N, double S10, double S20, double sig1,
                               double sig2, double rho, double r, double q1,
                               double q2, double X, double T, double *value,
                               double h)
{
    double  omega, A;
    double *mat1, *mat2, *z, *n1, *n2, *art, *G1, *G2;
    double  m, diffExp, siggy, term1, term2;
    double  a_0, b_0, a_N, b_N, wind;
    long    i, N2 = 2 * N;

    mat1 = (double *)R_alloc(N2, sizeof(double));
    mat2 = (double *)R_alloc(N2, sizeof(double));
    G1   = (double *)R_alloc(N2, sizeof(double));
    G2   = (double *)R_alloc(N2, sizeof(double));
    z    = (double *)R_alloc(N2 + 1, sizeof(double));
    n1   = (double *)R_alloc(N2 + 1, sizeof(double));
    n2   = (double *)R_alloc(N2 + 1, sizeof(double));
    art  = (double *)R_alloc(N2 + 1, sizeof(double));

    omega = sig1 * sqrt(T) / h;
    siggy = sig1 * sqrt(T);

    for (i = 0; i <= N; i++)
    { /* rows and columns of underlying Toeplitz matrix */

        wind  = (i - N) / omega;
        n1[i] = NDist(wind - siggy);
        n2[i] = NDist(wind);
    }

    for (i = N + 1; i <= N2; i++)
    {
        n1[i] = NDist((i - N) / omega - siggy);
        n2[i] = 1 - n2[N2 - i];
    }

    for (i = 0; i < N2; i++)
    {
        mat1[i] = n1[i + 1] - n1[i];
        mat2[i] = n2[i + 1] - n2[i];
    }

    m = log(S10) + (r - q1 - 0.5 * sig1 * sig1) * T;

    for (i = 0; i <= N2; i++)
    {
        art[i] = exp(m + h * (i - N)); /* grid points */
        z[i]   = Inner(art[i], S10, S20, sig1, sig2, rho, r, q1, q2, T, X);
    }

    for (i = 0; i < N2; i++)
    { /* linear interpolation between grid points */

        diffExp = art[i + 1] - art[i];
        G1[i]   = (z[i + 1] - z[i]) / diffExp;
        G2[i]   = (art[i + 1] * z[i] - art[i] * z[i + 1]) / diffExp;
    }

    A = (r * (1.0 - rho * sig2 / sig1) - (q2 - q1 * rho * sig2 / sig1) +
         0.5 * rho * sig2 * (sig1 - rho * sig2)) *
        T;
    a_N = exp(A) * S20 / pow(S10, rho * sig2 / sig1);
    b_N = z[N2] - a_N * pow(art[N2], rho * sig2 / sig1);
    a_0 = exp(A) * S20 / pow(S10, rho * sig2 / sig1);
    b_0 = z[0] - a_0 * pow(art[0], rho * sig2 / sig1);

    *value = 0.0;

    for (i = 0; i < N2; i++) *value += mat1[i] * G1[i];

    *value *= exp(m + 0.5 * sig1 * sig1 * T);

    for (i = 0; i < N2; i++) *value += mat2[i] * G2[i];

    term1 = a_0 * NDist(-N / omega - rho * sig2 * sqrt(T)) *
            exp(rho * sig2 * m / sig1 + 0.5 * rho * rho * sig2 * sig2 * T);
    term2   = b_0 * NDist(-N / omega);
    *value += term1 + term2;

    term1 = a_N * (1.0 - NDist(N / omega - rho * sig2 * sqrt(T))) *
            exp(rho * sig2 * m / sig1 + 0.5 * rho * rho * sig2 * sig2 * T);
    term2   = b_N * (1.0 - NDist(N / omega));
    *value += term1 + term2;

    *value *= exp(-r * T);
}

void delta_to_strike(const int *n, const double *T, const double *S,
                     const double *vol, const double *del, const char **opt,
                     const double *r, double *q, double *ans)
{
    double const   *end = del + *n;
    register double t1, t2, tau;

    for (; del < end; vol++, S++, T++, del++, ans++, opt++)
    {
        if (*T < 0 || *S < 0 || *vol < 0)
        {
            *ans = NA_REAL;
            continue;
        }

        tau = *T / 365.0;

        switch (**opt)
        {
            case 'c':

                t1 = -*vol * sqrt(tau) *
                     qnorm(exp(*q * tau) * (*del), 0, 1, 1, 0);
                t2   = (*r - *q + 0.5 * SQR(*vol)) * tau;
                *ans = *S * exp(t1 + t2);
                break;

            case 'p':

                t1 = -*vol * sqrt(tau) *
                     qnorm(1.0 + exp(*q * tau) * (*del), 0, 1, 1, 0);
                t2   = (*r - *q + 0.5 * SQR(*vol)) * tau;
                *ans = *S * exp(t1 + t2);
                break;

            case 's':

                t1 = -*vol * sqrt(tau) *
                     qnorm(0.5 + 0.5 * exp(*q * tau) * (*del), 0, 1, 1, 0);
                t2   = (*r - *q + 0.5 * SQR(*vol)) * tau;
                *ans = *S * exp(t1 + t2);
                break;

            default:

                *ans = NA_REAL;
                break;
        }
    }
}

#define ANS_LEN 12
SEXP dvHedge2_debug(SEXP posn, SEXP posnStrikes, SEXP optionTypes,
                    SEXP hedgePeriod, SEXP T, SEXP S, SEXP vols,
                    SEXP volStrikes, SEXP r, SEXP q, SEXP undTransCosts,
                    SEXP optTransCosts, SEXP isRelUndCosts, SEXP isRelOptCosts)
{
    int numprot = 0;

    int  numVolCols, posnLen;
    long T_len, S_len, i;
    SEXP strikesTmp, ans, tmp, undPL, strdPL, undTC, strdTC, tmp2, maxUndPosn,
        minUndPosn, maxStrdPosn, minStrdPosn, maxDelta, minDelta, maxVega,
        minVega, r0, q0;

    char *names[ANS_LEN] = {"undPL",       "strdPL",      "undTC",
                            "strdTC",      "maxUndPosn",  "minUndPosn",
                            "maxStrdPosn", "minStrdPosn", "maxDelta",
                            "minDelta",    "maxVega",     "minVega"};

    dvHedge2_chk(&posn, &posnStrikes, &optionTypes, &hedgePeriod, &T, &S, &vols,
                 &volStrikes, &r, &r, &undTransCosts, &optTransCosts,
                 &isRelUndCosts, &isRelOptCosts, &numprot);

    S_len      = length(S);
    T_len      = S_len - 1;
    posnLen    = length(posn);
    numVolCols = ncols(vols);

    PROT2(ans = NEW_LIST(ANS_LEN), numprot);
    PROT2(strikesTmp = NEW_NUMERIC(posnLen), numprot);
    PROT2(undPL = NEW_NUMERIC(T_len), numprot);
    PROT2(strdPL = NEW_NUMERIC(T_len), numprot);
    PROT2(undTC = NEW_NUMERIC(T_len), numprot);
    PROT2(strdTC = NEW_NUMERIC(T_len), numprot);
    PROT2(r0 = NEW_NUMERIC(1), numprot);
    PROT2(q0 = NEW_NUMERIC(1), numprot);

    PROT2(maxUndPosn = NEW_NUMERIC(T_len), numprot);
    PROT2(minUndPosn = NEW_NUMERIC(T_len), numprot);
    PROT2(maxStrdPosn = NEW_NUMERIC(T_len), numprot);
    PROT2(minStrdPosn = NEW_NUMERIC(T_len), numprot);

    PROT2(maxDelta = NEW_NUMERIC(T_len), numprot);
    PROT2(minDelta = NEW_NUMERIC(T_len), numprot);
    PROT2(maxVega = NEW_NUMERIC(T_len), numprot);
    PROT2(minVega = NEW_NUMERIC(T_len), numprot);

    set_names(ans, names);

    dsetna(REAL(undPL), T_len);
    dsetna(REAL(strdPL), T_len);
    dsetna(REAL(undTC), T_len);
    dsetna(REAL(strdTC), T_len);

    for (i = 0; i < T_len; i++)
    {
        dmatrow1(posnStrikes, i, REAL(strikesTmp));

        REAL(r0)[0] = REAL(r)[i];
        REAL(q0)[0] = REAL(q)[i];

        tmp = dvHedge(
            posn, strikesTmp, optionTypes, hedgePeriod,
            PROTECT(dSubVector(T, i, T_len - 1)),
            PROTECT(dSubVector(S, i, S_len - 1)),
            PROTECT(dSubMatrix(vols, i, T_len - 1, 0, numVolCols - 1)),
            PROTECT(dSubMatrix(volStrikes, i, T_len - 1, 0, numVolCols - 1)),
            r0, q0, undTransCosts, optTransCosts, isRelUndCosts, isRelOptCosts);

        PROTECT(tmp);

        REAL(undPL)[i]  = *REAL(getListElt(tmp, "undPL"));
        REAL(strdPL)[i] = *REAL(getListElt(tmp, "strdPL"));
        REAL(undTC)[i]  = *REAL(getListElt(tmp, "undTC"));
        REAL(strdTC)[i] = *REAL(getListElt(tmp, "strdTC"));

        /* Following output mainly for debuggin purposes */
        REAL(maxUndPosn)
        [i] = double_max(REAL(getListElt(tmp, "undPosn")), S_len - i);
        REAL(minUndPosn)
        [i] = double_min(REAL(getListElt(tmp, "undPosn")), S_len - i);
        REAL(maxStrdPosn)
        [i] = double_max(REAL(getListElt(tmp, "strdPosn")), T_len - i);
        REAL(minStrdPosn)
        [i] = double_min(REAL(getListElt(tmp, "strdPosn")), T_len - i);

        REAL(maxDelta)
        [i] = double_max(REAL(getListElt(tmp, "deltas")), S_len - i);
        REAL(minDelta)
        [i] = double_min(REAL(getListElt(tmp, "deltas")), S_len - i);
        REAL(maxVega)
        [i] = double_max(REAL(getListElt(tmp, "vegas")), T_len - i);
        REAL(minVega)
        [i] = double_min(REAL(getListElt(tmp, "vegas")), T_len - i);

        UNPROTECT(5);
    }

    SET_ELT(ans, 0, undPL);
    SET_ELT(ans, 1, strdPL);
    SET_ELT(ans, 2, undTC);
    SET_ELT(ans, 3, strdTC);

    /* For debugging */
    SET_ELT(ans, 4, maxUndPosn);
    SET_ELT(ans, 5, minUndPosn);
    SET_ELT(ans, 6, maxStrdPosn);
    SET_ELT(ans, 7, minStrdPosn);

    SET_ELT(ans, 8, maxDelta);
    SET_ELT(ans, 9, minDelta);
    SET_ELT(ans, 10, maxVega);
    SET_ELT(ans, 11, minVega);

    UNPROT2;

    return ans;
}
#undef ANS_LEN

static void optimGa3_chk(SEXP *initPar, const int parLen, const int basePopSize,
                         int *numprot)
{
    long    i = 0, j = 0;
    long    initLen = 0;
    double *dptr;
    SEXP    tmp;

    if (basePopSize < 1 || parLen < 1)
        error(
            "Must have positive population size and positive parameter length");

    PROT2(tmp = NEW_LIST(basePopSize), *numprot);

    GetRNGstate();

    for (j = 0; j < basePopSize; j++)
    {
        SET_ELT(tmp, j, NEW_NUMERIC(parLen));

        dptr = REAL(GET_ELT(tmp, j));

        for (i = 0; i < parLen; i++) *dptr++ = norm_rand();
    }

    if (isMatrix(*initPar))
    {
        if (parLen != nrows(*initPar))
            error("Inconsistent parameter length info");

        for (j = 0; j < MIN(ncols(*initPar), basePopSize); j++)
            memcpy(REAL(GET_ELT(tmp, j)), matcol1(*initPar, j),
                   parLen * sizeof(double));

    } else if (length(*initPar))
    {
        ENSURE_NUMERIC(*initPar, *numprot);

        if (length(*initPar) != parLen)
            error("Inconsistent parameter length info");

        SET_ELT(tmp, 0, *initPar);
    }

    *initPar = tmp;

    PutRNGstate();
}

static void mutate3(SEXP pop, const int basePopSize)

/*******************************************************************
 *
 *  Description: Mutates population of elts 0 through
 *    basePopSize - 1 in pop putting the results in elts
 *    basePopSize throught 2 * basePopSize - 1.
 *
 *******************************************************************/

{
    int              ok, i, j;
    const int        parLen = length(GET_ELT(pop, 0));
    register double *ptr1 = NULL, *ptr2 = NULL;
    SEXP             sxp1 = R_NilValue, sxp2 = R_NilValue;

    ok = isNewList(pop) && basePopSize > 0 && length(pop) >= 2 * basePopSize;

    if (!ok) error("bad argument passed to mutate2");

    GetRNGstate();

    for (j = 0; j < basePopSize; j++)
    {
        sxp1 = GET_ELT(pop, j);
        sxp2 = GET_ELT(pop, j + basePopSize);

        if (length(sxp1) != parLen || length(sxp2) != parLen)
            error("unequal population lengths in mutate3");

        if (!IS_NUMERIC(sxp1) || !IS_NUMERIC(sxp2))
            error("non-numeric population element in mutate3");

        ptr1 = REAL(sxp1);
        ptr2 = REAL(sxp2);

        for (i = 0; i < parLen; i++, ptr1++, ptr2++)
            *ptr2 =
                ISNA(*ptr1) ? norm_rand() : *ptr1 + norm_rand() * ABS(*ptr1);
    }

    PutRNGstate();
}

static void breed3(SEXP pop, const int basePopSize)

/*******************************************************************
 *
 *  Description: Breeds population of column vectors 0 through
 *
 *******************************************************************/

{
    const int bigPopSize = length(pop);
    const int parLen     = length(GET_ELT(pop, 0));

    int counter;

    register int     i, j, k, ok;
    register double  w;
    register double *ptr1, *ptr2, *kidPtr;

    ok = isNewList(pop) && parLen && bigPopSize == BIG_POP_SIZE(basePopSize);

    if (!ok) error("bad arguments in breed3");

    GetRNGstate();

    counter = 0;

    for (i = 0; i < 2 * basePopSize - 1; i++)
        for (j = i + 1; j < 2 * basePopSize; j++)
        {
            ptr1   = REAL(GET_ELT(pop, i));
            ptr2   = REAL(GET_ELT(pop, j));
            kidPtr = REAL(GET_ELT(pop, 2 * basePopSize + counter));

            for (k = 0; k < parLen; k++)
            {
                w = norm_rand() + 0.5;

                *kidPtr++ = w * (*ptr1++) + (1.0 - w) * (*ptr2++);
            }

            if (++counter > bigPopSize - 2 * basePopSize)
                error("breed3 ... something wrong here");
        }

    PutRNGstate();
}

static void getNextGen3(optimFun f, SEXP controlPar, SEXP pop, SEXP fvals,
                        SEXP ix, SEXP tmp, const int basePopSize,
                        const int minimize)
{
    const int    bigPopSize = length(pop);
    const int    parLen     = length(GET_ELT(pop, 0));
    const double sgn        = minimize ? 1.0 : -1.0;

    int    ok, *ixPtr;
    long   i, j;
    double val;

    ok = isNewList(pop) && isNewList(tmp) && IS_INTEGER(ix) &&
         IS_NUMERIC(fvals) && bigPopSize == length(ix) &&
         bigPopSize == length(fvals) &&
         bigPopSize == BIG_POP_SIZE(basePopSize) &&
         basePopSize == length(tmp) && parLen == length(GET_ELT(tmp, 0));

    if (!ok) error("bad arguments to getNextGen");

    for (i = 0; i < bigPopSize; i++) INTEGER(ix)[i] = i;

    mutate3(pop, basePopSize);
    breed3(pop, basePopSize);

    for (i = 0; i < bigPopSize; i++)
    {
        val = f(GET_ELT(pop, i), controlPar);

        REAL(fvals)[i] = R_FINITE(val) ? sgn * val : HUGE;
    }

    /* Move top basePopSize population members to front */
    rsort_with_index(REAL(fvals), INTEGER(ix), bigPopSize);

    ixPtr = INTEGER(ix);

    for (j = 0; j < basePopSize; j++, ixPtr++)
    {
        memcpy(REAL(GET_ELT(tmp, j)), REAL(GET_ELT(pop, *ixPtr)),
               parLen * sizeof(double));
    }

    for (j = 0; j < basePopSize; j++)
    {
        memcpy(REAL(GET_ELT(pop, j)), REAL(GET_ELT(tmp, j)),
               parLen * sizeof(double));
    }

    /* Done sorting population */

    /* Scale back 1st basePopSize function values */
    for (i = 0; i < basePopSize; i++) REAL(fvals)[i] *= sgn;
}

#define ANS_LEN 3
SEXP optimGa3(optimFun f, SEXP initPar, SEXP controlPar, const int parLen,
              const int basePopSize, const long maxit, const int minimize,
              const double tol, const int relTol)

/*******************************************************************
 *
 *  Description: Uses genetic algorithms with gradient descent
 *    at the end to optimize function x |-> f(x, controlPar).
 *
 *******************************************************************/

{
    int numprot = 0;

    const int bigPopSize = BIG_POP_SIZE(basePopSize);

    int j, gen = 1, notConvergedYet = 1, numAnsCols = 1;

    register double fmin, fmax;
    char           *names[] = {"par", "value", "convergence"};
    int            *dimPtr  = NULL;

    SEXP ans, pop, popMatrix, fvals, tmp, ix, popDim;

    /* Check arguments and perform necessary adjustments */
    optimGa3_chk(&initPar, parLen, basePopSize, &numprot);

    PROT2(fvals = NEW_NUMERIC(bigPopSize), numprot);
    PROT2(ix = NEW_INTEGER(bigPopSize), numprot);
    PROT2(ans = NEW_LIST(ANS_LEN), numprot);

    PROT2(pop = NEW_LIST(bigPopSize), numprot);
    PROT2(tmp = NEW_LIST(basePopSize), numprot);

    for (j = 0; j < basePopSize; j++)
    {
        SET_ELT(tmp, j, NEW_NUMERIC(parLen));
        SET_ELT(pop, j, NEW_NUMERIC(parLen));

        memcpy(REAL(GET_ELT(pop, j)), REAL(GET_ELT(initPar, j)),
               parLen * sizeof(double));
    }

    for (j = basePopSize; j < bigPopSize; j++)
    {
        SET_ELT(pop, j, NEW_NUMERIC(parLen));
    }

    GetRNGstate();

    for (gen = 0; gen < maxit && notConvergedYet; gen++)
    {
        getNextGen3(f, controlPar, pop, fvals, ix, tmp, basePopSize, minimize);

        fmin = REAL(fvals)[0];
        fmax = REAL(fvals)[basePopSize - 1];

        notConvergedYet = (ABS(fmax - fmin) >= tol);
    }

    PutRNGstate();

    numAnsCols = notConvergedYet ? basePopSize : 1;

    PROT2(popMatrix = allocMatrix(REALSXP, parLen, numAnsCols), numprot);

    for (j = 0; j < numAnsCols; j++)
        memcpy(matcol1(popMatrix, j), REAL(GET_ELT(pop, j)),
               parLen * sizeof(double));

    SET_ELT(ans, 0, popMatrix);
    SET_ELT(ans, 1, SET_LENGTH(fvals, numAnsCols));
    SET_ELT(ans, 2, ScalarInteger(notConvergedYet));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}
#undef ANS_LEN

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains some experimental weighted MLE type
**  of stuff.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#define WEIGHT_INDEX 6

static void ng11_w1(const long T, const double *x, const double *par, double w0,
                    int llonly, double *h, double *z, double *ll)

/********************************************************************
 *
 *   Description: Evaluation of loglikelihood, conditional variance
 *    and residuals for the given time series x[0], ..., x[T - 1]
 *    and the given coefficients par[0], ..., par[NUM_NG11_PAR-1].
 *
 ********************************************************************/

{
    int inc = 0;

    const double   *Lh   = NULL;
    double const   *xT   = x + T;
    double          a[1] = {0.0};
    double          b[1] = {0.0};
    register double lam, a0, a1, b1, gam, h1, w, wsum = 0;
    int             ok;

    inc = !(llonly = llonly || !h || !z);

    llonly ? h = (double *)a, z = (double *)b : 0;

    ok = x && par && ll;

    if (!ok) error("null pointer passed where not allowed");

    lam = par[LAMBDA_INDEX];
    a0  = par[A0_INDEX];
    a1  = par[A1_INDEX];
    b1  = par[B1_INDEX];
    gam = par[GAMMA_INDEX];
    h1  = par[H1_INDEX];

    ok = R_FINITE(T) && R_FINITE(lam) && R_FINITE(a0) && R_FINITE(b1) &&
         R_FINITE(a1) && R_FINITE(gam) && 0 <= w0 && w0 <= 1;

    if (!ok)
    {
        *ll = HUGE_NEGATIVE;
        warning("improper arguments, assuming loglikelihood = %.0f",
                HUGE_NEGATIVE);

        return;
    }

    if (!R_FINITE(h1))
    {
        warning("Assuming h[1] = E(h)");

        h1 = a0 / (1.0 - b1 - a1 * (1.0 + SQR(gam)));
    }

    *h   = h1;
    *z   = *x / sqrt(*h) - lam + 0.5 * sqrt(*h);
    wsum = w = R_pow_di(w0, T - 1);
    *ll      = w * (log(M_2PI) + log(*h) + SQR(*z));

    Lh  = h;
    x  += 1;
    h  += inc;

    for (; x < xT && R_FINITE(*ll); x++, h += inc, Lh += inc)
    {
        *h = a0 + Lh[0] * (b1 + a1 * SQR(*z - gam));

        z += inc;

        *z = *x / sqrt(*h) - lam + 0.5 * sqrt(*h);

        w    /= w0;
        wsum += w;

        *ll += w * (log(M_2PI) + log(*h) + SQR(*z));
    }

    *ll = R_FINITE(*ll) ? -0.5 * (T / wsum) * (*ll) : HUGE_NEGATIVE;
}

static double ll_ngarch11_w1(SEXP par, SEXP x)
{
    int    numprot = 0;
    int    ok;
    double ll = 0.0;

    ok = IS_NUMERIC(x) && IS_NUMERIC(par);

    if (!ok) error("bad args to wll_ngarch11");

    ENSURE_LENGTH(par, NUM_NG11_PAR, numprot);

    ng11_w1(length(x), REAL(x), REAL(par), REAL(par)[WEIGHT_INDEX],
            /* llonly = */ 1,
            /* h = */ (double *)NULL,
            /* z = */ (double *)NULL, &ll);

    UNPROTECT(numprot);

    return ll;
}

static void adjustInitPar_w1(SEXP *initPar, int parLen, int *numprot)
{
    SEXP par = R_NilValue;

    if (isNewList(*initPar))
    {
        PROT2(par = NEW_NUMERIC(parLen), *numprot);

        REAL(par)[LAMBDA_INDEX] = asReal(getListElt(*initPar, "lambda"));
        REAL(par)[A0_INDEX]     = asReal(getListElt(*initPar, "a0"));
        REAL(par)[A1_INDEX]     = asReal(getListElt(*initPar, "a1"));
        REAL(par)[B1_INDEX]     = asReal(getListElt(*initPar, "b1"));
        REAL(par)[GAMMA_INDEX]  = asReal(getListElt(*initPar, "gamma"));
        REAL(par)[WEIGHT_INDEX] = asReal(getListElt(*initPar, "weight"));

        if (H1_INDEX < parLen)
            REAL(par)[H1_INDEX] = asReal(getListElt(*initPar, "h1"));

        SET_NA_TO_ZERO(REAL(par)[LAMBDA_INDEX]);
        SET_NA_TO_ZERO(REAL(par)[A1_INDEX]);
        SET_NA_TO_ZERO(REAL(par)[B1_INDEX]);
        SET_NA_TO_ZERO(REAL(par)[GAMMA_INDEX]);

        if (!R_FINITE(REAL(par)[WEIGHT_INDEX])) REAL(par)[WEIGHT_INDEX] = 1;

        *initPar = par;
    }
}

#define ANS_LEN 7
SEXP fit_ngarch11_w1(SEXP x, SEXP initPar, SEXP fitInit, SEXP basePopSize,
                     SEXP tol, SEXP stopLags, SEXP minit, SEXP maxit,
                     SEXP options)
{
    const int    parLen       = NUM_NG11_PAR + !!asInteger(fitInit);
    const int    ibasePopSize = asInteger(basePopSize);
    const long   istopLags    = asInteger(stopLags);
    const long   imaxit       = asInteger(maxit);
    const long   iminit       = asInteger(minit);
    const double dtol         = asReal(tol);

    int   numprot = 0;
    int   gradToo;
    long  gradMaxIt;
    char *names[ANS_LEN] = {"par", "ll",  "convergence", "x",
                            "h",   "res", "gradconv"};
    SEXP  ans, fit1, fit2, tmp, dimNames, parNames;

    ENSURE_NUMERIC(x, numprot);

    adjustInitPar_w1(&initPar, parLen, &numprot);

    /* Do the fitting! */
    fit1 = optimGa2(ll_ngarch11_w1, initPar, x, parLen, ibasePopSize, istopLags,
                    iminit, imaxit, 0, dtol, 0);

    PROT2(fit1, numprot);
    PROT2(ans = NEW_LIST(ANS_LEN), numprot);
    PROT2(parNames = NEW_STRING(parLen), numprot);
    PROT2(dimNames = NEW_LIST(2), numprot);

    gradToo = asInteger(getListElt(options, "grad"));
    gradToo = SET_NA_TO_ZERO(gradToo);

    if (gradToo)
    {
        gradMaxIt = asInteger(getListElt(options, "maxit"));

        if (!R_FINITE((double)gradMaxIt) || gradMaxIt <= 0) gradMaxIt = imaxit;

        fit2 = optimGradient2(ll_ngarch11_w1, getListElt(fit1, "par"), x, dtol,
                              0, 0, gradMaxIt);
        PROT2(fit2, numprot);

    } else
    {
        fit2 = fit1;
    }

    PROT2(tmp = ngarch11(x, getListElt(fit2, "par"), ScalarInteger(0)),
          numprot);

    SET_ELT(ans, 0, getListElt(fit2, names[0]));
    SET_ELT(ans, 1, getListElt(fit2, "value"));
    SET_ELT(ans, 2, getListElt(fit1, names[2]));
    SET_ELT(ans, 3, x);
    SET_ELT(ans, 4, getListElt(tmp, names[4]));
    SET_ELT(ans, 5, getListElt(tmp, names[5]));
    SET_ELT(ans, 6, getListElt(fit2, "convergence"));

    set_names(ans, names);

    CHAR_PTR(parNames)[LAMBDA_INDEX] = mkChar("lambda");
    CHAR_PTR(parNames)[A0_INDEX]     = mkChar("a0");
    CHAR_PTR(parNames)[A1_INDEX]     = mkChar("a1");
    CHAR_PTR(parNames)[B1_INDEX]     = mkChar("b1");
    CHAR_PTR(parNames)[GAMMA_INDEX]  = mkChar("gamma");
    CHAR_PTR(parNames)[WEIGHT_INDEX] = mkChar("weight");

    if (NUM_NG11_PAR == parLen) CHAR_PTR(parNames)[H1_INDEX] = mkChar("h1");

    SET_ELT(dimNames, 0, parNames);
    SET_ELT(dimNames, 1, R_NilValue);

    setAttrib(GET_ELT(ans, 0), R_DimNamesSymbol, dimNames);

    UNPROTECT(numprot);

    return ans;
}
#undef ANS_LEN

#undef WEIGHT_INDEX

#endif /* YOU_SUCK */
