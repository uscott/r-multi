#ifdef FUCKOFF
#include "uri.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
** The following section deals with the NGARCH(1, 1) model
** augmented by a single regressor
**
**     x[t]  = m0 + lambda * sqrt(h[t]) - 0.5 * h[t]
**               + phi * x[t - 1]
**               + sqrt(h[t]) * z[t]
**
**     h[t]  = a00 + a01 * y[t] + h[t-1] * (b1 + a1 * (z[t-1] - gamma)^2)
**
**     where y is some other time series such that y[t] is in
**     the time t-1 info set.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

static int ng11_aug(const int *T, const double *x, const double *y,
                    const double *m0, const double *lambda, const double *a00,
                    const double *a01, const double *a1, const double *b1,
                    const double *gamma, const double *h1, double *h, double *z,
                    double *ll)
{
    const int llonly = !h || !z;
    const int inc    = !llonly;

    const double *Lx = NULL;
    const double *Lh = NULL;
    double const *xT = x + *T;

    double a[] = {0.0}; /* Dummy arrays basically */
    double b[] = {0.0};

    register double a0 = 0.0;

    int ok, uncond;

    llonly ? h = (double *)a, z = (double *)b : 0;

    ok = T && x && m0 && lambda && a00 && b1 && a1 && gamma && ll && a01 && y;

    if (!ok) error("null pointer passed where not allowed");

    ok = R_FINITE(*T) && R_FINITE(*m0) && R_FINITE(*lambda) && R_FINITE(*a00) &&
         R_FINITE(*b1) && R_FINITE(*a1) && R_FINITE(*gamma) && R_FINITE(*a01);

    if (!ok)
    {
        *ll = HUGE_NEGATIVE;

        warning("NA/NaN/infinite parameters");
        warning("assuming loglikelihood = HUGE_NEGATIVE");

        return 0;
    }

    uncond = !h1 || !R_FINITE(*h1);

    a0 = (*a00) + (*a01) * (*y++);

    if (uncond)
    {
        warning("NULL/NA/NaN/infinite h1 passed");
        warning("h(0) = E(h)");

        h[0] = a0 / (1.0 - *b1 - *a1 * (1.0 + SQR(*gamma)));

        if (h[0] <= 0) h[0] = sum_product(x, x, *T) / *T - dmean(x, *T);

    } else
    {
        h[0] = *h1;
    }

    *z = (*x - *m0) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);

    *ll = (*T) * log(M_2PI) + log(*h) + SQR(*z);

    Lx  = x;
    Lh  = h;
    x  += 1;
    h  += inc;

    for (0; x < xT && R_FINITE(*ll); x++, Lx++, h += inc, Lh += inc)
    {
        a0 = (*a00) + (*a01) * (*y++);

        *h = a0 + *Lh * (*b1 + *a1 * SQR(*z - *gamma));

        z += inc;

        *z = (*x - *m0) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);

        *ll += log(*h) + SQR(*z);
    }

    *ll = R_FINITE(*ll) ? -0.5 * (*ll) : HUGE_NEGATIVE;

    return 1;
}

#define N_SEXP_ARGS 3
#define NPAR        8
static void ngarch11_aug_chkargs(SEXP *x, SEXP *y, SEXP *par, int *numprot)
{
    int i;

    SEXP *args[N_SEXP_ARGS] = {x, y, par};

    for (i = 0; i < N_SEXP_ARGS; i++)
        if (!IS_NUMERIC(*args[i]))
            PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (length(*x) != length(*y)) error("length(x) != length(y)");

    if (NPAR == length(*par) + 1)
        PROT2(*par = SET_LENGTH(*par, NPAR), *numprot);

    if (NPAR != length(*par)) error("NPAR != length(*par)");
}
#undef N_SEXP_ARGS
#undef NPAR

SEXP ngarch11_aug(SEXP x, SEXP y, SEXP par, SEXP llonly)

/*******************************************************************
 *
 *  Description:
 *
 *    Returns list with elements named "h", "res", "ll".
 *
 *  Inputs:
 *    x
 *    - Univariate time series (i.e. numeric vector).
 *    y
 *    - Time series of regressors.
 *    par
 *    - Vector of length
 *
 *******************************************************************/

{
    double  ll;
    long    T;
    int     numprot = 0;
    double *zPtr = NULL, *hPtr = NULL;
    char   *names[THREE] = {"h", "res", "ll"};

    SEXP h = R_NilValue, z = R_NilValue, ans;

    ngarch11_aug_chkargs(&x, &y, &par, &numprot);

    T = length(x);

    if (!asInteger(llonly))
    {
        PROT2(h = NEW_NUMERIC(T), numprot);
        PROT2(z = NEW_NUMERIC(T), numprot);

        hPtr = REAL(h);
        zPtr = REAL(z);
    }

    PROT2(ans = NEW_LIST(THREE), numprot);

    ng11_aug(&T, REAL(x), REAL(y), REAL(par) + 0, REAL(par) + 1, REAL(par) + 2,
             REAL(par) + 3, REAL(par) + 4, REAL(par) + 5, REAL(par) + 6,
             REAL(par) + 7, hPtr, zPtr, &ll);

    SET_ELT(ans, 0, h);
    SET_ELT(ans, 1, z);
    SET_ELT(ans, 2, ScalarReal(ll));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

static void sim_aug_ng11_path(const long T, const double *y, const double *par,
                              const double *x0, const double *h0,
                              const double *z, double *x, double *h)
{
    const double *Lh = NULL;
    double const *hE = h + T;

    double          z0;
    register double a0;
    double          mu, lambda, a00, a01, a1, b1, gamma;

    int ok;

    ok = par && y && x && h && h0 && x0;

    if (!ok) error("sim_aug_ng11_path: Null pointer passed where not allowed");

    ok = R_FINITE(*x0) && R_FINITE(*h0);

    if (!ok) error("sim_aug_ng11_path: Improper x0 or h0");

    mu     = par[0];
    lambda = par[1];
    a00    = par[2];
    a01    = par[3];
    a1     = par[4];
    b1     = par[5];
    gamma  = par[6];

    GetRNGstate();

    z0  = *x0 - (mu - 0.5 * (*h0) + lambda * sqrt(*h0));
    z0 /= sqrt(*h0);

    a0 = a00 + a01 * y[0];
    *h = a0 + (*h0) * (b1 + a1 * SQR(z0 - gamma));
    *x = mu - 0.5 * (*h) + sqrt(*h) * (lambda + *z);

    h  += 1;
    x  += 1;
    y  += 1;
    Lh  = h - 1;

    for (0; h < hE && R_FINITE(*h); x++, h++, Lh++)
    {
        a0 = a00 + a01 * (*y++);

        *h = a0 + (*Lh) * (b1 + a1 * SQR(*z - gamma));

        z += 1;

        *x = mu - 0.5 * (*h) + sqrt(*h) * (lambda + *z);
    }

    PutRNGstate();
}

#define NPAR  7
#define NARGS 7
static void sim_aug_ng11_chkargs(SEXP *T, SEXP *npaths, SEXP *y, SEXP *par,
                                 SEXP *x0, SEXP *h0, SEXP *z, int *numprot)
{
    int ok, i;

    SEXP *args[NARGS] = {T, npaths, y, par, x0, h0, z};

    char *names[NARGS] = {"T", "npaths", "y", "par", "x0", "h0", "z"};

    for (i = 0; i < NARGS - 1; i++)
        if (isNull(*args[i]))
            error("Argument %s not allowed to be NULL", names[i]);

    ok = asInteger(*T) > 0 && asInteger(*npaths) > 0;

    if (!ok) error("Non-positive number of paths or path length.");

    if (NPAR != length(*par)) error("Wrong number of parameters passed");

    for (i = 2; i < NARGS; i++)
        if (!IS_NUMERIC(*args[i]))
            PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (1 == length(*y)) PROT2(*y = numrep(*y, asInteger(*T)), *numprot);

    if (asInteger(*T) != length(*y)) error("y has wrong length");
}

SEXP sim_aug_ng11(SEXP T, SEXP npaths, SEXP y, SEXP par, SEXP x0, SEXP h0,
                  SEXP z)
{
    int  zlen, prevolve;
    int  numprot = 0;
    long pathLen, numPaths;
    long j;

    double mu, lambda, a00, a01, a1, b1, gamma;

    double *uPtr = NULL, *xPtr = NULL, *hPtr = NULL;

    SEXP ans, x, h;

    char *names[2 + NARGS] = {"x",   "h",  "T",  "npaths", "y",
                              "par", "x0", "h0", "z"};

    sim_aug_ng11_chkargs(&T, &npaths, &y, &par, &x0, &h0, &z, &numprot);

    pathLen  = (long)asInteger(T);
    numPaths = (long)asInteger(npaths);

    zlen = length(z);

    PROT2(ans = NEW_LIST(2 + NARGS), numprot);

    PROT2(x = allocMatrix(REALSXP, pathLen, numPaths), numprot);
    PROT2(h = allocMatrix(REALSXP, pathLen, numPaths), numprot);

    GetRNGstate();

    uPtr = (double *)R_alloc(pathLen, sizeof(double));
    xPtr = (double *)R_alloc(pathLen, sizeof(double));
    hPtr = (double *)R_alloc(pathLen, sizeof(double));

    for (j = 0; j < numPaths; j++)
    {
        get_random_sample2(REAL(z), uPtr, zlen, pathLen);

        sim_aug_ng11_path(pathLen, REAL(y), REAL(par), REAL(x0), REAL(h0), uPtr,
                          xPtr, hPtr);

        memcpy(matcol1(x, j), xPtr, pathLen * sizeof(double));
        memcpy(matcol1(h, j), hPtr, pathLen * sizeof(double));
    }

    set_names(ans, names);

    PutRNGstate();

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, h);
    SET_ELT(ans, 2, T);
    SET_ELT(ans, 3, npaths);
    SET_ELT(ans, 4, y);
    SET_ELT(ans, 5, par);
    SET_ELT(ans, 6, x0);
    SET_ELT(ans, 7, h0);
    SET_ELT(ans, 8, z);

    UNPROTECT(numprot);

    return ans;
}
#undef NARGS
#undef NPAR

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
** The following section deals with the NGARCH(1, 1) model
** augmented by a single regressor with exponentially decaying
** influence
**
**     x[i]  = (m0 + lambda * sqrt(h[i]) - 0.5 * h[i]) * dt[i]
**               + sqrt(h[i] * dt[i]) * z[i]
**
**     h[i]  = a0[i] + h[i-1] * (b1 + a1 * (z[i-1] - gamma)^2)
**     a0[i] = a00 * (1 - beta^age[i]) + a01 * beta^age[i] * y[i]
**
**     where y is some other time series such that y[i] is in
**     the time i - age[i] info set.
**
**     dt[i] = amt of time in between x[i] and x[i-1].
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#define NPAR 9
static int ng11_decay(const int *xLen, const double *x, const double *dt,
                      const double *y, const double *age, const double *par,
                      double *h, double *z, double *ll)
{
    const int llonly = !h || !z;
    const int inc    = !llonly;

    const double *Lx = NULL;
    const double *Lh = NULL;
    double const *xT = x + *xLen;

    double m0, lambda, a00, a01, beta, a1, b1, gamma, h0;

    double a[] = {0.0}; /* General-purpose arrays */
    double b[] = {0.0};

    register double a0 = 0.0; /* Intercept term for h innovations */
    register double w1, w2;   /* Weights for intercept term       */

    int ok, i;

    llonly ? h = (double *)a, z = (double *)b : 0;

    ok = xLen && x && dt && par && ll && y;

    if (!ok) error("null pointer passed where not allowed");

    ok = R_FINITE(*xLen);

    for (i = 0; i < NPAR - 1; i++) ok *= R_FINITE(par[i]);

    if (!ok)
    {
        *ll = HUGE_NEGATIVE;

        warning("NA/NaN/infinite parameters, assuming ll = %f", HUGE_NEGATIVE);

        return 0;
    }

    /* Get individual parameters */
    m0     = par[0];
    lambda = par[1];
    a00    = par[2];
    a01    = par[3];
    beta   = par[4];
    a1     = par[5];
    b1     = par[6];
    gamma  = par[7];
    h0     = par[8];

    w2 = R_pow(beta, *age++);
    w1 = 1.0 - w2;

    a0 = w1 * a00 + w2 * a01 * (*y++);

    if (!R_FINITE(h0))
    {
        warning("NULL/NA/NaN/infinite h0 passed, assuming h(0) = E(h)");

        h[0] = a0 / (1.0 - b1 - a1 * (1.0 + SQR(gamma)));

        if (h[0] <= 0)
            h[0] = sum_product(x, x, *xLen) / xLen[0] - dmean(x, *xLen);

    } else
    {
        h[0] = h0;
    }

    *z  = *x - (m0 + lambda * sqrt(*h) - 0.5 * sqrt(*h)) * dt[0];
    *z /= sqrt(h[0] * dt[0]);

    *ll = xLen[0] * log(M_2PI) + log(*h) + SQR(*z);

    Lx  = x;
    Lh  = h;
    x  += 1;
    dt += 1;
    h  += inc;

    for (0; x < xT && R_FINITE(*ll); x++, Lx++, dt++, h += inc, Lh += inc)
    {
        w2 = R_pow(beta, *age++);
        w1 = 1.0 - w2;

        a0 = w1 * a00 + w2 * a01 * (*y++);

        *h = a0 + *Lh * (b1 + a1 * SQR(*z - gamma));

        z += inc;

        *z  = *x - (m0 + lambda * sqrt(*h) - 0.5 * h[0]) * dt[0];
        *z /= sqrt(h[0] * dt[0]);

        *ll += log(*h) + SQR(*z);
    }

    *ll = R_FINITE(*ll) ? -0.5 * (*ll) : HUGE_NEGATIVE;

    return 1;
}

#define N_SEXP_ARGS 5
static void ngarch11_decay_chkargs(SEXP *x, SEXP *dt, SEXP *y, SEXP *age,
                                   SEXP *par, int *numprot)
{
    long i;
    long xLen;

    SEXP *args[N_SEXP_ARGS] = {x, dt, y, age, par};

    xLen = length(*x);

    if (isNull(*dt))
    {
        warning("Null dt parameter, assuming all dt[i] == 1");

        PROT2(*dt = NEW_NUMERIC(xLen), *numprot);

        for (i = 0; i < xLen; i++) REAL(dt)[i] = 1.0;
    }

    for (i = 0; i < N_SEXP_ARGS; i++)
        if (!IS_NUMERIC(*args[i]))
            PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (xLen != length(*y)) error("length(x) != length(y)");

    if (NPAR == length(*par) + 1)
        PROT2(*par = SET_LENGTH(*par, NPAR), *numprot);

    if (NPAR != length(*par)) error("NPAR != length(*par)");
}
#undef N_SEXP_ARGS

SEXP ngarch11_decay(SEXP x, SEXP dt, SEXP y, SEXP age, SEXP par, SEXP llonly)

/*******************************************************************
 *
 *  Description:
 *
 *    Returns list with elements named "h", "res", "ll".
 *
 *  Inputs:
 *    x
 *    - Univariate time series (i.e. numeric vector).
 *    dt
 *    - Univariate time series giving amt of time since
 *    last measurement of x.
 *    y
 *    - Time series consisting of regressors.
 *    age
 *    - Time series consisting of age of corresponding value
 *    of y.
 *    par
 *    - Vector of length
 *
 *******************************************************************/

{
    double  ll;
    long    xLen;
    int     numprot = 0;
    double *zPtr = NULL, *hPtr = NULL;
    char   *names[THREE] = {"h", "res", "ll"};

    SEXP h = R_NilValue, z = R_NilValue, ans;

    ngarch11_decay_chkargs(&x, &dt, &y, &age, &par, &numprot);

    xLen = length(x);

    if (!asInteger(llonly))
    {
        PROT2(h = NEW_NUMERIC(xLen), numprot);
        PROT2(z = NEW_NUMERIC(xLen), numprot);

        hPtr = REAL(h);
        zPtr = REAL(z);
    }

    PROT2(ans = NEW_LIST(THREE), numprot);

    ng11_decay(&xLen, REAL(x), REAL(dt), REAL(y), REAL(age), REAL(par), hPtr,
               zPtr, &ll);

    SET_ELT(ans, 0, h);
    SET_ELT(ans, 1, z);
    SET_ELT(ans, 2, ScalarReal(ll));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

static void sim_decay_ng11_path(const long pathLen, const double *dt,
                                const double *y, const double *age,
                                const double *par, const double *z, double *x,
                                double *h)
{
    const double *Lh = NULL;
    double const *hE = h + pathLen;

    const double a[1] = {1.0}; /* General-purpose array */

    register double a0;
    register double w1, w2;

    double m0, lambda, a00, a01, beta, a1, b1, gamma, h0;

    int ok, inc;

    ok = par && y && age && b1 && gamma && x && h;

    if (!ok) error("sim_ng11: Null pointer passed where not allowed");

    if (!dt)
    {
        dt  = (double *)a;
        inc = 0;

    } else
        inc = 1;

    /* Get individual parameters */
    m0     = par[0];
    lambda = par[1];
    a00    = par[2];
    a01    = par[3];
    beta   = par[4];
    a1     = par[5];
    b1     = par[6];
    gamma  = par[7];
    h0     = par[8];

    GetRNGstate();

    /* Get initial states */
    *h = h0;

    *x  = (m0 + lambda * sqrt(*h) - 0.5 * h[0]) * dt[0]; /* mean */
    *x += sqrt(h[0] * dt[0]) * z[0];                     /* noise */

    h  += 1;
    x  += 1;
    y  += 1;
    dt += inc;
    Lh  = h - 1;

    for (0; h < hE && R_FINITE(*h); x++, h++, Lh++, dt += inc)
    {
        w2 = R_pow(beta, *age++);
        w1 = 1.0 - w2;

        a0 = w1 * a00 + w2 * a01 * (*y++);

        *h = a0 + (*Lh) * (b1 + a1 * SQR(*z - gamma));

        z += 1;

        *x  = (m0 + lambda * sqrt(*h) - 0.5 * h[0]) * dt[0]; /* mean */
        *x += sqrt(h[0] * dt[0]) * z[0];                     /* noise */
    }

    PutRNGstate();
}

#define NARGS 7
static void sim_decay_ng11_chkargs(SEXP *pathLen, SEXP *numPaths, SEXP *dt,
                                   SEXP *y, SEXP *age, SEXP *par, SEXP *z,
                                   int *numprot)
{
    int  ok, i;
    long len;

    SEXP *args[NARGS] = {pathLen, numPaths, dt, y, age, par, z};

    char *names[NARGS] = {"pathLen", "numPaths", "dt", "y", "age", "par", "z"};

    len = asInteger(*pathLen);

    for (i = 0; i < NARGS - 1; i++)
        if (isNull(*args[i]))
            error("Argument %s not allowed to be NULL", names[i]);

    ok = len > 0 && asInteger(*numPaths) > 0;

    if (!ok) error("Non-positive number of paths or path length.");

    for (i = 2; i < NARGS; i++)
        if (!IS_NUMERIC(*args[i]))
            PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    /* Make sure *y has correct length */
    if (1 == length(*y)) PROT2(*y = numrep(*y, len), *numprot);

    if (len != length(*y)) error("y has wrong length");

    /* Make sure *dt has correct length */
    if (!length(*dt))
    {
        PROT2(*dt = NEW_NUMERIC(len), *numprot);

        for (i = 0; i < len; i++) REAL(*dt)[i] = 1.0;

    } else if (1 == length(*dt))
        PROT2(*dt = numrep(*dt, len), *numprot);

    if (len != length(*dt)) error("dt has wrong length");

    /* Make sure *age has correct length */
    if (1 == length(*age))
    {
        PROT2(*age = SET_LENGTH(*age, len), *numprot);

        for (i = 1; i < len; i++)
            REAL(*age)[i] = REAL(*age)[i - 1] + REAL(*dt)[i];
    }

    if (len != length(*age)) error("age has wrong length");
}

SEXP sim_decay_ng11(SEXP T, SEXP numPaths, SEXP dt, SEXP y, SEXP age, SEXP par,
                    SEXP z)
{
    int  zlen, prevolve;
    int  numprot = 0;
    long pathLen, nPaths;
    long j;

    double *uPtr = NULL, *xPtr = NULL, *hPtr = NULL;

    SEXP ans, x, h;

    char *names[2 + NARGS] = {"x", "h",   "T",   "numPaths", "dt",
                              "y", "age", "par", "z"};

    sim_decay_ng11_chkargs(&T, &numPaths, &dt, &y, &age, &par, &z, &numprot);

    pathLen = (long)asInteger(T);
    nPaths  = (long)asInteger(numPaths);

    zlen = length(z);

    PROT2(ans = NEW_LIST(2 + NARGS), numprot);

    PROT2(x = allocMatrix(REALSXP, pathLen, nPaths), numprot);
    PROT2(h = allocMatrix(REALSXP, pathLen, nPaths), numprot);

    GetRNGstate();

    uPtr = (double *)malloc(pathLen * sizeof(double));
    xPtr = (double *)malloc(pathLen * sizeof(double));
    hPtr = (double *)malloc(pathLen * sizeof(double));

    !uPtr || !xPtr || !hPtr ? error("malloc failed") : 0;

    for (j = 0; j < nPaths; j++)
    {
        get_random_sample2(REAL(z), uPtr, zlen, pathLen);

        sim_decay_ng11_path(pathLen, REAL(dt), REAL(y), REAL(age), REAL(par),
                            uPtr, xPtr, hPtr);

        memcpy(matcol1(x, j), xPtr, pathLen * sizeof(double));
        memcpy(matcol1(h, j), hPtr, pathLen * sizeof(double));
    }

    set_names(ans, names);

    PutRNGstate();

    free(uPtr);
    free(xPtr);
    free(hPtr);

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, h);
    SET_ELT(ans, 2, T);
    SET_ELT(ans, 3, numPaths);
    SET_ELT(ans, 4, dt);
    SET_ELT(ans, 5, y);
    SET_ELT(ans, 6, age);
    SET_ELT(ans, 7, par);
    SET_ELT(ans, 8, z);

    UNPROTECT(numprot);

    return ans;
}
#undef NARGS
#undef NPAR

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
** The following section deals with the seasonal NGARCH(1, 1) model
** given by
**
**     x[t]  = m0 + lambda * sqrt(h[t]) - 0.5 * h[t]
**               + phi * x[t - 1]
**               + sqrt(h[t]) * z[t]
**
**     h[t]  = a0[t] + h[t-1] * (b1 + a1 * (z[t-1] - gamma)^2)
**     a0[t] = a00 * exp(a01 * cos(2*pi/365*(yday[t] - omega))
**
**      where yday[t] is the day of the year (0 - 365)
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

static int ng11_seasonal(const int *T, const double *x, const double *yday,
                         const double *m0, const double *lambda,
                         const double *a00, const double *a01,
                         const double *omega, const double *a1,
                         const double *b1, const double *gamma,
                         const double *h0, double *h, double *z, double *ll)
{
    const int llonly = !h || !z;
    const int inc    = !llonly;

    const double *Lx = NULL;
    const double *Lh = NULL;
    double const *xT = x + *T;

    double a[] = {0.0}; /* Dummy arrays basically */
    double b[] = {0.0};

    register double a0 = 0.0;

    int ok, uncond;

    llonly ? h = (double *)a, z = (double *)b : 0;

    ok = T && x && m0 && lambda && a00 && omega && b1 && a1 && gamma && ll &&
         a01 && yday;

    if (!ok) error("null pointer passed where not allowed");

    ok = R_FINITE(*T) && R_FINITE(*m0) && R_FINITE(*lambda) && R_FINITE(*a00) &&
         R_FINITE(*b1) && R_FINITE(*omega) && R_FINITE(*a1) &&
         R_FINITE(*gamma) && R_FINITE(*a01);

    if (!ok)
    {
        *ll = HUGE_NEGATIVE;

        warning("NA/NaN/infinite parameters");
        warning("assuming loglikelihood = HUGE_NEGATIVE");

        return 0;
    }

    uncond = !h0 || !R_FINITE(*h0);

    a0 = (*a00) * exp(*a01 * cos(M_2PI / 365.0 * (*yday++ - *omega)));

    if (uncond)
    {
        warning("NULL/NA/NaN/infinite xmin1 or h0 passed");
        warning("assuming x(-1) = E(x) and h(0) = E(h)");

        h[0] = a0 / (1.0 - *b1 - *a1 * (1.0 + SQR(*gamma)));

        if (h[0] <= 0) h[0] = sum_product(x, x, *T) / *T;

    } else
    {
        h[0] = *h0;
    }

    *z = (*x - *m0) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);

    *ll = (*T) * log(M_2PI) + log(*h) + SQR(*z);

    Lx  = x;
    Lh  = h;
    x  += 1;
    h  += inc;

    for (0; x < xT && R_FINITE(*ll); x++, Lx++, h += inc, Lh += inc)
    {
        a0 = (*a00) * exp(*a01 * cos(M_2PI / 365.0 * (*yday++ - *omega)));

        *h = a0 + *Lh * (*b1 + *a1 * SQR(*z - *gamma));

        z += inc;

        *z = (*x - *m0) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);

        *ll += log(*h) + SQR(*z);
    }

    *ll = R_FINITE(*ll) ? -0.5 * (*ll) : HUGE_NEGATIVE;

    return 1;
}

#define N_SEXP_ARGS 3
#define NPAR        9
static void ngarch11_seasonal_chkargs(SEXP *x, SEXP *yday, SEXP *par,
                                      int *numprot)
{
    int i;

    SEXP *args[N_SEXP_ARGS] = {x, yday, par};

    for (i = 0; i < N_SEXP_ARGS; i++)
        if (!IS_NUMERIC(*args[i]))
            PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    if (length(*x) != length(*yday)) error("length(x) != length(yday)");

    if (NPAR == length(*par) + 1)
        PROT2(*par = SET_LENGTH(*par, NPAR), *numprot);

    if (NPAR != length(*par)) error("NPAR != length(*par)");
}
#undef N_SEXP_ARGS
#undef NPAR

SEXP ngarch11_seasonal(SEXP x, SEXP yday, SEXP par, SEXP llonly)

/*******************************************************************
 *
 *  Description:
 *
 *    Returns list with elements named "h", "res", "ll".
 *
 *  Inputs:
 *    x
 *    - Univariate time series (i.e. numeric vector).
 *    yday
 *    - Time series consisting of days of year.
 *    par
 *    - Vector of length
 *
 *******************************************************************/

{
    double  ll;
    long    T;
    int     numprot = 0;
    double *zPtr = NULL, *hPtr = NULL;
    char   *names[THREE] = {"h", "res", "ll"};

    SEXP h = R_NilValue, z = R_NilValue, ans;

    ngarch11_seasonal_chkargs(&x, &yday, &par, &numprot);

    T = length(x);

    if (!asInteger(llonly))
    {
        PROT2(h = NEW_NUMERIC(T), numprot);
        PROT2(z = NEW_NUMERIC(T), numprot);

        hPtr = REAL(h);
        zPtr = REAL(z);
    }

    PROT2(ans = NEW_LIST(THREE), numprot);

    ng11_seasonal(&T, REAL(x), REAL(yday), REAL(par) + 0, REAL(par) + 1,
                  REAL(par) + 2, REAL(par) + 3, REAL(par) + 4, REAL(par) + 5,
                  REAL(par) + 6, REAL(par) + 7, REAL(par) + 8, hPtr, zPtr, &ll);

    SET_ELT(ans, 0, h);
    SET_ELT(ans, 1, z);
    SET_ELT(ans, 2, ScalarReal(ll));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
** The following section deals with the AR(1)-NGARCH(1, 1) model
** given by
**
**     x[t] = m0 + lambda * sqrt(h[t]) - 0.5 * h[t]
**               + phi * x[t - 1]
**               + sqrt(h[t]) * z[t]
**
**     h[t] = a0 + h[t-1] * (b1 + a1 * (z[t-1] - gamma)^2)
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void ar1ng11(const int *T, const double *x, const double *m0,
             const double *lambda, const double *phi, const double *a0,
             const double *a1, const double *b1, const double *gamma,
             const double *h0, const double *xmin1, double *h, double *z,
             double *ll)

/********************************************************************
 *
 *   Description: Evaluation of loglikelihood, conditional variance
 *    and residuals for the given time series x[0], ..., x[*T - 1]
 *    and the given coefficients *m0, *lambda, ...
 *
 ********************************************************************/

{
    const int llonly = !h || !z;
    const int inc    = !llonly;

    const double *Lx  = NULL;
    const double *Lh  = NULL;
    double const *xT  = x + *T;
    double        a[] = {0.0};
    double        b[] = {0.0};
    double        xm1, hin;
    int           ok, uncond;

    llonly ? h = (double *)a, z = (double *)b : 0;

    ok = T && x && m0 && lambda && phi && a0 && b1 && a1 && gamma && ll;

    if (!ok) error("null pointer passed where not allowed");

    ok = R_FINITE(*T) && R_FINITE(*m0) && R_FINITE(*lambda) && R_FINITE(*phi) &&
         R_FINITE(*a0) && R_FINITE(*b1) && R_FINITE(*a1) && R_FINITE(*gamma);

    if (ok)
    {
        uncond = !xmin1 || !h0 || !R_FINITE(*xmin1) || !R_FINITE(*h0);

        if (uncond)
        {
            warning("NULL/NA/NaN/infinite xmin1 or h0 passed");
            warning("assuming x(-1) = E(x) and h(0) = E(h)");
            xm1 = sum_array(x, *T) / (*T);
            hin = *a0 / (1.0 - *b1 - *a1 * (1.0 + SQR(*gamma)));

        } else
        {
            hin = *h0;
            xm1 = *xmin1;
        }

        *h  = hin;
        *z  = (*x - *m0 - *phi * xm1) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);
        *ll = *T * log(M_2PI) + log(*h) + SQR(*z);

        Lx  = x;
        Lh  = h;
        x  += 1;
        h  += inc;

        for (0; x < xT && R_FINITE(*ll); x++, Lx++, h += inc, Lh += inc)
        {
            *h = *a0 + *Lh * (*b1 + *a1 * SQR(*z - *gamma));

            z += inc;

            *z =
                (*x - *m0 - *phi * (*Lx)) / sqrt(*h) - *lambda + 0.5 * sqrt(*h);
            *ll += log(*h) + SQR(*z);
        }

        *ll = R_FINITE(*ll) ? -0.5 * (*ll) : HUGE_NEGATIVE;

    } else
    {
        *ll = HUGE_NEGATIVE;
        warning("NA/NaN/infinite parameters");
        warning("assuming loglikelihood = HUGE_NEGATIVE");
    }
}

#define LEN 14

static SEXP ar1ng11_matrix_apply(const SEXP x, const double *m0,
                                 const double *lambda, const double *phi,
                                 const double *a0, const double *a1,
                                 const double *b1, const double *gamma,
                                 const double *h0, const double *xmin1)
{
    int    T, n, j, numprot = 0, ok;
    double llsum      = 0.0;
    char  *names[LEN] = {"x",     "m0", "lambda", "phi", "a0",  "a1", "b1",
                         "gamma", "h0", "xmin1",  "h",   "res", "ll", "ll.sum"};

    SEXP ans, ll, h, res;

    if (!isMatrix(x)) error("C error: ng11_matrix_apply.  Non-matrix passed.");

    ok = m0 && lambda && phi && a0 && a1 && b1 && gamma && h0 && xmin1;

    if (!ok) error("Null pointer passed to ar1ng11_matrix_apply");

    T = nrows(x);
    n = ncols(x);

    PROT2(ans = NEW_LIST(LEN), numprot);
    PROT2(ll = NEW_NUMERIC(n), numprot);
    PROT2(h = allocMatrix(REALSXP, T, n), numprot);
    PROT2(res = allocMatrix(REALSXP, T, n), numprot);

    SET_DIMNAMES(h, GET_DIMNAMES(x));
    SET_DIMNAMES(res, GET_DIMNAMES(x));

    for (j = 0; j < n; j++)
    {
        ar1ng11(&T, matcol1(x, j), m0++, lambda++, phi++, a0++, b1++, a1++,
                gamma++, h0++, xmin1++, matcol1(h, j), matcol1(res, j),
                REAL(ll) + j);

        llsum += REAL(ll)[j];
    }

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, make_sexp_vec(m0 - n, n));
    SET_ELT(ans, 2, make_sexp_vec(lambda - n, n));
    SET_ELT(ans, 3, make_sexp_vec(phi - n, n));
    SET_ELT(ans, 4, make_sexp_vec(a0 - n, n));
    SET_ELT(ans, 5, make_sexp_vec(a1 - n, n));
    SET_ELT(ans, 6, make_sexp_vec(b1 - n, n));
    SET_ELT(ans, 7, make_sexp_vec(gamma - n, n));
    SET_ELT(ans, 8, make_sexp_vec(h0 - n, n));
    SET_ELT(ans, 9, make_sexp_vec(xmin1 - n, n));
    SET_ELT(ans, 10, h);
    SET_ELT(ans, 11, res);
    SET_ELT(ans, 12, ll);
    SET_ELT(ans, 13, ScalarReal(llsum));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

static SEXP ar1ng11_list_apply(const SEXP x, const double *m0,
                               const double *lambda, const double *phi,
                               const double *a0, const double *a1,
                               const double *b1, const double *gamma,
                               const double *h0, const double *xmin1)
{
    const int n = length(x);
    int       T, j, numprot = 0, ok;
    double    llsum      = 0.0;
    char     *names[LEN] = {"x",     "m0", "lambda", "phi", "a0",  "a1", "b1",
                            "gamma", "h0", "xmin1",  "h",   "res", "ll", "ll.sum"};

    SEXP ans, ll, h, res;

    if (!isNewList(x)) error("C error: ng11_lst.  Non-list passed.");

    ok = m0 && lambda && phi && a0 && a1 && b1 && gamma && h0 && xmin1;

    if (!ok) error("Null pointer passed to ar1ng11_list_apply");

    PROT2(ans = NEW_LIST(LEN), numprot);
    PROT2(ll = NEW_NUMERIC(n), numprot);
    PROT2(h = NEW_LIST(n), numprot);
    PROT2(res = NEW_LIST(n), numprot);

    SET_NAMES(h, GET_NAMES(x));
    SET_NAMES(res, GET_NAMES(x));

    for (j = 0; j < n; j++)
    {
        T = length(GET_ELT(x, j));

        SET_ELT(h, j, NEW_NUMERIC(T));
        SET_ELT(res, j, NEW_NUMERIC(T));

        SET_NAMES(GET_ELT(h, j), GET_NAMES(GET_ELT(x, j)));
        SET_NAMES(GET_ELT(res, j), GET_NAMES(GET_ELT(x, j)));

        ar1ng11(&T, REAL(GET_ELT(x, j)), m0++, lambda++, phi++, a0++, b1++,
                a1++, gamma++, h0++, xmin1++, REAL(GET_ELT(h, j)),
                REAL(GET_ELT(res, j)), REAL(ll) + j);

        llsum += REAL(ll)[j];
    }

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, make_sexp_vec(m0 - n, n));
    SET_ELT(ans, 2, make_sexp_vec(lambda - n, n));
    SET_ELT(ans, 3, make_sexp_vec(phi - n, n));
    SET_ELT(ans, 4, make_sexp_vec(a0 - n, n));
    SET_ELT(ans, 5, make_sexp_vec(a1 - n, n));
    SET_ELT(ans, 6, make_sexp_vec(b1 - n, n));
    SET_ELT(ans, 7, make_sexp_vec(gamma - n, n));
    SET_ELT(ans, 8, make_sexp_vec(h0 - n, n));
    SET_ELT(ans, 9, make_sexp_vec(xmin1 - n, n));
    SET_ELT(ans, 10, h);
    SET_ELT(ans, 11, res);
    SET_ELT(ans, 12, ll);
    SET_ELT(ans, 13, ScalarReal(llsum));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

#undef LEN

static void ar1ng11_apply_chkargs(SEXP *x, SEXP *par, SEXP *fn, SEXP *arity,
                                  SEXP *env, int *numprot)
{
    int ok;

    PROT2(*par = AS_NUMERIC(*par), *numprot);
    PROT2(*arity = AS_INTEGER(*arity), *numprot);

    ok = (isMatrix(*x) || isNewList(*x)) &&
         (length(par) == asInteger(*arity)) && isFunction(fn) &&
         isEnvironment(env);

    if (!ok) error("ar1ng11_apply.  Wrong argument types.");
}

SEXP ar1ng11_apply(const SEXP x, const SEXP par, const SEXP fn,
                   const SEXP arity, const SEXP env)

/********************************************************************
 *
 * Description: Computes the loglikelihoods, conditional variances
 *    and residuals of the
 *    NG(1, 1) model in ng11 for the columns of x
 *    (resp. elements of the list x).
 *    SEXP parameter is a function taking vectors of length
 *    asInteger(arity) and returning vectors of length
 *    6 * (isMatrix(x) ? ncols(x) : length(x)).
 *
 ********************************************************************/

{
    const int n       = isMatrix(x) ? ncols(x) : length(x);
    int       numprot = 0;
    double   *f       = NULL;
    SEXP      fvals, ans;

    ar1ng11_apply_chkargs(&x, &par, &fn, &arity, &env, &numprot);

    PROT2(fvals = feval(fn, par, env), numprot);

    if (length(fvals) != n * NINE)
        error("C error: ng11a.  Wrong return length for parameter fn.");

    f = REAL(fvals);

    if (isMatrix(x))
    {
        PROTECT(ans = ar1ng11_matrix_apply(x, f + 0 * n, f + 1 * n, f + 2 * n,
                                           f + 3 * n, f + 4 * n, f + 5 * n,
                                           f + 6 * n, f + 7 * n, f + 8 * n));
    } else
    {
        PROTECT(ans = ar1ng11_list_apply(x, f + 0 * n, f + 1 * n, f + 2 * n,
                                         f + 3 * n, f + 4 * n, f + 5 * n,
                                         f + 6 * n, f + 7 * n, f + 8 * n));
    }

    numprot += 1;

    UNPROTECT(numprot);

    return ans;
}

static void sim_ar1ng11_path(const long T, const double *m0,
                             const double *lambda, const double *phi,
                             const double *a0, const double *a1,
                             const double *b1, const double *gamma,
                             const double *h0, const double *xmin1,
                             const double *z, double *x, double *h)
{
    const double *Lh = NULL;
    const double *Lx = NULL;
    double const *hE = h + T;

    int ok;

    ok =
        m0 && lambda && phi && a0 && a1 && b1 && gamma && x && h && h0 && xmin1;

    if (!ok) error("sim_ng11: Null pointer passed where not allowed");

    GetRNGstate();

    *h = *h0;
    *x = *m0 + *phi * (*xmin1) - 0.5 * (*h) + sqrt(*h) * (*lambda + *z);

    h  += 1;
    x  += 1;
    Lx  = x - 1;
    Lh  = h - 1;

    for (0; h < hE && R_FINITE(*h); x++, h++, Lh++, Lx++)
    {
        *h = *a0 + (*Lh) * (*b1 + *a1 * SQR(*z - *gamma));

        z += 1;

        *x = *m0 + *phi * (*Lx) - 0.5 * (*h) + sqrt(*h) * (*lambda + *z);
    }

    PutRNGstate();
}

#define NARGS 12
static void sim_ar1ng11_chkargs(SEXP *T, SEXP *npaths, SEXP *m0, SEXP *lambda,
                                SEXP *phi, SEXP *a0, SEXP *a1, SEXP *b1,
                                SEXP *gamma, SEXP *h0, SEXP *xmin1, SEXP *z,
                                int *numprot)
{
    int ok, i;

    SEXP *args[] = {T,  npaths, m0,    lambda, phi,   a0,
                    a1, b1,     gamma, h0,     xmin1, z};

    char *names[] = {"T",  "npaths", "m0",    "lambda", "phi",   "a0",
                     "a1", "b1",     "gamma", "h0",     "xmin1", "z"};

    for (i = 0; i < NARGS - 4; i++)
        if (isNull(*args[i]))
            error("Argument %s not allowed to be NULL", names[i]);

    PROT2(*T = AS_INTEGER(*T), *numprot);
    PROT2(*npaths = AS_INTEGER(*npaths), *numprot);

    ok = asInteger(*T) > 0 && asInteger(*npaths) > 0;

    if (!ok) error("Non-positive number of paths or path length.");

    for (i = 2; i < NARGS; i++)
        PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    for (i = 2; i < NARGS - 1; i++)
        if (1 == length(*args[i]))
            PROT2(*args[i] = numrep(*args[i], asInteger(*npaths)), *numprot);

    for (i = 2; i < NARGS - 3; i++)
        if (asInteger(*npaths) != length(*args[i]))
            error("Argument %s has wrong length", names[i]);
}

SEXP sim_ar1ng11(SEXP T, SEXP npaths, SEXP m0, SEXP lambda, SEXP phi, SEXP a0,
                 SEXP a1, SEXP b1, SEXP gamma, SEXP h0, SEXP xmin1, SEXP z)
{
    int  zlen, prevolve;
    int  numprot = 0;
    long pathLen, numPaths, ulen;
    long j;

    double *u = NULL, *xx = NULL, *hh = NULL;

    double d; /* Random-task variable */

    SEXP ans, x, h;

    char *names[] = {"x",  "h",  "T",  "npaths", "m0", "lambda", "phi",
                     "a0", "a1", "b1", "gamma",  "h0", "xmin1",  "z"};

    sim_ar1ng11_chkargs(&T, &npaths, &m0, &lambda, &phi, &a0, &a1, &b1, &gamma,
                        &h0, &xmin1, &z, &numprot);

    pathLen  = (long)asInteger(T);
    numPaths = (long)asInteger(npaths);
    ulen     = pathLen;
    zlen     = length(z);

    PROT2(ans = NEW_LIST(2 + NARGS), numprot);

    PROT2(x = allocMatrix(REALSXP, pathLen, numPaths), numprot);
    PROT2(h = allocMatrix(REALSXP, pathLen, numPaths), numprot);

    GetRNGstate();

    /* See if prevolving is necessary */
    prevolve = numPaths != length(h0) || numPaths != length(xmin1) ||
               !IS_NUMERIC(h0) || !IS_NUMERIC(xmin1);

    for (j = 0; j < length(h0) && !prevolve; j++)
        prevolve += !R_FINITE(REAL(h0)[j]);

    for (j = 0; j < length(xmin1) && !prevolve; j++)
        prevolve += !R_FINITE(REAL(xmin1)[j]);

    if (prevolve)
    {
        warning("sim_ar1ng11:  doing unconditional simulation");

        ulen += NUM_PREVOLVE_STEPS;

        PROT2(h0 = NEW_NUMERIC(numPaths), numprot);
        PROT2(xmin1 = NEW_NUMERIC(numPaths), numprot);

        for (j = 0; j < numPaths; j++)
        {
            REAL(xmin1)[j] = 0.0;

            d = 1.0 - REAL(b1)[j] - REAL(a1)[j] * (1.0 + SQR(REAL(gamma)[j]));

            REAL(h0)[j] = REAL(a0)[j] / d;
        }
    }

    u  = (double *)malloc(ulen * sizeof(double));
    xx = (double *)malloc(ulen * sizeof(double));
    hh = (double *)malloc(ulen * sizeof(double));

    !u || !xx || !hh ? error("malloc failed") : 0;

    for (j = 0; j < numPaths; j++)
    {
        get_random_sample2(REAL(z), u, zlen, ulen);

        sim_ar1ng11_path(ulen, REAL(m0) + j, REAL(lambda) + j, REAL(phi) + j,
                         REAL(a0) + j, REAL(a1) + j, REAL(b1) + j,
                         REAL(gamma) + j, REAL(h0) + j, REAL(xmin1) + j, u, xx,
                         hh);

        memcpy(matcol1(x, j), xx + ulen - pathLen, pathLen * sizeof(double));
        memcpy(matcol1(h, j), hh + ulen - pathLen, pathLen * sizeof(double));
    }

    set_names(ans, names);

    PutRNGstate();

    free(u);
    free(xx);
    free(hh);

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, h);
    SET_ELT(ans, 2, T);
    SET_ELT(ans, 3, npaths);
    SET_ELT(ans, 4, m0);
    SET_ELT(ans, 5, lambda);
    SET_ELT(ans, 6, phi);
    SET_ELT(ans, 7, a0);
    SET_ELT(ans, 8, a1);
    SET_ELT(ans, 9, b1);
    SET_ELT(ans, 10, gamma);
    SET_ELT(ans, 11, h0);
    SET_ELT(ans, 12, xmin1);
    SET_ELT(ans, 13, z);

    UNPROTECT(numprot);

    return ans;
}
#undef NARGS

#endif /* FUCK_OFF */
