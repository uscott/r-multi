
#include "uri.h"

SEXP fillNAfwd(SEXP x)
{
    int    numprot = 0;
    long   xLen = 0, i = 0;
    double tmp = NA_REAL, *yPtr = NULL;
    SEXP   y = R_NilValue;

    ENSURE_NUMERIC(x, numprot);

    xLen = length(x);

    PROT2(y = NEW_NUMERIC(xLen), numprot);

    memcpy(REAL(y), REAL(x), xLen * sizeof(double));

    yPtr = REAL(y);

    for (i = 0; i < xLen; i++, yPtr++)
    {
        if (ISNAN(*yPtr) && !ISNAN(tmp))
            *yPtr = tmp;
        else
            tmp = *yPtr;
    }

    UNPROTECT(numprot);

    return y;
}

SEXP fillNAbkwd(SEXP x)
{
    int    numprot = 0;
    long   xLen = 0, i = 0;
    double tmp = NA_REAL, *yPtr = NULL;
    SEXP   y = R_NilValue;

    ENSURE_NUMERIC(x, numprot);

    xLen = length(x);

    PROT2(y = NEW_NUMERIC(xLen), numprot);

    memcpy(REAL(y), REAL(x), xLen * sizeof(double));

    yPtr = REAL(y) + xLen - 1;

    for (i = 0; i < xLen; i++, yPtr--)
    {
        if (ISNAN(*yPtr) && !ISNAN(tmp))
            *yPtr = tmp;
        else
            tmp = *yPtr;
    }

    UNPROTECT(numprot);

    return y;
}

SEXP varRatio(SEXP x, SEXP maxOrd)
{
    int     numprot = 0;
    int     ord, ordMax;
    long    xLen, i;
    double  var1;
    double *returns = NULL;
    SEXP    ans;

    PROT2(x = AS_NUMERIC(x), numprot);

    xLen   = length(x);
    ordMax = MIN(asInteger(maxOrd), xLen - 1);

    PROT2(ans = NEW_NUMERIC(ordMax), numprot);

    returns = (double *)R_alloc(xLen, sizeof(double));

    /* Get order 1 returns and 1st order variance: */

    for (i = 1; i < xLen; i++) returns[i] = REAL(x)[i] - REAL(x)[i - 1];

    var1 = dvar(returns + 1, xLen - 1);

    REAL(ans)[0] = 1.0;

    /* Get higher-order variances: */
    for (ord = 2; ord <= ordMax; ord++)
    {
        for (i = ord; i < xLen; i++) returns[i] = REAL(x)[i] - REAL(x)[i - ord];

        REAL(ans)[ord - 1] = dvar(returns + ord, xLen - ord) / (ord * var1);
    }

    UNPROTECT(numprot);

    return ans;
}

SEXP roll_apply(SEXP x, const SEXP fn, const SEXP order, const SEXP env)
{
    const int n   = length(x);
    const int ord = asInteger(order) - 1;

    int     t;
    double *xp, *ap, *tp;
    SEXP    ans, tmp;

    if (!isFunction(fn)) error("Bad function input in roll_apply");
    if (!isEnvironment(env)) error("Bad environment input in roll_apply");

    PROTECT(x = AS_NUMERIC(x));
    PROTECT(ans = NEW_NUMERIC(n));
    PROTECT(tmp = NEW_NUMERIC(ord + 1));

    for (t = 0, ap = REAL(ans); t < ord; t++) *ap++ = NA_REAL;

    for (t = n, xp = REAL(x), tp = REAL(tmp); t > ord; t--)
    {
        memcpy(REAL(tmp), xp++, (ord + 1) * sizeof(double));
        *ap++ = *REAL(feval(tmp, fn, env));
    }

    UNPROTECT(3);

    SET_NAMES(ans, GET_NAMES(x));

    return ans;
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  The following two functions are for stratified sampling
**  of uniform and normal random variables.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

SEXP unifuri(SEXP n, SEXP jiggle, SEXP permute, SEXP replace)
{
    const long m = (long)asInteger(n);
    SEXP       ans;

    PROTECT(ans = NEW_NUMERIC(m));

    unif_rand_uri(REAL(ans), m, asInteger(jiggle), asInteger(permute),
                  asInteger(replace));

    UNPROTECT(1);

    return ans;
}

SEXP normuri(SEXP n, SEXP jiggle, SEXP permute, SEXP replace)
{
    const long m = (long)asInteger(n);
    SEXP       ans;

    PROTECT(ans = NEW_NUMERIC(m));

    norm_rand_uri(REAL(ans), m, asInteger(jiggle), asInteger(permute),
                  asInteger(replace));

    UNPROTECT(1);

    return ans;
}

SEXP bvt_normuri(SEXP n, SEXP rho, SEXP jiggle)
{
    const long      m   = (long)asInteger(n);
    const long      N   = m * m;
    const int       jig = asInteger(jiggle);
    const double    mi  = 1.0 / (double)m;
    const double    r   = asReal(rho);
    register double u1, u2;
    long            i, j;
    double         *x1, *x2;

    SEXP ans;

    if (m < 1 || 1.0 < ABS(r)) error("bad args in bvt_normuri");

    PROTECT(ans = allocMatrix(REALSXP, N, 2));

    x1 = REAL(ans);
    x2 = REAL(ans) + N;

    GetRNGstate();

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            u1 = mi * ((double)i + (jig ? unif_rand() : 0.5));
            u2 = mi * ((double)j + (jig ? unif_rand() : 0.5));

            x1[i * m + j] = qnorm(u1, 0.0, 1.0, 1, 0);
            x2[i * m + j] = r * x1[i * m + j] +
                            sqrt(1.0 - SQR(r)) * qnorm(u2, 0.0, 1.0, 1, 0);
        }
    }

    PutRNGstate();

    UNPROTECT(1);

    return ans;
}

SEXP simulate_garch1(SEXP coefs, SEXP days, SEXP r0, SEXP h0)
/* Last modified 22-Aug-2002. */
{
    int          conditional = !isNull(r0) & !isNull(h0);
    int          t;
    const int    T         = *INTEGER(days);
    double      *r         = NULL;
    const double a0        = REAL(coefs)[0];
    const double a1        = REAL(coefs)[1];
    const double b1        = REAL(coefs)[2];
    const int    init_days = 100;
    double       r_init, h_init;

    SEXP ans, ans_names, y, h;

    if (!isNumeric(coefs) || !isInteger(days) || 3 != length(coefs))
        error("bad inputs in simulate_garch1");

    PROTECT(ans = NEW_LIST(2));
    PROTECT(y = NEW_NUMERIC(T));
    PROTECT(h = NEW_NUMERIC(T));

    r = (double *)R_alloc(T, sizeof(double));

    GetRNGstate();

    if (conditional)
    {
        r[0]       = asReal(r0);
        REAL(h)[0] = asReal(h0);

    } else
    {
        r_init = 0;
        h_init = a0 / (1 - a1 - b1);

        for (t = 0; t < init_days; t++)
        {
            h_init = a0 + a1 * r_init * r_init + b1 * h_init;
            r_init = sqrt(h_init) * norm_rand();
        }

        *r       = r_init;
        *REAL(h) = h_init;
    }

    *REAL(y) = 0;

    for (t = 1; t < T; t++)
    {
        REAL(h)[t] = a0 + a1 * r[t - 1] * r[t - 1] + b1 * REAL(h)[t - 1];
        r[t]       = sqrt(REAL(h)[t]) * norm_rand();
        REAL(y)[t] = REAL(y)[t - 1] + r[t];
    }
    PutRNGstate();

    SET_VECTOR_ELT(ans, 0, y);
    SET_VECTOR_ELT(ans, 1, h);

    PROTECT(ans_names = NEW_STRING(2));
    CHARACTER_POINTER(ans_names)[0] = mkChar("y");
    CHARACTER_POINTER(ans_names)[1] = mkChar("h");
    SET_NAMES(ans, ans_names);

    UNPROTECT(4);

    return (ans);
}
