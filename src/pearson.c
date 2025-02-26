#include "uri.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains code copied from Dan Mahoney
**  for evaluation of spread options
**  using the "Pearson" method.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* Inner integration for spread options (from Pearson 1995) */
#define NDist NORM_DIST
static double Inner(double S, double S1, double S2, double sig1, double sig2,
                    double rho, double r, double q1, double q2, double tau,
                    double X)
{
    const double A1  = r * (1.0 - rho * sig2 / sig1);
    const double A2  = q2 - q1 * rho * sig2 / sig1;
    const double A3  = 0.5 * rho * sig2 * (sig1 - rho * sig2);
    const double A   = (A1 - A2 + A3) * tau;
    const double sig = sig2 * sqrt((1.0 - rho * rho) * tau);
    const double m1  = log(S1) + (r - q1 - 0.5 * sig1 * sig1) * tau;
    const double m2  = log(S2) + (r - q2 - 0.5 * sig2 * sig2) * tau;
    const double M2  = m2 + (rho * sig2 / sig1) * (log(S) - m1);
    double       x1 = 0, x2 = 0;

    if (S + X > 0)
    {
        if (1 == rho * rho)
            return (M2 > log(S + X)
                        ? exp(A) * S2 * pow(S / S1, rho * sig2 / sig1) - (S + X)
                        : 0);
        else
        {
            x1 = (M2 + sig * sig - log(S + X)) / sig;
            x2 = x1 - sig;
            return exp(A) * S2 * pow(S / S1, rho * sig2 / sig1) * NDist(x1) -
                   (S + X) * NDist(x2);
        }
    } else
        return exp(A) * S2 * pow(S / S1, rho * sig2 / sig1) - (S + X);
}

static double OneStepGridSpread(long nGrdPts, double S1, double S2, double vol1,
                                double vol2, double rho, double r, double q1,
                                double q2, double X, double T)
{
    const long   N = nGrdPts, N2 = 2 * N;
    const double sig1 = vol1 > 0 ? vol1 : 1e-10;
    const double sig2 = vol2 > 0 ? vol2 : 1e-10;
    const double sigT = sig1 * sqrt(T);
    const double h    = 5 * sigT / N;
    const double w    = sigT / h;
    const double m    = log(S1) + T * (r - q1 - 0.5 * sig1 * sig1);
    const double c    = rho * sig2;
    const double A    = T * (r * (1.0 - c / sig1) - (q2 - q1 * c / sig1) +
                          0.5 * rho * sig2 * (sig1 - c));

    double u1, u2, v1, v2;
    double z1, z2, x1, x2;
    double a, b, aN, bN, a0, b0, term1, term2, term3, term4;
    long   i;
    double value = 0;

    if (0.0 >= T) return 0.0 > T ? NA_REAL : MAX(S2 - S1, 0.0);

    if (vol1 < 0 || vol2 < 0 || S1 < 0 || S2 < 0 || rho < -1 || 1 < rho ||
        N <= 0)
        return NA_REAL;

    x2 = exp(m - N * h);
    z2 = Inner(x2, S1, S2, sig1, sig2, rho, r, q1, q2, T, X);
    u2 = NORM_DIST(-N / w - sigT);
    v2 = NORM_DIST(-N / w);

    for (i = -N; i < N; i++)
    {
        x1 = x2;
        z1 = z2;
        u1 = u2;
        v1 = v2;

        x2 = exp(m + h * (i + 1));
        z2 = Inner(x2, S1, S2, sig1, sig2, rho, r, q1, q2, T, X);
        u2 = NORM_DIST((i + 1) / w - sigT);
        v2 = NORM_DIST((i + 1) / w);

        a = (z2 - z1) / (x2 - x1);
        b = (x2 * z1 - x1 * z2) / (x2 - x1);

        value += a * (u2 - u1) * S1 * exp((r - q1) * T) + b * (v2 - v1);
    }

    x1 = exp(m - N * h);

    aN = exp(A) * S2 / pow(S1, c / sig1);
    bN = z2 - aN * pow(x2, c / sig1);
    a0 = exp(A) * S2 / pow(S1, c / sig1);
    b0 = Inner(x1, S1, S2, sig1, sig2, rho, r, q1, q2, T, X) -
         a0 * pow(x1, c / sig1);

    term1 = a0 * NORM_DIST(-N / w - c * sqrt(T)) *
            exp(c * m / sig1 + 0.5 * c * c * T);
    term2 = b0 * NORM_DIST(-N / w);
    term3 = aN * (1.0 - NORM_DIST(N / w - c * sqrt(T))) *
            exp(c * m / sig1 + 0.5 * c * c * T);
    term4 = bN * (1.0 - NORM_DIST(N / w));

    value += term1 + term2 + term3 + term4;

    return exp(-r * T) * value;
}

static void margrabeApprox(double tau, double S1, double S2, double K,
                           double vol1, double vol2, double rho, double r,
                           double q1, double q2, char optionType,
                           int calcDeltas, double *valG)
{
    /*
    ** There are some errors in the way this is currently implemented.
    ** For now, simply set valG[i] = NA.
    */

    double tmp[MARGRABE0_PTR_LEN] = {0.0};

    if (!valG) return;

    *valG = NA_REAL;
    if (calcDeltas)
        valG[PEARSON0_PRICE_INDEX]  = valG[PEARSON0_D1_INDEX] =
            valG[PEARSON0_D2_INDEX] = NA_REAL;

    if (0.0 == K)
    {
        margrabe0(tau, S1, S2, vol1, vol2, rho, r, q1, q2, optionType, tmp);

        valG[PEARSON0_PRICE_INDEX] = tmp[MARGRABE0_PRICE_INDEX];

        if (calcDeltas)
        {
            valG[PEARSON0_D1_INDEX] = tmp[MARGRABE0_D1_INDEX];
            valG[PEARSON0_D2_INDEX] = tmp[MARGRABE0_D2_INDEX];
        }
    }
    /*
    ** Old code below.
    */

#undef OLD_CODE_BELOW
#ifdef OLD_CODE_BELOW

    double S3 = 0, vol3 = 0, rho2 = 0, a = 0, b = 0, c = 0, d = 0;
    /* For partial derivatives (p="prime"): */
    double ap = 0, bp = 0, cp = 0, dp = 0, vol3p = 0, rho2p = 0;

    if (!valG) error("null pointer passed to margrabeApprox");

    if (calcDeltas) valG[PEARSON0_D1_INDEX] = valG[PEARSON0_D2_INDEX] = 0;

    if (K > 0)
    {
        S3    = S1 + exp(-q1 * tau) * K;
        a     = S1 * S1 * exp(vol1 * vol1 * tau) + 2 * K * S1 + K * K;
        b     = S3 * S3;
        ap    = 2 * S1 * exp(vol1 * vol1 * tau) + 2 * K;
        bp    = 2 * S3;
        vol3  = sqrt(log(a / b) / tau);
        vol3p = (ap / a - bp / b) / (tau * 2 * vol3);
        c     = S1 * exp(rho * vol1 * vol2 * tau) + K;
        d     = S3;
        cp    = exp(rho * vol1 * vol2 * tau);
        dp    = 1;
        rho2  = log(c / d) / (vol2 * vol3 * tau);
        rho2p =
            -(vol3p / vol3) * rho2 + (cp / c - dp / d) / (vol2 * vol3 * tau);

        margrabe0(tau, S3, S2, vol3, vol2, rho2, r, q1, q2, optionType, tmp);

        if (calcDeltas)
            valG[PEARSON0_D1_INDEX] += tmp[MARGRABE0_V1_INDEX] * vol3p +
                                       tmp[MARGRABE0_ETA_INDEX] * rho2p;

    } else if (K < 0)
    {
        S3    = S2 - exp(-q2 * tau) * K;
        a     = S2 * S2 * exp(vol2 * vol2 * tau) - 2 * K * S2 + K * K;
        b     = S3 * S3;
        ap    = 2 * S2 * exp(vol2 * vol2 * tau) - 2 * K;
        bp    = 2 * S3;
        vol3  = sqrt(log(a / b) / tau);
        vol3p = (ap / a - bp / b) / (tau * 2 * vol3);
        c     = S2 * exp(rho * vol1 * vol2 * tau) - K;
        d     = S3;
        cp    = exp(rho * vol1 * vol2 * tau);
        dp    = 1;
        rho2  = log(c / d) / (vol1 * vol3 * tau);
        rho2p =
            -(vol3p / vol3) * rho2 + (cp / c - dp / d) / (vol2 * vol3 * tau);

        margrabe0(tau, S1, S3, vol1, vol3, rho2, r, q1, q2, optionType, tmp);

        if (calcDeltas)
            valG[PEARSON0_D2_INDEX] += tmp[MARGRABE0_V2_INDEX] * vol3p +
                                       tmp[MARGRABE0_ETA_INDEX] * rho2p;

    } else
        margrabe0(tau, S1, S2, vol1, vol2, rho, r, q1, q2, optionType, tmp);

    valG[PEARSON0_PRICE_INDEX] = tmp[MARGRABE0_PRICE_INDEX];

    if (calcDeltas)
    {
        valG[PEARSON0_D1_INDEX] += tmp[MARGRABE0_D1_INDEX];
        valG[PEARSON0_D2_INDEX] += tmp[MARGRABE0_D2_INDEX];
    }

#endif /* OLD_CODE_BELOW */
}

/*******************************************************************
 *
 *  NB in pearson0 tau is assumed to be in years (or whatever
 *    units vol1, vol2, r, q1, q2 are in).
 *
 *******************************************************************/

void pearson0(double tau, double S1, double S2, double K, double vol1,
              double vol2, double rho, double r, double q1, double q2,
              char optionType, long nGrdPts, int calcDeltas, double *valG)
{
    const double u = 0.01;
    double       valUp, valDn;

    if (!valG) error("Null pointer.");

    if (!R_FINITE(tau + S1 + S2 + K + vol1 + vol2 + rho + r + q1 + q2))
    {
        valG[PEARSON0_PRICE_INDEX] = NA_REAL;

        if (calcDeltas)
            valG[PEARSON0_D1_INDEX] = valG[PEARSON0_D2_INDEX] = NA_REAL;

        return;
    }

    if (0.0 == K || ABS(rho) > 1 || nGrdPts <= 0 || !R_FINITE((double)nGrdPts))
    {
        margrabeApprox(tau, S1, S2, K, vol1, vol2, rho, r, q1, q2, optionType,
                       calcDeltas, valG);
        return;
    }

    *valG =
        OneStepGridSpread(nGrdPts, S1, S2, vol1, vol2, rho, r, q1, q2, K, tau);

    if (calcDeltas)
    {
        valUp = OneStepGridSpread(nGrdPts, (1 + u) * S1, S2, vol1, vol2, rho, r,
                                  q1, q2, K, tau);
        valDn = OneStepGridSpread(nGrdPts, (1 - u) * S1, S2, vol1, vol2, rho, r,
                                  q1, q2, K, tau);

        valG[PEARSON0_D1_INDEX] = (valUp - valDn) / (2 * u * S1);

        valUp = OneStepGridSpread(nGrdPts, S1, (1 + u) * S2, vol1, vol2, rho, r,
                                  q1, q2, K, tau);
        valDn = OneStepGridSpread(nGrdPts, S1, (1 - u) * S2, vol1, vol2, rho, r,
                                  q1, q2, K, tau);

        valG[PEARSON0_D2_INDEX] = (valUp - valDn) / (2 * u * S2);
    }

    switch (optionType)
    {
        case 'c':
            break;
        case 'p':
            valG[PEARSON0_PRICE_INDEX] -=
                S2 * exp(-q2 * tau) - S1 * exp(-q1 * tau) - K * exp(-r * tau);
            if (calcDeltas)
            {
                valG[PEARSON0_D1_INDEX] += exp(-q1 * tau);
                valG[PEARSON0_D2_INDEX] -= exp(-q2 * tau);
            }
            break;
        case 's':
            valG[PEARSON0_PRICE_INDEX] *= 2;
            valG[PEARSON0_PRICE_INDEX] -=
                (S2 * exp(-q2 * tau) - S1 * exp(-q1 * tau) - K * exp(-r * tau));
            if (calcDeltas)
            {
                valG[PEARSON0_D1_INDEX] =
                    2 * valG[PEARSON0_D1_INDEX] + exp(-q1 * tau);
                valG[PEARSON0_D2_INDEX] =
                    2 * valG[PEARSON0_D2_INDEX] - exp(-q2 * tau);
            }
            break;
        default:
            valG[PEARSON0_PRICE_INDEX] = NA_REAL;
            if (calcDeltas)
                valG[PEARSON0_D1_INDEX] = valG[PEARSON0_D2_INDEX] = NA_REAL;
            break;
    }
}
#undef NDist

static void pearson_chkargs(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, SEXP *vol1,
                            SEXP *vol2, SEXP *rho, SEXP *r, SEXP *q1, SEXP *q2,
                            SEXP *opt, int *numprot)
{
#define N_NUMARGS 10
    int i, N = 0;

    SEXP *args[N_NUMARGS + 1]  = {tau, S1, S2, K,  vol1, vol2,
                                  rho, r,  q1, q2, opt};
    char *names[N_NUMARGS + 1] = {"tau", "S1", "S2", "K",  "vol1", "vol2",
                                  "rho", "r",  "q1", "q2", "opt"};
    /* Make sure you haven't passed any null pointers */
    for (i = 0; i < N_NUMARGS + 1; i++)
        if (0 == args[i]) error("Null pointer passed.");

    if (0 == numprot) error("Null pointer passed.");

    /* Coerce types  */
    for (i = 0; i < N_NUMARGS; i++)
        PROT2(*args[i] = AS_NUMERIC(*args[i]), *numprot);

    /* Find maximal length */
    for (i = 0, N = 0; i < N_NUMARGS + 1; i++) N = MAX(N, length(*args[i]));

    if (N > 1)
    { /* Repeat arguments of length == 1. */

        for (i = 0; i < N_NUMARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = numrep(*args[i], N), *numprot);

        if (1 == length(*opt)) PROT2(*opt = charrep(*opt, N), *numprot);
    }

    /* Check lengths */

    for (i = 0; i < N_NUMARGS + 1; i++)
        if (N != length(*args[i]))
            error("Argument %s of wrong length", names[i]);

#undef N_NUMARGS
}

SEXP pearson(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho,
             SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP calcDeltas, SEXP grdPts)
{
    int        numprot = 0;
    const int  calcDel = asInteger(calcDeltas);
    const int  retLen  = calcDel ? PEARSON0_PTR_LEN : 1;
    const long nGrdPts = asInteger(grdPts);
    long       i, N;
    SEXP       ans, dimnames, colnames;
    carray     ansMatrix;
    char      *names[PEARSON0_PTR_LEN] = PEARSON0_RET_NAMES;

    pearson_chkargs(&tau, &S1, &S2, &K, &vol1, &vol2, &rho, &r, &q1, &q2, &opt,
                    &numprot);

    N = length(tau); /* By now all arguments should have same length */

    ansMatrix = make_zero_matrix(N, retLen);

    for (i = 0; i < N; i++)
    {
        pearson0(REAL(tau)[i] / 365.0, REAL(S1)[i], REAL(S2)[i], REAL(K)[i],
                 REAL(vol1)[i], REAL(vol2)[i], REAL(rho)[i], REAL(r)[i],
                 REAL(q1)[i], REAL(q2)[i], *CHAR(GET_ELT(opt, i)), nGrdPts,
                 /* calcDeltas = */ calcDel, ARRAY2(ansMatrix)[i]);
    }

    PROT2(ans = carray_to_sexp(ansMatrix), numprot);
    PROT2(colnames = NEW_STRING(retLen), numprot);
    PROT2(dimnames = NEW_LIST(2), numprot);

    for (i = 0; i < retLen; i++) SET_ELT(colnames, i, mkChar(names[i]));

    SET_ELT(dimnames, 0, GET_NAMES(tau));
    SET_ELT(dimnames, 1, colnames);

    setAttrib(ans, R_DimNamesSymbol, dimnames);

    UNPROTECT(numprot);

    return ans;
}

/*******************************************************************
 *
 *  Description: Arguments are tau, S1, S2, K, vol1, vol2, rho,
 *    r, q1, q1, opt are lists where each element is of
 *    the proper type for function pearson above.  Returns list
 *    where each element is the output from applying pearson()
 *    to the corresponding arguments.
 *
 *******************************************************************/

#define NARGS 11
static void listPearson_chk(SEXP *tau, SEXP *S1, SEXP *S2, SEXP *K, SEXP *vol1,
                            SEXP *vol2, SEXP *rho, SEXP *r, SEXP *q1, SEXP *q2,
                            SEXP *opt, int *numprot)
{
    long  i, N = 0;
    SEXP *args[NARGS] = {tau, S1, S2, K, vol1, vol2, rho, r, q1, q2, opt};

    for (i = 0; i < NARGS; i++)
    {
        if (!isNewList(*args[i])) PROT2(*args[i] = tolist(*args[i]), *numprot);

        N = MAX(N, length(*args[i]));
    }

    if (1 < N)
        for (i = 0; i < NARGS; i++)
            if (1 == length(*args[i]))
                PROT2(*args[i] = listrep(*args[i], N), *numprot);
}
#undef NARGS

SEXP listPearson(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2,
                 SEXP rho, SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP calcDeltas,
                 SEXP grdPts)
{
    int  numprot = 0;
    long i, N;

    SEXP ans, tmp;

    listPearson_chk(&tau, &S1, &S2, &K, &vol1, &vol2, &rho, &r, &q1, &q2, &opt,
                    &numprot);

    N = length(tau);
    PROT2(ans = NEW_LIST(N), numprot);

    for (i = 0; i < N; i++)
    {
        tmp = pearson(GET_ELT(tau, i), GET_ELT(S1, i), GET_ELT(S2, i),
                      GET_ELT(K, i), GET_ELT(vol1, i), GET_ELT(vol2, i),
                      GET_ELT(rho, i), GET_ELT(r, i), GET_ELT(q1, i),
                      GET_ELT(q2, i), GET_ELT(opt, i), calcDeltas, grdPts);

        SET_ELT(ans, i, tmp);
    }

    SET_NAMES(ans, GET_NAMES(tau));

    UNPROT2;

    return ans;
}
