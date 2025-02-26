
#include "multi.h"

#define STATE_LEN 9
static void algo3_chk(SEXP *currentState, int *numprot)
{
    int   i;
    char *names[STATE_LEN] = {"alpha", "delta",  "sigv", "m", "M",
                              "N",     "lambda", "h",    "r"};
    SEXP  args[STATE_LEN];

    for (i = 0; i < STATE_LEN; i++)
        PROTECT(args[i] = getListElt(*currentState, names[i]));

    for (i = 0; i < STATE_LEN; i++)
        if (!IS_NUMERIC(args[i]))
        {
            PROT2(args[i] = AS_NUMERIC(args[i]), *numprot);
            SET_ELT(*currentState, i, args[i]);
        }

    if (length(getListElt(*currentState, "h")) !=
        length(getListElt(*currentState, "r")))
        error("length(h) != length(r)");

    UNPROTECT(STATE_LEN);
}

#define RIGAMMA(N, lambda) (1 / rgamma((N) / 2.0, 2.0 / (lambda)))
#define DIGAMMA(x, phi, theta) \
    (dgamma(1.0 / (x), (phi) + 2.0, 1.0 / (theta), 0))

#define STATE_LEN     9
#define NUM_INIT_PARS STATE_LEN - 1
static SEXP getInitState(SEXP r, SEXP initState, int diffusePrior)
{
    int     numprot        = 0;
    int     initStateGiven = 0;
    long    T, t, i;
    double  alpha, delta, sigv, N, lambda;
    double *hPtr, *loghPtr;
    carray  m, M;
    SEXP    h, logh, rNew, initVals[NUM_INIT_PARS];
    char   *names[STATE_LEN] = {"alpha", "delta",  "sigv", "m", "M",
                                "N",     "lambda", "h",    "r"};

    /* Replace r by (NA, r, NA) */
    PROT2(r = AS_NUMERIC(r), numprot);
    T = length(r);
    PROT2(rNew = NEW_NUMERIC(T + 2), numprot);
    REAL(rNew)[0] = REAL(rNew)[T + 1] = NA_REAL;
    memcpy(REAL(rNew) + 1, REAL(r), T * sizeof(double));
    r = rNew;
    T = length(r);

    PROT2(h = NEW_NUMERIC(T), numprot);
    PROT2(logh = NEW_NUMERIC(T), numprot);
    hPtr    = REAL(h);
    loghPtr = REAL(logh);

    GetRNGstate();

    /* See if initial parameters given */
    if (!isNull(initState))
    {
        for (i = 0, initStateGiven = 1; i < NUM_INIT_PARS; i++)
        {
            PROTECT(initVals[i] = AS_NUMERIC(getListElt(initState, names[i])));
            numprot        += 1;
            initStateGiven *= !isNull(initVals[i]);
        }

        if (!initStateGiven)
            warning(
                "some initial given value was null, using uninformed priors");
    }

    if (initStateGiven)
    {
        alpha  = REAL(initVals[0])[length(initVals[0]) - 1];
        delta  = REAL(initVals[1])[length(initVals[1]) - 1];
        sigv   = REAL(initVals[2])[length(initVals[2]) - 1];
        m      = sexp_to_carray(initVals[3], 1);
        M      = sexp_to_carray(initVals[4], 1);
        N      = asReal(initVals[5]);
        lambda = asReal(initVals[6]);
        PROT2(h = initVals[7], numprot);

        if (!is_matrix(m) || !is_matrix(M))
            error("Non-matrices provided for m, M");
        if (2 != NROW(m) || 2 != NROW(M) || 1 != NCOL(m) || 2 != NCOL(M))
            error("M, m of wrong size");
        if (length(h) != T)
            error("given conditional variances have wrong length");

    } else
    {
        /* Draw alpha, delta, sigv from priors */
        N      = 4;
        lambda = 0.005;
        sigv   = sqrt(RIGAMMA(N, lambda));

        alpha = norm_rand();
        delta = unif_rand();
        m     = make_zero_matrix(2, 1);
        M     = make_identity_matrix(2);
        scalar_op(M, sigv * sigv, '/', M);
        ARRAY1(m)[1]  = 0.5;
        ARRAY1(M)[3] /= 12;

        /* Simulate log(h) */
        loghPtr[0] = rnorm(/*mu=*/alpha / (1 - delta),
                           /*sigma=*/sigv / sqrt(1 - delta * delta));

        for (t = 1; t < T; t++)
            loghPtr[t] = delta * loghPtr[t - 1] + rnorm(alpha, sigv);

        for (t = 0; t < T; t++) hPtr[t] = exp(loghPtr[t]);
    }

    PutRNGstate();

    PROT2(initState = NEW_LIST(STATE_LEN), numprot);

    SET_ELT(initState, 0, ScalarReal(alpha));
    SET_ELT(initState, 1, ScalarReal(delta));
    SET_ELT(initState, 2, ScalarReal(sigv));
    SET_ELT(initState, 3, carray_to_sexp(m));
    SET_ELT(initState, 4, carray_to_sexp(M));
    SET_ELT(initState, 5, ScalarReal(N));
    SET_ELT(initState, 6, ScalarReal(lambda));
    SET_ELT(initState, 7, h);
    SET_ELT(initState, 8, r);

    set_names(initState, names);

    UNPROTECT(numprot);

    return initState;
}
#undef NUM_INIT_PARS

#define PDENS(h, r, mu, sigma)                 \
    (R_pow(h, -1.5) * exp(-0.5 * SQR(r) / h) * \
     exp(-0.5 * SQR(log(h) - mu) / SQR(sigma)))

#define ANSLEN 9
#define MAX_IT 10
SEXP algo3(SEXP r, SEXP info, SEXP initState)
{
    int     numprot = 0;
    int     maxit, drawSuccess, diffusePrior;
    long    i, it, t, numLoops, initLoops, T, n;
    double *hPtr = NULL, *rPtr = NULL;
    double alpha, delta, sigv, N, lambda, hnew, c, phi, mu, theta, qMode, sigma;
    double rho, sig1, sig2, e1, e2;
    double p, q, pnew, qnew, fnew, f;
    char  *ansNames[ANSLEN] = {"alpha", "delta",  "sigv", "m",    "M",
                               "N",     "lambda", "h",    "hTerm"};
    carray m, M;
    carray y, X, Xt, Mi, Mnew, res, b, Xb, XtX, XtXi, Xty;
    carray tmp_2x2, tmp_1x2, tmp_2x1, tmp_1x1;
    SEXP   ans, alphaOut, deltaOut, sigvOut, hOut;
    SEXP   newState, h;

    diffusePrior = asInteger(getListElt(info, "diffusePrior"));

    PROT2(initState = getInitState(r, initState, diffusePrior), numprot);

    if (!isNewList(initState) || STATE_LEN != length(initState))
        error("initState must be a list");

    algo3_chk(&initState, &numprot);

    numLoops  = asInteger(getListElt(info, "numLoops"));
    initLoops = asInteger(getListElt(info, "initLoops"));
    maxit     = asInteger(getListElt(info, "maxit"));
    alpha     = asReal(getListElt(initState, "alpha"));
    delta     = asReal(getListElt(initState, "delta"));
    sigv      = asReal(getListElt(initState, "sigv"));
    m         = sexp_to_carray(getListElt(initState, "m"), 1);
    M         = sexp_to_carray(getListElt(initState, "M"), 1);
    N         = asReal(getListElt(initState, "N"));
    lambda    = asReal(getListElt(initState, "lambda"));
    PROT2(h = getListElt(initState, "h"), numprot);
    PROT2(r = getListElt(initState, "r"), numprot);
    T    = length(h);
    hPtr = REAL(h);
    rPtr = REAL(r);

    PROT2(alphaOut = NEW_NUMERIC(numLoops), numprot);
    PROT2(deltaOut = NEW_NUMERIC(numLoops), numprot);
    PROT2(sigvOut = NEW_NUMERIC(numLoops), numprot);
    PROT2(hOut = NEW_NUMERIC(numLoops), numprot);
    PROT2(ans = NEW_LIST(ANSLEN), numprot);
    PROT2(newState = NEW_LIST(STATE_LEN), numprot);

    y       = make_zero_matrix(T - 1, 1);
    res     = make_zero_matrix(T - 1, 1);
    X       = make_zero_matrix(T - 1, 2);
    Xb      = make_zero_matrix(T - 1, 1);
    Xt      = make_zero_matrix(2, T - 1);
    Mi      = make_zero_matrix(2, 2);
    Mnew    = make_zero_matrix(2, 2);
    b       = make_zero_matrix(2, 1);
    XtX     = make_zero_matrix(2, 2);
    XtXi    = make_zero_matrix(2, 2);
    Xty     = make_zero_matrix(2, 1);
    tmp_2x2 = make_zero_matrix(2, 2);
    tmp_1x2 = make_zero_matrix(1, 2);
    tmp_2x1 = make_zero_matrix(2, 1);
    tmp_1x1 = make_zero_matrix(1, 1);

    for (t = 0; t < T - 1; t++) ARRAY2(Xt)[0][t] = 1;

    GetRNGstate();

    for (n = -initLoops; n < numLoops; n++)
    {
        /*
        **  Update h
        */

        sigma = sigv / sqrt(1 + delta * delta);
        phi   = exp(sigv * sigv) / expm1(sigv * sigv) + 1 + 0.5;

        for (i = 1; i < T - 1; i++)
        {
            mu = alpha * (1 - delta) +
                 delta * (log(hPtr[i - 1]) + log(hPtr[i + 1]));
            mu /= (1 + delta * delta);
            theta = 0.5 * SQR(rPtr[i]) + (phi - 1) * exp(mu + 0.5 * SQR(sigma));

            hnew = RIGAMMA(2 * phi, theta / 2);

            pnew = PDENS(hnew, rPtr[i], mu, sigma);
            /* pnew  = R_pow(hnew, -1.5) * exp(-0.5*SQR(rPtr[i])/hnew); */
            /* pnew *= exp(-0.5 * SQR(log(hnew) - mu) / SQR(sigma)); */
            qnew = DIGAMMA(hnew, phi, theta);
            p    = PDENS(hPtr[i], rPtr[i], mu, sigma);
            /* p = R_pow(hPtr[i], -1.5) * exp(-0.5*SQR(rPtr[i])/hPtr[i]); */
            /* p *= exp(-0.5 * SQR(log(hPtr[i]) - mu) / SQR(sigma)); */
            q = DIGAMMA(hPtr[i], phi, theta);

            hPtr[i] = unif_rand() < (pnew / qnew) / (p / q) ? hnew : hPtr[i];
        }

        hPtr[0] = exp(alpha + delta * log(hPtr[1]) + sigv * norm_rand());
        hPtr[T - 1] =
            exp(alpha + delta * log(hPtr[T - 2]) + sigv * norm_rand());

        /*
        **  Update alpha, delta, sigv
        */

        for (t = 0; t < T - 1; t++)
        {
            ARRAY2(Xt)[1][t] = log(hPtr[t]);
            ARRAY1(y)[t]     = log(hPtr[t + 1]);
        }

        transpose_matrix(Xt, X);
        matrix_prod_ns(X, X, 1, 0, XtX);
        matrix_inverse_2x2(XtX, XtXi);
        matrix_prod_ns(X, y, 1, 0, Xty);
        matrix_prod_ns(XtXi, Xty, 0, 0, b);
        matrix_prod_ns(X, b, 0, 0, Xb);

        /* Update N and lambda */
        N = N + T - 1;

        array_op(y, Xb, '-', res);

        if (!matrix_inverse_2x2(M, Mi)) set_array_to_zero(Mi);

        array_op(XtX, Mi, '+', Mnew);            /* Mnew <- XtX + M^-1 */
        matrix_inverse_2x2(Mnew, Mnew);          /* Mnew <- (XtX + M^-1)^-1 */
        matrix_prod_ns(Mi, Mnew, 0, 0, tmp_2x2); /* tmp_2x2 <- Mi * Mnew*/
        matrix_prod(tmp_2x2, XtX, 0, 0,
                    tmp_2x2); /* tmp_2x2 <- Mi * Mnew * X'X */
        array_op(b, m, '-', tmp_2x1);
        /* tmp_1x2 <- (b-m)' * Mi * Mnew * X'X : */
        matrix_prod_ns(tmp_2x1, tmp_2x2, 1, 0, tmp_1x2);
        /* tmp_1x1 <- (b-m)' * Mi * Mnew * X'X * (b-m) : */
        matrix_prod_ns(tmp_1x2, tmp_2x1, 0, 0, tmp_1x1);

        lambda += sumsq(res) + ARRAY1(tmp_1x1)[0];

        /* Draw new sigv */
        sigv = sqrt(RIGAMMA(N, lambda));

        /* Draw new alpha, delta conditional on new sigma.v */
        matrix_prod(Mi, m, 0, 0, m);
        array_op(m, Xty, '+', m);
        matrix_prod(Mnew, m, 0, 0, m);

        for (i = 0; i < 4; i++) ARRAY1(M)[i] = ARRAY1(Mnew)[i];

        if (ABS(ARRAY1(M)[1] - ARRAY1(M)[2]) > DOUBLE_EPS)
            error("M not symmetric");

        if (ARRAY1(M)[0] < 0 || ARRAY1(M)[3] < 0)
            error("Negative variances in M");

        matrix_inverse_2x2(M, Mi);

        /* Draw new alpha, delta conditional on new sigv */

        sig1 = sigv * sqrt(ARRAY1(M)[0]);
        sig2 = sigv * sqrt(ARRAY1(M)[3]);
        rho  = ARRAY1(M)[1] / sqrt(ARRAY1(M)[0] * ARRAY1(M)[3]);
        it   = 0;

        do
        {
            e1 = norm_rand();
            e2 = rho * e1 + sqrt(1 - rho * rho) * norm_rand();

            alpha = ARRAY1(m)[0] + sig1 * e1;
            delta = ARRAY1(m)[1] + sig2 * e2;

            it += 1;

            /* Continue until |delta| < 1 for stationarity */

        } while (ABS(delta) >= 1 && it < MAX_IT);

        if (ABS(delta) >= 1 || !R_FINITE(alpha + delta))
            error("stationary draw of delta unsuccessful");

        if (n >= 0)
        {
            REAL(alphaOut)[n] = alpha;
            REAL(deltaOut)[n] = delta;
            REAL(sigvOut)[n]  = sigv;
            REAL(hOut)[n]     = hPtr[T - 2];
        }
    }

    PutRNGstate();

    SET_ELT(ans, 0, alphaOut);
    SET_ELT(ans, 1, deltaOut);
    SET_ELT(ans, 2, sigvOut);
    SET_ELT(ans, 3, carray_to_sexp(m));
    SET_ELT(ans, 4, carray_to_sexp(M));
    SET_ELT(ans, 5, ScalarReal(N));
    SET_ELT(ans, 6, ScalarReal(lambda));
    SET_ELT(ans, 7, h);
    SET_ELT(ans, 8, hOut);

    set_names(ans, ansNames);

    UNPROT2;

    return ans;
}
#undef STATE_LEN
#undef ANSLEN
#undef MAX_IT

#undef RIGAMMA
#undef DIGAMMA
