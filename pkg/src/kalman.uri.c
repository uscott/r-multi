
#include "multi.h"

#define N_SEXP_ARGS 7
static void kalman1_chkargs(SEXP *y1, SEXP *H1tr, SEXP *F1, SEXP *Q1, SEXP *R1,
                            SEXP *x1, SEXP *A1tr, int *numprot)
{
    int  i, n, r, ok, k;
    long T;
    SEXP tmp;

    SEXP *args[N_SEXP_ARGS]  = {y1, H1tr, F1, Q1, R1, x1, A1tr};
    char *names[N_SEXP_ARGS] = {"y1", "H1tr", "F1", "Q1", "R1", "x1", "A1tr"};

    if (!isMatrix(*y1))
    {
        PROTECT(tmp = AS_NUMERIC(*y1));

        T = length(tmp);

        PROT2(*y1 = allocMatrix(REALSXP, T, 1), *numprot);

        memcpy(REAL(*y1), REAL(tmp), T * sizeof(double));

        UNPROTECT_PTR(tmp);
    }

    T = nrows(*y1);
    n = ncols(*y1);
    r = ncols(*H1tr);

    if (isNull(*x1) || length(*x1) <= 0)
    {
        k = 1;

        PROT2(*x1 = allocMatrix(REALSXP, T, k), *numprot);
        PROT2(*A1tr = allocMatrix(REALSXP, n, k), *numprot);

        dsetzero(REAL(*x1), T * k);
        dsetzero(REAL(*A1tr), n * k);
    }

    k = ncols(*x1);

    for (i = 0; i < N_SEXP_ARGS; i++)
        if (!isMatrix(*args[i])) error("%s is not a matrix", names[i]);

    ok = r == nrows(*F1) && r == ncols(*F1) && n == nrows(*H1tr) &&
         r == nrows(*Q1) && r == ncols(*Q1) && n == nrows(*R1) &&
         n == ncols(*R1) && n == nrows(*A1tr) && k == ncols(*A1tr);

    if (!ok) error("Some argument does not have the proper dimension");
}
#undef N_SEXP_ARGS

static void kalman_recursion(carray y, carray Htr, carray F, carray Q, carray R,
                             carray Atr, carray x, carray xi, carray yfit,
                             carray res, double *ll)
{
    long   t;
    int    r, k, n, ok;
    long   T;
    carray H, A, P, K;
    /* For use as "generic" matrices: */
    carray rxr, nxn_a, nxn_b, rx1_a, rx1_b, nx1_a, nx1_b, kx1, r2xr2, r2x1;

    ok = is_matrix(y) && is_matrix(Htr) && is_matrix(F) && is_matrix(Q) &&
         is_matrix(R) && is_matrix(Atr) && is_matrix(x) && is_matrix(xi) &&
         is_matrix(yfit) && is_matrix(res);

    ok *= (ll != NULL);

    if (!ok) error("bad args in kalman_recursion");

    r = NCOL(Htr);
    k = NCOL(x);
    n = NCOL(y);
    T = NROW(y);

    H = make_zero_matrix(r, n); /* For transpose of Htr     */
    A = make_zero_matrix(k, n); /* For transpose of Atr     */

    P = make_zero_matrix(r, r); /* Covariance matrix        */
    K = make_zero_matrix(r, n); /* Gain matrix              */

    rxr   = make_zero_matrix(r, r); /* Generic r x r matrix     */
    nxn_a = make_zero_matrix(n, n); /* Generic n x n matrix     */
    nxn_b = make_zero_matrix(n, n); /* Generic n x n matrix     */
    rx1_a = make_zero_matrix(r, 1); /* Generic r x 1 matrix     */
    rx1_b = make_zero_matrix(r, 1); /* Generic r x 1 matrix     */
    nx1_a = make_zero_matrix(n, 1); /* Generic n x 1 matrix     */
    nx1_b = make_zero_matrix(n, 1); /* Generic n x 1 matrix     */
    kx1   = make_zero_matrix(k, 1); /* Generic k x 1 matrix     */

    r2xr2 = make_zero_matrix(r * r, r * r); /* Generic r^2 x r^2 matrix */
    r2x1  = make_zero_matrix(r * r, 1);     /* Generic r^2 x 1   matrix */

    transpose_matrix(Htr, H);
    transpose_matrix(Atr, A);

    kronecker(F, F, r2xr2); /* r2xr2 <- kronecker product of F w/ itself */

    array_op(make_identity_matrix(r * r), r2xr2, '-', r2xr2);

    matrix_prod(r2xr2, vec(Q), 0, 0, r2x1);

    memcpy(ARRAY1(P),
           ARRAY1(transpose_matrix_b(make_matrix(ARRAY1(r2x1), r, r))),
           r * r * sizeof(double));

    matrix_inverse(P, P); /* P <- P(1|0) */

    *ll = n * T * log(M_2PI);

    for (t = 1; t < T && R_FINITE(*ll); t++)
    {
        memcpy(ARRAY1(rx1_b), ARRAY1(rx1_a), r * sizeof(double));
        memcpy(ARRAY1(kx1), MATRIX(x)[t - 1], k * sizeof(double));
        memcpy(ARRAY1(nx1_a), MATRIX(y)[t - 1], n * sizeof(double));
        /*
        **  rx1_b <- xi(t|t-1) (= 0 if t == 1)
        **  kx1   <- x[t, ]
        **  nx1_a <- y[t, ]
        */

        matrix_prod_3(H, P, H, 1, 0, 0, nxn_a);
        /*
        **  nxn_a <- H'P(t|t-1)H
        */

        array_op(nxn_a, R, '+', nxn_a);
        /*
        **  nxn_a <- H'P(t|t-1)H + R
        */

        *ll += log(det(nxn_a));
        /*
        **  *ll += log det (H'P(t|t-1)H + R)
        */

        matrix_inverse(nxn_a, nxn_b);
        /*
        **  nxn_b <- (H'P(t|t-1)H + R)^(-1)
        */

        matrix_prod_4(F, P, H, nxn_b, 0, 0, 0, 0, K);
        /*
        **  K <- FPH(H'P(t|t-1)H + R)^(-1)
        **      = gain matrix at time t
        */

        matrix_prod(F, rx1_a, 0, 0, rx1_a);
        /*
        **  rx1_a <- F * xi(t|t-1)
        */

        array_op(nx1_a, matrix_prod_b(Atr, kx1, 0, 0), '-', nx1_a);
        array_op(nx1_a, matrix_prod_b(Htr, rx1_b, 0, 0), '-', nx1_a);
        /*
        **  nx1_a <- y(t) - A'x(t) - H'xi(t|t - 1)
        **          = residual at time t
        */

        *ll += *ARRAY1(matrix_prod_3b(nx1_a, nxn_b, nx1_a, 1, 0, 0));
        /*
        **  *ll += z(t)' (H'P(t|t-1)H + R)^(-1) z(t)
        */

        memcpy(MATRIX(res)[t - 1], ARRAY1(nx1_a), n * sizeof(double));
        /*
        **  MATRIX(res)[t - 1] <- residual at time t
        */
        matrix_prod(K, nx1_a, 0, 0, rx1_b);

        array_op(rx1_a, rx1_b, '+', rx1_a);
        /*
        **  rx1_a <- xi(t+1|t)
        */

        memcpy(MATRIX(xi)[t], ARRAY1(rx1_a), r * sizeof(double));
        /*
        **  MATRIX(xi)[t] <- xi(t+1|t)
        */
        matrix_prod(K, Htr, 0, 0, rxr);
        array_op(F, rxr, '-', rxr);
        /*
        **  rxr <- F - KH'
        */

        matrix_prod_3(rxr, P, rxr, 0, 0, 1, P);
        /*
        **  P <- (F - KH')P(t|t-1)(F - KH')'
        */

        matrix_prod_3(K, R, K, 0, 0, 1, rxr);
        /*
        **  rxr <- KRK'
        */

        array_op(P, rxr, '+', P);
        array_op(P, Q, '+', P);
        /*
        **  P <- P(t+1|t)
        */
    }
    /*
    **    Time T stuff:
    */

    memcpy(ARRAY1(rx1_b), ARRAY1(rx1_a), r * sizeof(double));
    memcpy(ARRAY1(kx1), MATRIX(x)[t - 1], k * sizeof(double));
    memcpy(ARRAY1(nx1_a), MATRIX(y)[t - 1], n * sizeof(double));

    matrix_prod_3(H, P, H, 1, 0, 0, nxn_a);
    /*
    **  nxn_a <- H'P(t|t-1)H
    */

    array_op(nxn_a, R, '+', nxn_a);
    /*
    **  nxn_a <- H'P(t|t-1)H + R
    */

    *ll += log(det(nxn_a));
    /*
    **  *ll += log det (H'P(T|T-1)H + R)
    */

    matrix_inverse(nxn_a, nxn_b);
    /*
    **  nxn_b <- (H'P(T|T-1)H + R)^(-1)
    */

    array_op(nx1_a, matrix_prod_b(Atr, kx1, 0, 0), '-', nx1_a);
    array_op(nx1_a, matrix_prod_b(Htr, rx1_a, 0, 0), '-', nx1_a);

    *ll += *ARRAY1(matrix_prod_3b(nx1_a, nxn_b, nx1_a, 1, 0, 0));
    /*
    **  *ll += z(T)' (H'P(T|T-1)H + R)^(-1) z(T)
    */

    memcpy(MATRIX(res)[t - 1], ARRAY1(nx1_a), n * sizeof(double));

    *ll *= -0.5;
}

SEXP kalman1(SEXP y1, SEXP H1tr, SEXP F1, SEXP Q1, SEXP R1, SEXP x1, SEXP A1tr)

/*******************************************************************
 *
 *  Description: Kalman Filter with constant structural coeficients.
 *    Based on treatment in "Time Series Analysis" by James Hamilton,
 *    chapter 13.
 *
 *    The state-space representation is assumed to be given by
 *
 *    xi[t + 1, ] = F %*% xi[t, ] + v[t + 1, ]
 *    y[t, ]      = t(A) %*% x[t, ] + t(H) * xi[t, ] + w[t, ]
 *
 *  Parameters:
 *    y1
 *    - Matrix representing a multivariate time series y.
 *    H1tr
 *    - Matrix representing t(H)
 *    F1
 *    - Matrix representing F.
 *    Q1
 *    - Covariance matrix of v.
 *    R1
 *    - Covariance matrix of w.
 *    x1
 *    - Matrix of endogonous variable observations.
 *    A1tr
 *    - Represents t(A).
 *
 *******************************************************************/

{
    SEXP   ans, ans_names;
    int    numprot = 0;
    int    n, r, k = 0;
    long   T;
    double ll = 0.0;
    carray y, yfit, Htr, F, Q, R, x, Atr, xi, res;

    kalman1_chkargs(&y1, &H1tr, &F1, &Q1, &R1, &x1, &A1tr, &numprot);

    r = ncols(H1tr);
    k = ncols(x1);
    n = ncols(y1);
    T = nrows(y1);

    y   = sexp_to_carray(y1, 1); /* 2nd arg == 1 -> 1st arg is duplicated */
    Htr = sexp_to_carray(H1tr, 1);
    F   = sexp_to_carray(F1, 1);
    Q   = sexp_to_carray(Q1, 1);
    R   = sexp_to_carray(R1, 1);
    x   = sexp_to_carray(x1, 1);
    Atr = sexp_to_carray(A1tr, 1);

    yfit = make_zero_matrix(T, n); /* fitted values of y       */
    xi   = make_zero_matrix(T, r); /* estimated states         */
    res  = make_zero_matrix(T, n); /* Matrix of residuals      */

    /* The actual run */
    kalman_recursion(y, Htr, F, Q, R, Atr, x, xi, yfit, res, &ll);

    PROT2(ans_names = NEW_STRING(4), numprot);
    PROT2(ans = NEW_LIST(4), numprot);

    CHAR_PTR(ans_names)[0] = mkChar("xi.fitted");
    CHAR_PTR(ans_names)[1] = mkChar("y.fitted");
    CHAR_PTR(ans_names)[2] = mkChar("residuals");
    CHAR_PTR(ans_names)[3] = mkChar("ll");

    SET_ELT(ans, 0, carray_to_sexp(xi));
    SET_ELT(ans, 1, carray_to_sexp(array_op_b(y, res, '-')));
    SET_ELT(ans, 2, carray_to_sexp(res));
    SET_ELT(ans, 3, ScalarReal(ll));

    SET_NAMES(ans, ans_names);

    UNPROTECT(numprot);

    return ans;
}

#define NPAR 8
SEXP kal1(SEXP y, SEXP par)

/*******************************************************************
 *
 *  Description:
 *    Computes residuals and loglikelihood for model
 *
 *    y[t]  = mu[t] + sqrt(h[t]) * z[t]
 *    h[t]  = a0 + a1 * h[t-1] * (z[t-1] - gamma)^2 + b1 * h[t-1]
 *    mu[t] = phi * mu[t-1] + c + sigma * eta[t]
 *
 *    Starting values of h and mu are assumed to be
 *        E(h)  = a0 / (1 - b1 - a1 * (1 + gamma^2))
 *        E(mu) = c  / (1 - phi)
 *
 *  Parameters:
 *    y
 *    - univariate time series.
 *    par
 *    - vector of length NPAR with elements corresponding to
 *    a0, a1, b1, gamma, phi, c, sigma, rho.  rho is the
 *    correlation between z[t] and eta[t].
 *
 *******************************************************************/

{
    SEXP            ans, h, m, res;
    long            t, T;
    int             ok;
    int             numprot = 0;
    double          ll      = 0.0;
    double          z;
    double          a0, a1, b1, gamma, phi, c, sigma, rho;
    double         *yp, *nu, *mp, *hp;
    register double P, f;
    char           *names[4] = {"m", "h", "res", "ll"};

    PROT2(y = AS_NUMERIC(y), numprot);
    PROT2(par = AS_NUMERIC(par), numprot);

    T = length(y);

    PROT2(h = NEW_NUMERIC(T), numprot);
    PROT2(m = NEW_NUMERIC(T), numprot);
    PROT2(res = NEW_NUMERIC(T), numprot);
    PROT2(ans = NEW_LIST(4), numprot);

    ok = NPAR == length(par) && T > 0;

    if (!ok) error("Bad args passed to kal1");

    a0    = REAL(par)[0];
    a1    = REAL(par)[1];
    b1    = REAL(par)[2];
    gamma = REAL(par)[3];
    phi   = REAL(par)[4];
    c     = REAL(par)[5];
    sigma = REAL(par)[6];
    rho   = REAL(par)[7];

    yp = REAL(y);
    nu = REAL(res);
    mp = REAL(m);
    hp = REAL(h);

    ll = T * log(M_2PI);

    /* Do t = 1 stuff */

    hp[0] = a0 / (1.0 - b1 - a1 * (1.0 + SQR(gamma))); /* h(1)   */
    mp[0] = phi * (c / (1.0 - phi)) + c;               /* m(1|0) */
    nu[0] = yp[0] - mp[0];

    P = SQR(sigma); /* P(1|0) */

    /* Loop */
    for (t = 0; t < T - 1; t++)
    {
        nu[t] = yp[t] - mp[t];

        f = P + 2.0 * sigma * rho * sqrt(hp[t]) + hp[t];

        ll += log(f) + SQR(nu[t]) / f;

        mp[t + 1] =
            c + phi * (mp[t] + (P + sigma * rho * sqrt(hp[t])) / f * nu[t]);

        P = SQR(phi) * P - SQR(phi) * SQR(P + sigma * rho * sqrt(hp[t])) / f +
            SQR(sigma);

        z = nu[t] / hp[t];

        hp[t + 1] = a0 + a1 * hp[t] * SQR(z - gamma) + b1 * hp[t];
    }

    nu[T - 1] = yp[T - 1] - mp[T - 1];
    f         = P + 2.0 * sigma * rho * sqrt(hp[T - 1]) + hp[T - 1];

    ll += log(f) + SQR(nu[t]) / f;

    ll *= -0.5;

    SET_ELT(ans, 0, m);
    SET_ELT(ans, 1, h);
    SET_ELT(ans, 2, res);
    SET_ELT(ans, 3, ScalarReal(ll));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}
#undef NPAR
