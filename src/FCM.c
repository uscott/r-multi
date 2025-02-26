#include "multi.h"

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This file contains implementations of the Fast Convolution Method
**  from Dan & Alex's paper.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*
**
**  This first implementation is a very rudimentary one
**  for European call options.  It is not intended to be
**  numerically efficient.  For simplicity it is
**  assumed that r = q = 0.  It only returns the price.
**
*/

SEXP fcm_EuroCall0(SEXP tau, SEXP price, SEXP strike, SEXP vol, SEXP nGrdPts,
                   SEXP nTimeSteps)
{
    const double T       = asReal(tau) / 365.0;
    const double S       = asReal(price);
    const double K       = asReal(strike);
    const double sig     = asReal(vol);
    const long   M       = asInteger(nGrdPts);
    const long   N       = asInteger(nTimeSteps);
    const double mu      = -0.5 * sig * sig;
    const double discFac = 1.0;
    const double dt      = T / N;
    const double h       = 5 * sig * sqrt(T) / M;
    const double w       = sig * sqrt(dt) / h;

    long            i, j, n;
    register double x1, x2, y1, y2, x;
    double         *z, *V, *a, *b;

    z = (double *)R_alloc(2 * M + 1, sizeof(double));
    V = (double *)R_alloc(2 * M + 1, sizeof(double));
    a = (double *)R_alloc(2 * M, sizeof(double));
    b = (double *)R_alloc(2 * M, sizeof(double));

    /* Initiate z and V */
    for (i = 0; i < 2 * M + 1; i++)
    {
        z[i] = mu * T + (-M + i) * h + log(S);
        V[i] = DMAX(exp(z[i]) - K, 0);
    }

    /* Work backwords through time */
    for (n = 0; n < N; n++)
    {
        /* Get linear interpolation coefficients */
        for (i = 0; i < 2 * M; i++)
        {
            y1 = V[i];
            y2 = V[i + 1];
            x1 = exp(z[i]);
            x2 = exp(z[i + 1]);

            a[i] = (y2 - y1) / (x2 - x1);
            b[i] = (x2 * y1 - x1 * y2) / (x2 - x1);
        }

        for (j = 0; j < 2 * M + 1; j++)
        {
            V[j] = 0;
            x    = exp(z[j] + sig * sig * dt / 2);

            for (i = 0; i < 2 * M; i++)
            {
                V[j] += a[i] * x *
                        (NORM_DIST((i - j + 1) / w - sig * sqrt(dt)) -
                         NORM_DIST((i - j) / w - sig * sqrt(dt)));
                V[j] += b[i] *
                        (NORM_DIST((i - j + 1) / w) - NORM_DIST((i - j) / w));
            }

            V[j] *= discFac;
        }

        /* Update z */
        for (i = 0; i < 2 * M + 1; i++) z[i] -= mu * dt;
    }

    return ScalarReal(V[M]);
}
