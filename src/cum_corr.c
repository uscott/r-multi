#include "uri.h"

SEXP cum_corr(SEXP x, SEXP y, SEXP reverse)
{
    int    numprot = 0;
    long   i, n;
    double dum;
    SEXP   corr, sx, sy;

    PROT2(corr = cum_cov(x, y, reverse), numprot);
    PROT2(sx = cum_var(x, reverse), numprot);
    PROT2(sy = cum_var(y, reverse), numprot);

    n = length(corr);
    if (length(sx) != n || length(sy) != n) error("huh?");

    for (i = 0; i < n; i++)
    {
        if (0 < (dum = REAL(sx)[i] * REAL(sy)[i]))
            REAL(corr)[i] /= sqrt(dum);
        else
            REAL(corr)[i] = NA_REAL;
    }

    UNPROT2;
    return corr;
}
