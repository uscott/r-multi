#include "uri.h"

SEXP cum_var(SEXP x, SEXP reverse)
{
    int       numprot = 0;
    const int rev     = asInteger(reverse);
    long      i, n;
    SEXP      rval = R_NilValue;
    double   *px = 0, *xr = 0;

    ENSURE_NUMERIC(x, numprot);

    n = length(x);
    if (!n)
    {
        UNPROT2;
        return R_NilValue;
    }

    PROT2(rval = NEW_NUMERIC(n), numprot);

    if (rev)
    {
        if (0 == (xr = (double *)malloc(n * DSIZE))) error("Malloc failed.");
        drev(n, REAL(x), xr);
    }

    px = rev ? xr : REAL(x);

    for (i = 0; i < n; i++) REAL(rval)[rev ? n - i - 1 : i] = dvar(px++, n - i);

    if (rev) free(xr);
    SET_NAMES(rval, GET_NAMES(x));
    UNPROT2;
    return rval;
}
