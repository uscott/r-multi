#include "uri.h"

/*******************************************************************
 *
 *  Description: Computes the maximum drawdown where
 *    x[0], ..., x[n - 1] is a series of _cumulative_ P&Ls.
 *
 *******************************************************************/

void maxdd(double *x, long n, double *mdd, long *b, long *e)
{
    long    i, j;
    double *y, *z;

    *mdd = 0.0;

    *e = *b = 0;

    for (i = 0; i < n - 1; i++)
    {
        y = x + i;
        z = y + 1;

        for (j = i + 1; j < n; j++, z++)
        {
            if (*y - *z > *mdd)
            {
                *mdd = *y - *z;
                *b   = i;
                *e   = j;
            }
        }
    }
}

/*******************************************************************
 *
 *  Description: Computes the maximal drawdown of the cumulative
 *    or incremental P&L series x.
 *
 *******************************************************************/

SEXP maxDrawdown(SEXP x, SEXP incremental)
{
#define ANS_LEN 3
    SEXP      ans, y = R_NilValue;
    double    mdd;
    int       numprot = 0;
    long      b, e, n, i;
    const int incr    = asInteger(incremental);
    char     *names[] = {"maxdd", "beginning", "end"};

    n = length(x);

    if (incr)
    {
        y = x;

        PROT2(x = NEW_NUMERIC(n + 1), numprot);

        for (i = 1; i < n + 1; i++)
            REAL(x)[i] = REAL(y)[i - 1] + REAL(x)[i - 1];
    }

    n = length(x);

    maxdd(REAL(x), n, &mdd, &b, &e);

    PROT2(ans = NEW_NUMERIC(ANS_LEN), numprot);

    REAL(ans)[0] = mdd;
    REAL(ans)[1] = (double)b + 1;
    REAL(ans)[2] = (double)e + 1;

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
#undef ANS_LEN
}
