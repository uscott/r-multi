
#include "uri.h"

/* Sorts [x] and applies same permutation to [y] */
void sortWithIndex2(double *x, double *y, long n)
{
    long    i;
    int    *ind   = (int *)R_alloc(n, sizeof(int));
    double *ycopy = (double *)R_alloc(n, sizeof(double));

    for (i = 0; i < n; i++)
    {
        ind[i]   = i;
        ycopy[i] = y[i];
    }

    rsort_with_index(x, ind, n);

    for (i = 0; i < n; i++) y[i] = ycopy[ind[i]];
}
