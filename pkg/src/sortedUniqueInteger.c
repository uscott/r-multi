
#include "multi.h"

SEXP sortedUniqueInteger(SEXP t)
{
    int  numprot = 0;
    long tlen, i, j;
    SEXP tcopy, tunique;

    PROT2(t = AS_INTEGER(t), numprot);
    tlen = length(t);
    PROT2(tcopy = NEW_INTEGER(tlen), numprot);
    PROT2(tunique = NEW_INTEGER(tlen), numprot);

    memcpy(INTEGER(tcopy), INTEGER(t), tlen * sizeof(int));

    R_isort(INTEGER(tcopy), tlen);

    *INTEGER(tunique) = *INTEGER(tcopy);

    for (i = j = 1; i < tlen && INTEGER(tcopy)[i] != NA_INTEGER; i++)
        if (INTEGER(tcopy)[i] != INTEGER(tcopy)[i - 1])
            INTEGER(tunique)[j++] = INTEGER(tcopy)[i];

    PROT2(SET_LENGTH(tunique, j), numprot);

    UNPROT2;

    return tunique;
}
