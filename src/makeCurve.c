
#include "uri.h"

static void makeCurve_chk(SEXP *xlist, SEXP *t, SEXP *mthNum, SEXPTYPE **types,
                          int *numprot)
{
    long i, n, listlen;
    SEXP tmp;

    if (!xlist || !t || !mthNum || !types || !numprot) error("null pointer!");

    PROT2(*xlist = tolist(*xlist), *numprot);

    listlen = length(*xlist);
    *types  = (SEXPTYPE *)R_alloc(listlen, sizeof(SEXPTYPE));

    for (i = 0; i < listlen; i++)
    {
        tmp         = GET_ELT(*xlist, i);
        (*types)[i] = TYPEOF(tmp);

        SET_ELT(*xlist, i, AS_CHAR(tmp));
    }

    ENSURE_INTEGER(*t, *numprot);
    ENSURE_INTEGER(*mthNum, *numprot);

    n = length(GET_ELT(*xlist, 0));

    for (i = 1; i < length(*xlist); i++)
        if (n != length(GET_ELT(*xlist, i))) error("mismatched lengths");

    if (n != length(*t) || n != length(*mthNum)) error("mismatched lengths");
}

SEXP makeCurve(SEXP xlist, SEXP t, SEXP mthNum)
{
#define MAX_CRV_MTHS 36
    int  numprot = 0;
    long i, j, k, h, tlen = 0, listlen = 0, t0;
    long nrows, ncols, coli;
    int  index;
    SEXPTYPE
    *types;
    SEXP tUnique, xCrv, names, *xtmp;

    makeCurve_chk(&xlist, &t, &mthNum, &types, &numprot);

    tlen    = length(t);
    listlen = length(xlist);

    PROT2(tUnique = sortedUniqueInteger(t), numprot);

    if (!(nrows = length(tUnique))) return R_NilValue;
    if (!(ncols = imax(INTEGER(mthNum), tlen))) return R_NilValue;
    if (MIN(INTEGER(mthNum), tlen) <= 0)
        error("month numbers must be positive");

    ncols = MIN(ncols, MAX_CRV_MTHS);

    PROT2(xCrv = NEW_LIST(listlen), numprot);

    for (i = 0; i < listlen; ++i)
        SET_ELT(xCrv, i, allocMatrix(STRSXP, nrows, ncols));

    SET_NAMES(xCrv, GET_NAMES(xlist));

    xtmp = (SEXP *)R_alloc(listlen, sizeof(SEXP));

    for (i = 0; i < listlen; ++i) PROT2(xtmp[i] = NEW_STRING(ncols), numprot);

    for (i = 0; i < nrows; ++i)
    {
        t0 = INTEGER(tUnique)[i];

        for (h = 0; h < listlen; ++h)
            for (j = 0; j < ncols; ++j) SET_ELT(xtmp[h], j, mkChar("NA"));

        for (k = 0; k < tlen; k++)
            if (INTEGER(t)[k] == t0 && (coli = INTEGER(mthNum)[k] - 1) < ncols)
                for (h = 0; h < listlen; h++)
                    SET_ELT(xtmp[h], coli,
                            duplicate(STRING_ELT(GET_ELT(xlist, h), k)));

        for (h = 0; h < listlen; h++)
            for (j = 0; j < ncols; j++)
                SET_ELT(GET_ELT(xCrv, h), i + j * nrows,
                        STRING_ELT(xtmp[h], j));
    }

    for (i = 0; i < listlen; i++)
        SET_ELT(xCrv, i, coerceVector(GET_ELT(xCrv, i), types[i]));

    PROT2(names = GET_NAMES(xCrv), numprot);
    PROT2(SET_LENGTH(names, listlen + 1), numprot);
    PROT2(SET_LENGTH(xCrv, listlen + 1), numprot);

    SET_ELT(xCrv, listlen, tUnique);
    SET_ELT(names, listlen, mkChar("t"));

    SET_NAMES(xCrv, names);

    UNPROT2;

    return xCrv;
#undef MAX_CRV_MTHS
}

#ifdef NOT_YET_oaisejfd

#define MAX_CRV_MTHS 36
SEXP makeCurve2(SEXP xlist, SEXP dates, SEXP expDates)
{
    int       numprot = 0;
    long      i, j, k, h, tlen = 0, listlen = 0, t0;
    long      nrows, ncols, coli;
    int       index, *mthNumTmp, **mthNumRowPtr;
    SEXPTYPE *types;
    SEXP      tUnique, xCrv, names, *xtmp, t;
    uritm    *t, *expt;

    makeCurve_chk(&xlist, &dates, &expDates, &types, &numprot);

    tlen    = length(dates);
    listlen = length(xlist);

    t    = (uritm *)R_alloc(tlen, sizeof(uritm));
    expt = (uritm *)R_alloc(tlen, sizeof(uritm));

    if (!(nrows = length(tUnique))) return R_NilValue;
    if (!(ncols = imax(INTEGER(mthNum), tlen))) return R_NilValue;
    if (imin(INTEGER(mthNum), tlen) <= 0)
        error("month numbers must be positive");

    ncols = LMIN(ncols, MAX_CRV_MTHS);

    PROT2(xCrv = NEW_LIST(listlen), numprot);
    for (i = 0; i < listlen; i++)
        SET_ELT(xCrv, i, allocMatrix(STRSXP, nrows, ncols));

    SET_NAMES(xCrv, GET_NAMES(xlist));

    xtmp = (SEXP *)R_alloc(listlen, sizeof(SEXP));
    for (i = 0; i < listlen; i++) PROT2(xtmp[i] = NEW_STRING(ncols), numprot);

    for (i = 0; i < nrows; i++)
    {
        t0 = INTEGER(tUnique)[i];

        for (h = 0; h < listlen; h++)
            for (j = 0; j < ncols; j++) SET_ELT(xtmp[h], j, mkChar("NA"));

        for (k = 0; k < tlen; k++)
            if (INTEGER(t)[k] == t0 && (coli = INTEGER(mthNum)[k] - 1) < ncols)
                for (h = 0; h < listlen; h++)
                    SET_ELT(xtmp[h], coli,
                            duplicate(STRING_ELT(GET_ELT(xlist, h), k)));

        for (h = 0; h < listlen; h++)
            for (j = 0; j < ncols; j++)
                SET_ELT(GET_ELT(xCrv, h), i + j * nrows,
                        STRING_ELT(xtmp[h], j));
    }

    for (i = 0; i < listlen; i++)
        SET_ELT(xCrv, i, coerceVector(GET_ELT(xCrv, i), types[i]));

    PROT2(names = GET_NAMES(xCrv), numprot);
    PROT2(SET_LENGTH(names, listlen + 1), numprot);
    PROT2(SET_LENGTH(xCrv, listlen + 1), numprot);

    SET_ELT(xCrv, listlen, tUnique);
    SET_ELT(names, listlen, mkChar("t"));

    SET_NAMES(xCrv, names);

    UNPROT2;

    return xCrv;
}
#undef MAX_CRV_MTHS

#endif
