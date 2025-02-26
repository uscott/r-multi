/*
** Seems to no longer be in use.
*/
#include <time.h>

#include "multi.h"

#define THREE           3
#define BIG_POP_SIZE(n) (2 * (n) + (2 * (n) - 1) * (n))
#define C2(n)           (((n) * ((n) - 1)) / 2)

#define OPTIM_ANS_LEN   3
#define OPTIM_ANS_NAMES {"par", "value", "convergence"}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section is for a GA optimizer in which integer-only
**  parameters are accomodated.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

static void optimGa3_chk(SEXP *initPar, SEXP *parType, const int parLen,
                         const int basePopSize, int *numprot)
{
    register long i       = 0;
    long          initLen = 0;
    register char type;
    SEXP          tmp;

    if (basePopSize < 1 || parLen < 1)
        error(
            "Must have positive population size and positive parameter length");

    GetRNGstate();

    ENSURE_CHAR(*parType, *numprot);

    if (1 == length(*parType))
        PROT2(*parType = charrep(*parType, parLen), *numprot);

    if (parLen != length(*parType))
        error("parType must have length %d", parLen);

    if (isNull(*initPar) || !length(*initPar))
    {
        PROT2(*initPar = allocMatrix(REALSXP, parLen, basePopSize), *numprot);

        for (i = 0; i < parLen * basePopSize; i++)
        {
            type = *CHAR(GET_ELT(*parType, i % parLen));

            if ('i' == type || 'I' == type)
                REAL(*initPar)[i] = fprec(rnorm(0, 10), 0);
            else
                REAL(*initPar)[i] = norm_rand();
        }
    } else if (isMatrix(*initPar))
    {
        if (parLen != nrows(*initPar))
            error("Inconsistent parameter length info");

        PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), *numprot);

        memcpy(REAL(tmp), REAL(*initPar),
               parLen * basePopSize * sizeof(double));

        for (i = 0; i < length(*initPar); i++)
        {
            type = *CHAR(GET_ELT(*parType, i % parLen));

            if ('i' == type || 'I' == type)
                REAL(tmp)[i] = fprec(REAL(tmp)[i], 0);
        }

        for (i = length(*initPar); i < parLen * basePopSize; i++)
        {
            type = *CHAR(GET_ELT(*parType, i % parLen));

            if ('i' == type || 'I' == type)
                REAL(tmp)[i] = fprec(rnorm(0, 10), 0);
            else
                REAL(tmp)[i] = norm_rand();
        }

        *initPar = tmp;

    } else
    {
        PROT2(*initPar = AS_NUMERIC(*initPar), *numprot);

        initLen = length(*initPar);

        if (initLen != parLen) error("Inconsistent parameter length info");

        PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), *numprot);

        memcpy(REAL(tmp), REAL(*initPar),
               parLen * basePopSize * sizeof(double));

        for (i = 0; i < initLen; i++)
        {
            type = *CHAR(GET_ELT(*parType, i % parLen));

            if ('i' == type || 'I' == type)
                REAL(tmp)[i] = fprec(REAL(tmp)[i], 0);
        }

        for (i = initLen; i < parLen * basePopSize; i++)
        {
            type = *CHAR(GET_ELT(*parType, i % parLen));

            if ('i' == type || 'I' == type)
                REAL(tmp)[i] = fprec(rnorm(0, 10), 0);
            else
                REAL(tmp)[i] = norm_rand();
        }

        *initPar = tmp;
    }

    PutRNGstate();
}

static void mutate3(SEXP pop, SEXP parType, const int basePopSize)

/*******************************************************************
 *
 *  Description: Mutates population of column vectors 0 throught
 *    basePopSize - 1 in par putting the results in column vectors
 *    basePopSize throught 2 * basePopSize - 1.
 *
 *******************************************************************/

{
    int              ok, j;
    const int        parLen = nrows(pop);
    register double *ptr1 = NULL, *ptr2 = NULL;
    register char    type;

    ok = isMatrix(pop) && basePopSize > 0 && ncols(pop) >= 2 * basePopSize &&
         length(parType) == parLen && IS_CHARACTER(parType);

    if (!ok) error("bad argument passed to mutate");

    GetRNGstate();

    ptr1 = REAL(pop);
    ptr2 = REAL(pop) + parLen * basePopSize;

    for (j = 0; j < parLen * basePopSize; j++, ptr1++, ptr2++)
    {
        type = *CHAR(GET_ELT(parType, j % basePopSize));

        *ptr2 = ISNAN(*ptr1) ? norm_rand() : rnorm(*ptr1, 1 + ABS(*ptr1));

        if ('i' == type || 'I' == type) *ptr2 = fprec(*ptr2, 0);
    }

    PutRNGstate();
}

static void breed3(SEXP pop, SEXP parType, const int basePopSize)

/*******************************************************************
 *
 *  Description: Breeds population of column vectors 0 through
 *
 *******************************************************************/

{
    const int        parLen     = nrows(pop);
    const int        bigPopSize = ncols(pop);
    register int     counter;
    register int     i, j, k, ok;
    register double  w;
    register double *ptr1, *ptr2, *kidPtr;
    register char    type;

    ok = isMatrix(pop) && parLen && bigPopSize == BIG_POP_SIZE(basePopSize) &&
         IS_CHARACTER(parType) && length(parType) == parLen;

    if (!ok) error("bad arguments in breed");

    GetRNGstate();

    kidPtr = REAL(pop) + 2 * parLen * basePopSize;

    counter = 0;

    for (i = 0; i < 2 * basePopSize - 1; i++)
        for (j = i + 1; j < 2 * basePopSize; j++)
        {
            ptr1 = matcol1(pop, i);
            ptr2 = matcol1(pop, j);

            for (k = 0; k < parLen; k++, kidPtr++)
            {
                type = *CHAR(GET_ELT(parType, k));
                w    = rnorm(0.5, 1);

                *kidPtr = w * (*ptr1++) + (1.0 - w) * (*ptr2++);

                if ('i' == type || 'I' == type) *kidPtr = fprec(*kidPtr, 0);
            }

            if (++counter >= bigPopSize)
                error("breed ... something wrong here");
        }

    PutRNGstate();
}

static void getNextGen3(optimFun f, SEXP controlPar, SEXP pop, SEXP parType,
                        SEXP fargs, SEXP fvals, SEXP ix, SEXP tmp,
                        const int basePopSize, const int minimize)
{
    const int    bigPopSize = ncols(pop);
    const int    parLen     = nrows(pop);
    const double sgn        = minimize ? 1.0 : -1.0;

    register int ok, *ixPtr;
    long         i, j;
    double       val;

    ok = isMatrix(pop) && isMatrix(tmp) && IS_INTEGER(ix) &&
         IS_NUMERIC(fargs) && IS_NUMERIC(fvals) && bigPopSize == length(ix) &&
         bigPopSize == length(fvals) &&
         bigPopSize == BIG_POP_SIZE(basePopSize) && parLen == length(fargs) &&
         parLen == nrows(tmp) && basePopSize == ncols(tmp);

    if (!ok) error("bad arguments to getNextGen");

    for (i = 0; i < bigPopSize; i++) INTEGER(ix)[i] = i;

    mutate3(pop, parType, basePopSize);
    breed3(pop, parType, basePopSize);

    for (i = 0; i < bigPopSize; i++)
    {
        memcpy(REAL(fargs), matcol1(pop, i), parLen * sizeof(double));

        val = f(fargs, controlPar);

        REAL(fvals)[i] = R_FINITE(val) ? sgn * val : HUGE;
    }

    /* Move top basePopSize population members to front */
    rsort_with_index(REAL(fvals), INTEGER(ix), bigPopSize);

    ixPtr = INTEGER(ix);

    for (j = 0; j < basePopSize; j++)
        memcpy(matcol1(tmp, j), matcol1(pop, *ixPtr++),
               parLen * sizeof(double));

    memcpy(REAL(pop), REAL(tmp), parLen * basePopSize * sizeof(double));
    /* Done sorting population */

    /* Scale back 1st basePopSize function values */
    for (i = 0; i < basePopSize; i++) REAL(fvals)[i] *= sgn;
}

SEXP optimGa3(optimFun f, SEXP initPar, SEXP controlPar, SEXP parType,
              const int parLen, const int basePopSize, const long stopLag,
              const long minit, const long maxit, const int minimize,
              const double tol, const int relTol)

/*******************************************************************
 *
 *  Description: Uses genetic algorithms to optimize function x |-> f(x,
 *controlPar).
 *
 *******************************************************************/

{
    int             numprot    = 0;
    const int       bigPopSize = BIG_POP_SIZE(basePopSize);
    int             gen = 1, notConvergedYet = 1, numAnsCols = 1;
    register double fmin, fmax;
    register long   lagNum  = 0;
    char           *names[] = OPTIM_ANS_NAMES;
    int            *dimPtr  = NULL;

    SEXP fargs, ans, pop, fvals, tmp, ix, popDim;

    /* Check arguments and perform necessary adjustments */
    optimGa3_chk(&initPar, &parType, parLen, basePopSize, &numprot);

    PROT2(fvals = NEW_NUMERIC(bigPopSize), numprot);
    PROT2(fargs = NEW_NUMERIC(parLen), numprot);
    PROT2(ix = NEW_INTEGER(bigPopSize), numprot);
    PROT2(ans = NEW_LIST(OPTIM_ANS_LEN), numprot);

    PROT2(pop = allocMatrix(REALSXP, parLen, bigPopSize), numprot);
    PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), numprot);

    GetRNGstate();

    memcpy(REAL(pop), REAL(initPar), parLen * basePopSize * sizeof(double));

    for (gen = 0; gen < maxit && notConvergedYet; gen++)
    {
        getNextGen3(f, controlPar, pop, parType, fargs, fvals, ix, tmp,
                    basePopSize, minimize);

        fmin = REAL(fvals)[0];
        fmax = REAL(fvals)[basePopSize - 1];

        notConvergedYet = gen < minit || (ABS(fmax - fmin) >= tol);

        lagNum = notConvergedYet ? 0 : lagNum + 1;

        if (!notConvergedYet && lagNum < stopLag && gen < maxit - 1)
            notConvergedYet = 1;
    }

    PutRNGstate();

    numAnsCols = notConvergedYet ? basePopSize : 1;

    PROT2(popDim = NEW_INTEGER(2), numprot);
    dimPtr = INTEGER(popDim);

    dimPtr[0] = parLen;
    dimPtr[1] = numAnsCols;

    PROT2(SET_LENGTH(pop, dimPtr[0] * dimPtr[1]), numprot);
    setAttrib(pop, R_DimSymbol, popDim);

    SET_ELT(ans, 0, pop);
    SET_ELT(ans, 1, SET_LENGTH(fvals, numAnsCols));
    SET_ELT(ans, 2, ScalarInteger(notConvergedYet));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}

#undef OPTIM_ANS_LEN
#undef OPTIM_ANS_NAMES

#undef BIG_POP_SIZE
#undef C2

#define DEBUG_THIS_SHIT
#ifdef DEBUG_THIS_SHIT

static double testOptimFun(SEXP x, SEXP controlPar)
{
    int    numprot = 0;
    long   i;
    double ans = 0;

    ENSURE_NUMERIC(x, numprot);

    for (i = 0; i < length(x); i++) ans += SQR(REAL(x)[i]);

    UNPROT2;

    return ans;
}

SEXP testOptim3(SEXP initPar, SEXP parType, SEXP controlPar, SEXP parLen,
                SEXP basePopSize, SEXP stopLag, SEXP minit, SEXP maxit,
                SEXP minimize, SEXP tol, SEXP relTol, SEXP optimType)
{
    int  type = asInteger(optimType);
    SEXP ans;

    switch (type)
    {
        case 0:

            ans =
                optimGa3(testOptimFun, initPar, controlPar, parType,
                         asInteger(parLen), asInteger(basePopSize),
                         asInteger(stopLag), asInteger(minit), asInteger(maxit),
                         asInteger(minimize), asReal(tol), asInteger(relTol));
            break;

        default:

            ans = R_NilValue;

            break;
    }

    return ans;
}
#undef DEBUG_THIS_SHIT
#endif
