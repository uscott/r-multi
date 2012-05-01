
#include "uri.h"

#include <time.h>


#define THREE 3
#define BIG_POP_SIZE(n) (2 * (n) + (2 * (n) - 1) * (n))
#define C2(n) (((n) * ((n) - 1)) / 2)








/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains the code for the gradient-descent
**  optimizer to called inside other C functions.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



#define NUM_TRIES 2
static int numDeriv2(optimFun f, SEXP x, SEXP controlPar, 
                     SEXP gradient, SEXP eps, SEXP xNew)

  /*******************************************************************
   *
   *  Description: Numerically estimates gradient vector at point x
   *    of function x |-> f(x, controlPar).
   *
   *******************************************************************/

{

  int ok = 0;
  int tryNum = 0;
  long i;
  const long xLen = length(x);
  double val = 0.0, valNew = 0.0;
  double *xPtr = REAL(x), *xNewPtr = REAL(xNew), 
    *gradPtr = REAL(gradient), *epsPtr = REAL(eps);

  ok = IS_NUMERIC(eps) && IS_NUMERIC(gradient) && IS_NUMERIC(xNew) && IS_NUMERIC(x) &&
    (xLen == length(gradient)) && (xLen == length(xNew)) && (xLen == length(eps));

  if (!ok) error ("bad params passed to numDeriv2");

  val = f(x, controlPar);

  if (!R_FINITE(val)) return 0;
  
  dsetna(gradPtr, xLen); /* Initialize elements of gradient to NA */
  memcpy(xNewPtr, xPtr, xLen * sizeof(double));

  GetRNGstate();

  for (i = 0; i < xLen; i++, gradPtr++, xNewPtr++, epsPtr++, xPtr++) {

    if (0.0 == ABS(*epsPtr) || !R_FINITE(*epsPtr))
      *epsPtr = sqrt(DOUBLE_EPS);

    if (unif_rand() < 0.5)
      *epsPtr *= -1;

    do {

      *xNewPtr = *xPtr + *epsPtr;
      *epsPtr *= runif(0.9, 1);
      tryNum += 1;

      valNew = f(xNew, controlPar);

    } while (!R_FINITE(valNew) && tryNum < NUM_TRIES);

    if (!R_FINITE(valNew))
      return 0;

    *gradPtr = (valNew - val) / epsPtr[0];
    *xNewPtr = *xPtr;
    tryNum = 0;

    if (0.0 == *gradPtr) /* perturbation might be too small */
      *epsPtr *= runif(1, 2);

  }

  PutRNGstate();

  return 1;

}
#undef NUM_TRIES









#define OPTIM_ANS_LEN 3
#define OPTIM_ANS_NAMES {"par","value","convergence"}
SEXP optimGradient1(optimFun f, SEXP initPar, SEXP controlPar,
                    const double tol, const int relTol, 
                    const int minimize, const long maxit)

  /*******************************************************************
   *
   *  Description: Basic gradient-descent optimizer intended for
   *    calling from inside other C functions and not from R.
   *    The function being optimized is x |-> f(x, controlPar).
   *
   *  Status: Not finished.
   *
   *******************************************************************/

{
	
  const double sgn  = minimize ? 1.0 : -1.0;

  int numprot = 0;

  long i, it;
  int convergence = 1;
  long parLen = 0;
  double stepsize = 1e-1;
  double val1, val2, tolMult;
  register double delta;

  register double *xPtr, *xNewPtr, *gradPtr;
  char *ansNames[OPTIM_ANS_LEN] = OPTIM_ANS_NAMES;

  SEXP ans = R_NilValue, par1 = R_NilValue, par2 = R_NilValue, 
    gradient = R_NilValue, eps = R_NilValue, 
    *x = &par1, *xNew = &par2, *tmp = NULL;

  ENSURE_NUMERIC(initPar, numprot);

  parLen = length(initPar);

  PROT2(ans      = NEW_LIST(OPTIM_ANS_LEN),   numprot);
  PROT2(gradient = NEW_NUMERIC(parLen), numprot);
  PROT2(par1     = NEW_NUMERIC(parLen), numprot);
  PROT2(par2     = NEW_NUMERIC(parLen), numprot);
  PROT2(eps = NEW_NUMERIC(parLen), numprot);

  x = &par1;
  xNew = &par2;

  set_names(ans, ansNames);

  val1 = f(initPar, controlPar);

  if (!R_FINITE(val1)) {

    SET_ELT(ans, 0, initPar);
    SET_ELT(ans, 1, ScalarReal(val1));
    SET_ELT(ans, 2, ScalarInteger(convergence));

    UNPROTECT(numprot);

    return ans;

  }

  memcpy(REAL(par1), REAL(initPar), parLen * sizeof(double));
  memcpy(REAL(par2), REAL(initPar), parLen * sizeof(double));

  for (i = 0; i < parLen; i++)
    REAL(eps)[i] = sqrt(DOUBLE_EPS);

  GetRNGstate();

  for (it = 0, convergence = 1; it < maxit && convergence; it++) {

    if (!numDeriv2(f, *x, controlPar, gradient, eps, *xNew))
      continue;

    stepsize *= runif(1, 1.1);
    xPtr = REAL(*x);
    xNewPtr = REAL(*xNew);
    gradPtr = REAL(gradient);

    for (i = 0; i < parLen; i++) {

      delta = -sgn * (*gradPtr++) * stepsize;
      *xNewPtr++ = *xPtr++  + delta;

    }

    val2 = f(*xNew, controlPar);

    if (R_FINITE(val2) && sgn * val2 < sgn * val1 ) {

      tolMult = relTol ? tol + ABS(val1) : 1.0;

      convergence = (ABS(val2-val1) >= tol * tolMult);

      val1 = val2;

      tmp = x;
      x = xNew;
      xNew = tmp;

    }
    else
      stepsize /= runif(1, 3);


  }

  PutRNGstate();

  SET_ELT(ans, 0, *x);
  SET_ELT(ans, 1, ScalarReal(val1));
  SET_ELT(ans, 2, ScalarInteger(convergence));


  UNPROTECT(numprot);

  return ans;

}






SEXP optimGradient2(optimFun f, SEXP initPar, SEXP controlPar,
                    const double tol, const int relTol, 
                    const int minimize, const long maxit)
{
  
  int numprot = 0;
  long j, m, n;
  char *ansNames[OPTIM_ANS_LEN] = OPTIM_ANS_NAMES;
  SEXP par = R_NilValue, val = R_NilValue, conv = R_NilValue,
    ans = R_NilValue, tmp1 = R_NilValue, tmp2 = R_NilValue,
    tmp3 = R_NilValue;
  
  PROT2(ans = NEW_LIST(OPTIM_ANS_LEN), numprot);
  set_names(ans, ansNames);

  if (IS_MATRIX(initPar)) {

    m = nrows(initPar);
    n = ncols(initPar);
    
    PROT2(par = allocMatrix(REALSXP, m, n), numprot);
    PROT2(val = NEW_NUMERIC(n), numprot);
    PROT2(conv = NEW_INTEGER(n), numprot);
    PROT2(tmp1 = NEW_NUMERIC(m), numprot);

    for (j = 0; j < n; j++) {

      memcpy(REAL(tmp1), matcol1(initPar, j), m * sizeof(double));
      
      tmp2 = optimGradient1(f, tmp1, controlPar, tol, relTol, minimize, maxit);
      PROTECT(tmp2);
      PROTECT(tmp3 = getListElt(tmp2, "par"));

      if (IS_NUMERIC(tmp3) && length(tmp3) == m)
        memcpy(matcol1(par, j), REAL(tmp3), m * sizeof(double));
      else
        dsetna(matcol1(par, j), m);
      
      REAL(val)[j] = asReal(getListElt(tmp2, "value"));
      INTEGER(conv)[j] = asInteger(getListElt(tmp2,"convergence"));

      UNPROTECT(2);

    }

    SET_ELT(ans, 0, par);
    SET_ELT(ans, 1, val);
    SET_ELT(ans, 2, conv);

  } else if (isNewList(initPar)) {

    n = length(initPar);
    PROT2(par = NEW_LIST(n), numprot);
    PROT2(val = NEW_NUMERIC(n), numprot);
    PROT2(conv = NEW_INTEGER(n), numprot);
    
    for (j = 0; j < n; j++) {

      PROTECT(tmp1 = GET_ELT(initPar, j));
      
      tmp2 = optimGradient1(f, tmp1, controlPar, tol, relTol, minimize, maxit);
      PROTECT(tmp2);

      SET_ELT(par, j, getListElt(tmp2, "par"));
      REAL(val)[j] = asReal(getListElt(tmp2, "value"));
      INTEGER(conv)[j] = asInteger(getListElt(tmp2,"convergence"));

      UNPROTECT(2);

    }

    SET_ELT(ans, 0, par);
    SET_ELT(ans, 1, val);
    SET_ELT(ans, 2, conv);

  } else {

    ans = optimGradient1(f, initPar, controlPar, tol, relTol, minimize, maxit);
    PROT2(ans, numprot);

  }

  UNPROT2;

  return ans;

}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following functions are for genetic algo optimizer intended for
**  being called from inside C functions and not from R.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/






static void optimGa2_chk(SEXP *initPar, const int parLen, const int basePopSize, int *numprot)
{

  long i = 0;
  long initLen = 0;
  SEXP tmp;

  if (basePopSize < 1 || parLen < 1)
    error ("Must have positive population size and positive parameter length");

  GetRNGstate();


  if (isNull(*initPar) || !length(*initPar)) {

    PROT2(*initPar = allocMatrix(REALSXP, parLen, basePopSize), *numprot);

    for (i = 0; i < parLen * basePopSize; i++)
      REAL(*initPar)[i] = norm_rand();

  }
  else if (isMatrix(*initPar)) {

    if (parLen != nrows(*initPar))
      error ("Inconsistent parameter length info");

    PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), *numprot);

    memcpy(REAL(tmp), REAL(*initPar), parLen * basePopSize * sizeof(double));

    for (i = length(*initPar); i < parLen * basePopSize; i++)
      REAL(tmp)[i] = norm_rand();

    *initPar = tmp;

  }
  else {

    PROT2(*initPar = AS_NUMERIC(*initPar), *numprot);

    initLen = length(*initPar);

    if (initLen != parLen)
      error ("Inconsistent parameter length info");

    PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), *numprot);
    
    memcpy(REAL(tmp), REAL(*initPar), parLen * basePopSize * sizeof(double));

    for (i = initLen; i < parLen * basePopSize; i++)
      REAL(tmp)[i] = norm_rand();

    *initPar = tmp;
    
  }
  

  PutRNGstate();

}





static void mutate(SEXP pop, const int basePopSize)

  /*******************************************************************
   *
   *  Description: Mutates population of column vectors 0 throught
   *    basePopSize - 1 in par putting the results in column vectors
   *    basePopSize throught 2 * basePopSize - 1.
   *
   *******************************************************************/

{
  int ok, j;
  const int parLen = nrows(pop);
  register double *ptr1 = NULL, *ptr2 = NULL;

  ok = isMatrix(pop) && basePopSize > 0 && ncols(pop) >= 2*basePopSize;

  if (!ok)
    error ("bad argument passed to mutate");

  GetRNGstate();

  ptr1 = REAL(pop);
  ptr2 = REAL(pop) + parLen * basePopSize;


  for (j = 0; j < parLen * basePopSize; j++, ptr1++, ptr2++)
    *ptr2 = ISNAN(*ptr1) ? norm_rand() : rnorm(*ptr1, 1+ABS(*ptr1));

  PutRNGstate();

}




static void breed(SEXP pop, const int basePopSize)

  /*******************************************************************
   *
   *  Description: Breeds population of column vectors 0 through
   *
   *******************************************************************/

{
  const int parLen = nrows(pop);
  const int bigPopSize = ncols(pop);

  register int counter;

  register int i, j, k, ok;
  register double w;
  register double *ptr1, *ptr2, *kidPtr;

  ok = 
    isMatrix(pop) && parLen &&
    bigPopSize == BIG_POP_SIZE(basePopSize);

  if (!ok)
    error ("bad arguments in breed");

  GetRNGstate();

  kidPtr = REAL(pop) + 2 * parLen * basePopSize;

  counter = 0;

  for (i = 0; i < 2*basePopSize - 1; i++)
    for (j = i + 1; j < 2*basePopSize; j++) {

      ptr1 = matcol1(pop, i);
      ptr2 = matcol1(pop, j);

      for(k = 0; k < parLen; k++) {

        w = norm_rand() + 0.5;

        *kidPtr++ = w * (*ptr1++) + (1.0 - w) * (*ptr2++);

      }

      if (++counter >= bigPopSize)
        error ("breed ... something wrong here");

    }

  PutRNGstate();

}





static void getNextGen(optimFun f, SEXP controlPar, SEXP pop,
                       SEXP fargs, SEXP fvals, SEXP ix,
                       SEXP tmp, const int basePopSize, const int minimize)
{

  const int bigPopSize = ncols(pop);
  const int parLen = nrows(pop);
  const double sgn = minimize ? 1.0 : -1.0;

  register int ok, *ixPtr;
  long i, j;
  double val;

  ok = 
    isMatrix(pop) && isMatrix(tmp) && IS_INTEGER(ix) &&
    IS_NUMERIC(fargs) && IS_NUMERIC(fvals) && bigPopSize == length(ix) &&
    bigPopSize == length(fvals) && bigPopSize == BIG_POP_SIZE(basePopSize) && 
    parLen == length(fargs) &&
    parLen == nrows(tmp) &&
    basePopSize == ncols(tmp);

  if (!ok)
    error ("bad arguments to getNextGen");
  
  for(i = 0; i < bigPopSize; i++) 
    INTEGER(ix)[i] = i;

  mutate(pop, basePopSize);
  breed (pop, basePopSize);

  for (i = 0; i < bigPopSize; i++) {

    memcpy(REAL(fargs), matcol1(pop, i), parLen * sizeof(double));

    val = f(fargs, controlPar);

    REAL(fvals)[i] = R_FINITE(val) ? sgn * val : HUGE;

  }

  /* Move top basePopSize population members to front */
  rsort_with_index(REAL(fvals), INTEGER(ix), bigPopSize);

  ixPtr = INTEGER(ix);
		
  for (j = 0; j < basePopSize; j++)
    memcpy(matcol1(tmp, j), matcol1(pop, *ixPtr++), parLen * sizeof(double));

  memcpy(REAL(pop), REAL(tmp), parLen * basePopSize * sizeof(double));
  /* Done sorting population */


  /* Scale back 1st basePopSize function values */
  for (i = 0; i < basePopSize; i++)
    REAL(fvals)[i] *= sgn;
  
}






SEXP optimGa2(optimFun f, SEXP initPar, SEXP controlPar,
              const int parLen, const int basePopSize, 
              const long stopLags, const long minit, const long maxit, 
              const int minimize, const double tol, const int relTol)

  /*******************************************************************
   *
   *  Description: Uses genetic algorithms with gradient descent
   *    at the end to optimize function x |-> f(x, controlPar).
   *
   *******************************************************************/

{

  int numprot = 0;  
  const int bigPopSize = BIG_POP_SIZE(basePopSize);
  int gen = 1, notConvergedYet = 1, numAnsCols = 1;
  register double fmin, fmax;
  register long lagNum;
  char *names[] = OPTIM_ANS_NAMES;
  int *dimPtr = NULL;

  SEXP fargs, ans, pop, fvals, tmp, ix, popDim;

  /* Check arguments and perform necessary adjustments */
  optimGa2_chk(&initPar, parLen, basePopSize, &numprot);


  PROT2(fvals = NEW_NUMERIC(bigPopSize), numprot);
  PROT2(fargs = NEW_NUMERIC(parLen),     numprot);
  PROT2(ix    = NEW_INTEGER(bigPopSize), numprot);
  PROT2(ans   = NEW_LIST(OPTIM_ANS_LEN), numprot);
  PROT2(pop = allocMatrix(REALSXP, parLen, bigPopSize),  numprot);
  PROT2(tmp = allocMatrix(REALSXP, parLen, basePopSize), numprot);

  GetRNGstate();

  memcpy(REAL(pop), REAL(initPar), parLen * basePopSize * sizeof(double));

  for (gen = 0; gen < maxit && notConvergedYet; gen++) {

    getNextGen(f, controlPar, pop, fargs, fvals, ix, tmp, basePopSize, minimize);

    fmin = REAL(fvals)[0];
    fmax = REAL(fvals)[basePopSize - 1];

    notConvergedYet = gen < minit || (ABS(fmax - fmin) >= tol);

    lagNum = notConvergedYet ? 0 : lagNum + 1;

    if (!notConvergedYet && lagNum < stopLags && gen < maxit-1)
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

  int numprot = 0;
  long i;
  double ans = 0;

  ENSURE_NUMERIC(x, numprot);

  for (i = 0; i < length(x); i++)
    ans += SQR(REAL(x)[i]);

  UNPROT2;

  return ans;

}




SEXP testOptim(SEXP initPar, SEXP controlPar, SEXP parLen, SEXP basePopSize, 
               SEXP stopLags, SEXP minit, SEXP maxit,
               SEXP minimize, SEXP tol, SEXP relTol, SEXP optimType)
{
  int type = asInteger(optimType);
  SEXP ans;

  switch (type) {

  case 0:

    ans = optimGa2(testOptimFun, initPar, controlPar,
                   asInteger(parLen), asInteger(basePopSize),
                   asInteger(stopLags), asInteger(minit), asInteger(maxit), 
                   asInteger(minimize), asReal(tol), asInteger(relTol));
    break;
    
  case 1:

    ans = optimGradient2(testOptimFun, initPar, controlPar,
                         asReal(tol), asInteger(relTol),
                         asInteger(minimize), asInteger(maxit));

    break;


  default:

    ans = R_NilValue;

    break;

  }

  return ans;

}
#undef DEBUG_THIS_SHIT
#endif
