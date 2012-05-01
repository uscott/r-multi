
#include "uri.h"

#include <time.h>


#define THREE 3
#define BIG_POP_SIZE(n) (2 * (n) + (2 * (n) - 1) * (n))
#define C2(n) (((n) * ((n) - 1)) / 2)






/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains the basic gradient-descent optimizer
**  callable from R using .Call.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/






SEXP optim_gradient1(const SEXP fn,    
                     const SEXP init,
                     const SEXP tol, 
                     const SEXP minimize,
                     const SEXP maxit,
                     const SEXP env)

  /*******************************************************************
   *
   *  Description: Basic gradient-descent optimizer.
   *
   *******************************************************************/

{
  const int n       = length(init);
  const int mit     = MAX(asInteger(maxit), 1);
  const double dtol = asReal(tol);
  const double sgn  = asInteger(minimize) ? 1.0 : -1.0;


  int i, it, ok;
  int convergence;
  int numprot = 0;
  double stepsize = 1e-1;
  double gradlen = R_PosInf;
  register double *pp1, *pp2, *vp1, *vp2, *gp;

  SEXP val1, val2, ans, R_fcall, par1, par2, gradient, ans_names;

  ok = 
    isFunction(fn) &&
    isEnvironment(env) &&
    isInteger(maxit) &&
    isNumeric(init);

  if(!ok)
    error ("error in optim_gradient1");

  PROTECT(ans       = NEW_LIST(4));
  PROTECT(val1      = NEW_NUMERIC(1));
  PROTECT(val2      = NEW_NUMERIC(1));
  PROTECT(par1      = NEW_NUMERIC(n));
  PROTECT(par2      = NEW_NUMERIC(n));
  PROTECT(gradient  = NEW_NUMERIC(n));
  PROTECT(ans_names = NEW_STRING(length(ans)));
  PROTECT(R_fcall   = lang2(fn, R_NilValue));

  numprot += 8;

  vp1 = REAL(val1);
  vp2 = REAL(val2);
  
  if (R_finite(*REAL(feval(fn, init, env)))) {

    memcpy(REAL(par1), REAL(init), n * sizeof(double));
    memcpy(REAL(par2), REAL(init), n * sizeof(double));

    GetRNGstate();

    for (it = 0, convergence = 1; it < mit && convergence; it++) {

      pp1 = REAL(par1);
      pp2 = REAL(par2);

      stepsize *= 1.0 + 0.1 * unif_rand();

      PROTECT(gradient = num_deriv(fn, par1, env));

      gradlen = *REAL(vector_norm(gradient));
      gp      = REAL(gradient);

      dnormalize(gp, n); /* Normalize gradient vector */

      *vp1 = *REAL(feval(fn, par1, env));

      for (i = 0; i < n; i++)
	*pp2++ = *pp1++ - sgn * (*gp++) * MAX(1.0, fabs(*vp1)) * stepsize;

      *vp2 = *REAL(feval(fn, par2, env));

      pp1 = REAL(par1);
      pp2 = REAL(par2);

      if (R_FINITE(*vp2) && sgn * (*vp2) < sgn * (*vp1)) {

        convergence = (ABS(*vp2 - *vp1) >= dtol * (dtol + ABS(*vp1)));

        *vp1 = *vp2;
	memcpy(pp1, pp2, n * sizeof(double));
      
      }
      else
	stepsize /= (1.0 + 2.0 * unif_rand());


      UNPROTECT(1);
    }

    PutRNGstate();

    SET_ELT(ans, 0, par1);
    SET_ELT(ans, 1, ScalarReal(*vp1));
    SET_ELT(ans, 2, ScalarInteger(convergence));
    SET_ELT(ans, 3, ScalarReal(gradlen));

  } else {

    for (i = 0; i < length(ans); i++)
      SET_ELT(ans, i, ScalarReal(NA_REAL));
  
  }

  CHAR_PTR(ans_names)[0] = mkChar("par");
  CHAR_PTR(ans_names)[1] = mkChar("value");
  CHAR_PTR(ans_names)[2] = mkChar("convergence");
  CHAR_PTR(ans_names)[3] = mkChar("gradient.length");

  SET_NAMES(ans, ans_names);

  UNPROTECT(numprot);
	
  return ans;
}









/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  The following section contains functions for optimization using
**  genetic algorithms.  The main optimization function optim_ga1
**  is intended to be called from R.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



static void breed1(SEXP par, SEXP ans)

  /*******************************************************************
   *
   *  Description: Breeds population of column vectors given by
   *    par and puts answer in ans.
   *
   *******************************************************************/

{
  const int m = nrows(par);
  const int n = ncols(par);

  register int i, j, k, ok;
  register double w;
  register double *p1, *p2, *kp;

  ok = 
    isMatrix(par) && 
    isMatrix(ans) && 
    0 < nrows(par) &&
    nrows(par) == nrows(ans) && 
    ncols(ans) == ncols(par) + C2(ncols(par));

  if (!ok)
    error ("bad arguments in breed1");

  GetRNGstate();

  memcpy(REAL(ans), REAL(par), m * n * sizeof(double));

  kp = REAL(ans) + m * n;

  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {

      p1 = matcol1(par, i);
      p2 = matcol1(par, j);

      for(k = 0; k < m; k++) {
        w = norm_rand() + 0.5;
        *kp++ = w * (*p1++) + (1 - w) * (*p2++);
      }
    }

  PutRNGstate();

}




static void mutate1(SEXP par, SEXP mut)

  /*******************************************************************
   *
   *  Description: Mutates population of column vectors in par
   *    and puts answer in ans.
   *
   *******************************************************************/

{
  int ok, j;
  const int         m = nrows(par);
  const int         n = ncols(par);
  register double *pp = NULL;
  register double *mp = NULL;

  ok = 
    isMatrix(par) &&
    isMatrix(mut) &&
    ncols(par) > 0 &&
    ncols(mut) == 2 * ncols(par) &&
    MAX(nrows(par), nrows(mut)) > 0 &&
    nrows(par) == nrows(mut);


  if (!ok)
    error ("bad argument passed to mutate1");

  GetRNGstate();

  memcpy(REAL(mut), REAL(par), m * n * sizeof(double));

  pp = REAL(par);
  mp = REAL(mut) + m * n;

  for (j = 0; j < m * n; j++, pp++, mp++)
    *mp = ISNA(*pp) ? norm_rand() : *pp + norm_rand() * ABS(*pp);

  PutRNGstate();

}






SEXP optim_ga1(SEXP fn, SEXP init, SEXP maxit, 
               SEXP minimize, SEXP tol, SEXP env)
{
  const int m          = nrows(init);
  const int n          = ncols(init);	
  const int mg         = asInteger(maxit);
  const int len[THREE] = {n, 2 * n, BIG_POP_SIZE(n)};
  const double sgn     = asInteger(minimize) ? 1.0 : -1.0;
  const double dtol    = asReal(tol);
  int i, j, gen, ok, convergence;
  int numprot = 0;	
  register int *ix;
  register double f, fmin, fmax;
  char *names[] = {"par", "value", "convergence"};
  SEXP fvals, fargs, ans, pop[THREE];

  ok = isFunction(fn) && isEnvironment(env) & isMatrix(init);

  if (!ok)
    error ("optim_ga1: bad args");

  ix = (int *) malloc(len[2] * sizeof(int));

  PROT2(fvals = NEW_NUMERIC(len[2]), numprot);
  PROT2(fargs = NEW_NUMERIC(m),      numprot);
  PROT2(ans   = NEW_LIST(THREE),     numprot);

  for (i = 0; i < THREE; i++)
    PROT2(pop[i] = allocMatrix(REALSXP, m, len[i]), numprot);


  GetRNGstate();

  memcpy(REAL(pop[0]), REAL(init), m * n * sizeof(double));


  for (gen = 0, convergence = 1; gen < mg && convergence; gen++) {

    for(i = 0; i < len[2]; i++) 
      ix[i] = i;


    mutate1(pop[0], pop[1]);
    breed1 (pop[1], pop[2]);


    for (i = 0; i < len[2]; i++) {

      memcpy(REAL(fargs), matcol1(pop[2], i), m * sizeof(double));

      f = *REAL(feval(fn, fargs, env));

      REAL(fvals)[i] = R_FINITE(f) ? sgn * f : HUGE;

    }


    rsort_with_index(REAL(fvals), ix, len[2]);
		
    for (j = 0; j < len[0]; j++)
      memcpy(matcol1(pop[0], j), matcol1(pop[2], ix[j]), m * sizeof(double));

    fmin = REAL(fvals)[0];
    fmax = REAL(fvals)[len[0] - 1];

    convergence = (fmax - fmin >= dtol);
     
     
  }

  PutRNGstate();

  for (i = 0; i < n; i++)
    REAL(fvals)[i] *= sgn;

  SET_ELT(ans, 0, pop[0]);
  SET_ELT(ans, 1, SET_LENGTH(fvals, len[0]));  
  SET_ELT(ans, 2, ScalarInteger(convergence));

  set_names(ans, names);

  free(ix);

  UNPROTECT(numprot);

  return(ans);

}




#undef BIG_POP_SIZE
#undef C2






