#ifndef _OPTIM0_H_
#define _OPTIM0_H_

SEXP OptimGa0(OptimFunc f, SEXP initPar, const long parlen, void *context,
              const int basePopSize, const long stopLags, const long minit,
              const long maxit, const int minimize, const double tol,
              const int relTol);

SEXP OptimGradient01(OptimFunc f, SEXP initpar, long parlen, void *context,
                     const double tol, const int relTol, const int minimize,
                     const long maxit);

#endif
