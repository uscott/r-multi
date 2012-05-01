
#include "uri.h"


static void cum_cov_chk(SEXP *x, SEXP *y, int *numprot)
{
  ENSURE_NUMERIC(*x, *numprot);
  ENSURE_NUMERIC(*y, *numprot);
  if (length(*x) != length(*y))
    error ("Unequal lengths.");
}

SEXP cum_cov(SEXP x, SEXP y, SEXP reverse)
{
  int numprot = 0;
  const int rev = asInteger(reverse);
  long i, n;
  double *xy = 0, *xr = 0, *yr = 0, *px = 0, *py = 0;
  double exy, ex, ey;
  SEXP ans = R_NilValue;

  cum_cov_chk(&x, &y, &numprot);

  n = length(x);
  if (!n) {
    UNPROT2;
    return ans;
  }

  if (0 == (xy = (double *) malloc(n * DSIZE)))
    error ("Malloc failed.");
  PROT2(ans = NEW_NUMERIC(n), numprot);

  if (rev) {
    if (0 == (xr = (double *) malloc(n * DSIZE)))
      error ("Malloc failed.");
    if (0 == (yr = (double *) malloc(n * DSIZE))) {
      free(xr);
      error ("Malloc failed.");
    }
    drev(n, REAL(x), xr);
    drev(n, REAL(y), yr);
  }

  px = rev ? xr : REAL(x);
  py = rev ? yr : REAL(y);
  dmult(n, px, py, xy);

  for (i = 0; i < n; i++) {
    ex = dmean(px++, n - i);
    ey = dmean(py++, n - i);
    exy = dmean(xy + i, n - i);
    REAL(ans)[rev ? n - i - 1 : i] = exy - ex * ey;
  }

  free(xy);
  if (rev) {free(xr); free(yr);}
  SET_NAMES(ans, isNull(GET_NAMES(x)) ? GET_NAMES(y) : GET_NAMES(x));
  UNPROT2;
  return ans;
}
