#include "uri.h"
#include <assert.h>
#include <time.h>



SEXP asRealMatrix( SEXP x )
{
    int  numprot = 0;
    long xlen    = 0;
    SEXP 
        y        = R_NilValue, 
        dim      = R_NilValue;

    if ( IS_MATRIX( x ))
        return x;

    ENSURE_NUMERIC( x, numprot );
    xlen = length(x);
 
    PROT2( dim = NEW_INTEGER( 2 ),    numprot);
    PROT2( y   = NEW_NUMERIC( xlen ), numprot);

    REAL( dim )[ 0 ] = xlen;
    REAL( dim )[ 1 ] = 1;
    
    memcpy( REAL(y), REAL(x), xlen * sizeof( double ));
    
    setAttrib( y, R_DimSymbol, dim );

    UNPROTECT( numprot );

    return y;
}



SEXP dSubMatrix(SEXP x, long row1, long row2, long col1, long col2)
{
    int ok, numprot=0;
    long j, mx, nx, mAns, nAns;
    double *xPtr, *ansPtr;
    SEXP ans;

    ENSURE_NUMERIC(x,numprot);

    mx = nrows(x);
    nx = ncols(x);
    mAns = row2 - row1 + 1;
    nAns = col2 - col1 + 1;

    ok = isMatrix(x) && 
        0 <= row1 && row1 <= row2 && row2 < mx && 
        0 <= col1 && col1 <= col2 && col2 < nx;

    if (!ok)
        error ("bad args dSubMatrix");

    PROT2( ans = allocMatrix( REALSXP, mAns, nAns ), numprot );

    for (j = 0; j < nAns; j++) 
    {
        ansPtr = REAL(ans) + j * mAns;
        xPtr = REAL(x) + (col1 + j) * mx + row1;
        memcpy(ansPtr, xPtr, mAns * sizeof(double));
    }

    UNPROT2;

    return ans;
}



SEXP dSubVector(SEXP x, long ind1, long ind2)
{
    int numprot=0;
    SEXP ans;

    ENSURE_NUMERIC(x,numprot);

    if ( !(0 <= ind1 && ind1 <= ind2 && ind2 < length(x)) )
        error ("bad args dSubVector");
  
    PROT2(ans = NEW_NUMERIC(ind2 - ind1 + 1),numprot);

    memcpy(REAL(ans), REAL(x) + ind1, length(ans) * sizeof(double));

    UNPROT2;

    return ans;
}



void dnormalize(double *x, long len)
{
    register long i;
    register double norm = 0.0; 
    register double *y;

    for (i = 0, y = x; i < len; i++, y++)
        norm += SQR(*y);

    norm = sqrt(norm);

    for (i = 0; i < len; i++)
        *x++ /= norm;
}



double dmean(const double *x, long n)
{
    register long i;
    register double mean = 0.0;

    for ( i = 0; i < n; i++ )
        mean += *x++;

    return mean / n;
}



double dvar(const double *x, long n)
{
    register long i;
    register double var = 0.0, m;

    for (i = 0; i < n; i++)
        var += DSQR(x[i]);

    var /= (double)n;

    m = dmean( x, n );

    return var - m * m;
}



void dswap(double *x, double *y)
{
    double tmp;

    tmp = *x;
    *x  = *y;
    *y  = tmp;
}



SEXP numrep(SEXP x, const long n)
{
    long i, j, m, numprot=0;

    if (n < 1)
        return R_NilValue;

    ENSURE_NUMERIC(x,numprot);

    m = length(x);

    PROT2( SET_LENGTH( x, m * n ), numprot);

    for (j = 1; j < n; j++)
        for (i = 0; i < m; i++)
            REAL(x)[i + m * j] = REAL(x)[i];

    UNPROT2;

    return x;
}



SEXP intrep(SEXP x, const long n)
{
    long i, j, m;

    if (n < 1)
        return R_NilValue;

    PROTECT(x = AS_INTEGER(x));

    m = length(x);

    PROTECT(SET_LENGTH(x, m * n));
    
    for (j = 1; j < n; j++)
        for (i = 0; i < m; i++)
            INTEGER(x)[i + m * j] = INTEGER(x)[i];

    UNPROTECT(2);

    return x;
}



SEXP charrep(SEXP x, const long n)
{
  long i, j, m;

  if (n < 1)
    return R_NilValue;

  PROTECT(x = AS_CHAR(x));

  m = length(x);

  PROTECT(SET_LENGTH(x, m * n));

  for (j = 1; j < n; j++)
    for (i = 0; i < m; i++)
      CHAR_PTR(x)[i + m * j] = CHAR_PTR(x)[i];

  UNPROTECT(2);

  return x;
  
}


SEXP listrep(SEXP x, const long n)
{
  long i, j, m;

  if (n < 1)
    return R_NilValue;

  PROTECT(x = AS_LIST(x));

  m = length(x);

  PROTECT(SET_LENGTH(x, m * n));

  for (j = 1; j < n; j++)
    for (i = 0; i < m; i++)
      SET_ELT(x, i + m * j, GET_ELT(x, i));

  UNPROTECT(2);

  return x;
}



SEXP rep(SEXP x, SEXPTYPE type, const long n)
{
  SEXP ans = R_NilValue;
  
  switch(type) {
    
  case REALSXP:
    PROTECT(ans = numrep(x, n));
    break;
  case INTSXP:
    PROTECT(ans = intrep(x, n));
    break;
  case STRSXP:
    PROTECT(ans = charrep(x, n));
    break;
  case VECSXP:
    PROTECT(ans = listrep(x, n));
    break;
  default:
    PROTECT(ans = R_NilValue);
    warning ("unimplemented SEXPTYPE in rep");
    break;

  }

  UNPROTECT(1);

  return ans;
}


SEXP tolist(SEXP x)
{
  SEXP ans;

  if (isNewList(x))
    return x;

  ans = PROTECT(NEW_LIST(1));

  SET_ELT(ans, 0, x);

  UNPROTECT(1);

  return ans;

}



void cvect_to_string(SEXP src, char *dest)
{
  const int n = length(src);
  int i;
 
  if (!IS_CHARACTER(src))
    error ("!IS_CHARACTER(src)");
 
  for (i = 0; i < n; i++)
    *dest++ = *CHAR(GET_ELT(src, i));

}


void dsetval(double *x, long len, double val)
  /*******************************************************************
   *
   *  Description: Sets x[0], ..., x[len - 1] to val.
   *
   *******************************************************************/
{
  while (len-- > 0)
    *x++ = val;
}



SEXP na_omit(SEXP x)

  /*******************************************************************
   *
   *  Description: Returns AS_NUMERIC(x) with NA's removed.
   *
   *******************************************************************/

{
  long  alen, xlen;
  register long    i;
  register double *xp, *ap;
  SEXP ans;

  PROTECT(x = AS_NUMERIC(x));

  xlen = length(x);

  PROTECT(ans = NEW_NUMERIC(xlen));
  alen = 0;
  xp   = REAL(x);
  ap   = REAL(ans);

  for (i = 0; i < xlen; i++, xp++)
    if (!ISNAN(*xp)) {
      *ap++ = *xp;
      alen += 1;
    }

  UNPROTECT( 2 );

  return SET_LENGTH( ans, alen );
}



SEXP make_sexp_vec(const double *x, const int len)
{
    SEXP ans;

    PROTECT( ans = NEW_NUMERIC(len) );

    memcpy( REAL( ans ), x, len * sizeof( double ));

    UNPROTECT(1);

    return ans;
}



void set_names(SEXP dest, char **names)
{
  int i, n;
  SEXP nameS;

  n     = length(dest);
  PROTECT(nameS = NEW_STRING(n));

  for (i = 0; i < n; i++)
    CHAR_PTR(nameS)[i] = mkChar(*names++);

  SET_NAMES(dest, nameS);

  UNPROTECT(1);
}



double *matcol1( const SEXP mat, const int j )
/********************************************************************
 *
 *   Description: Returns pointer to jth column of matrix mat
 *     (where the column indexing starts at 0).
 *
 ********************************************************************/
{
    double *ans = NULL;

    if ( !isMatrix( mat ) || j < 0 || ncols( mat ) <= j )
        error( "matcol1 :  Non-matrix or index out of range" );

    ans = (double *) ( REAL( mat ) + nrows( mat ) * j);

    return ans;
}


SEXP matcol2(const SEXP mat, const int j)
/*******************************************************************
 *
 *  Description: Returns a numeric SEXP that is a *COPY* of the jth
 *    column of mat.  Column indexing starts at 0.
 *
 *******************************************************************/
{
    int m;
    SEXP ans;

    m = nrows(mat);

    PROTECT( ans = NEW_NUMERIC( m ));

    memcpy( REAL(ans), REAL(mat) + m * j, m * sizeof(double));

    UNPROTECT( 1 );

    return ans;
}



void dmatrow1(SEXP mat, const int i, double *ans)
/*******************************************************************
 *
 *  Description: Puts the values of the i^th row of numeric matrix mat
 *    into ans[0], ..., ans[ncols(mat) - 1].  Row indexing starts
 *    at 0.
 *
 *******************************************************************/
{
    const int m = nrows(mat);
    const int n = ncols(mat);
    int j;
    double *p;

    if (!isMatrix(mat) || i < 0 || m <= i)
        error ("bad argument passed to dmatrow1");

    for (j = 0, p = REAL(mat) + i; j < n; j++, p += m)
        *ans++ = *p;
}



void dmatrow2(SEXP mat, const long i, double **ans)

  /*******************************************************************
   *
   *  Description: Makes ans[j] point to the jth element of the
   *    ith row of mat.
   *
   *******************************************************************/

{
    const long m = nrows(mat);
    const long n = ncols(mat);
    long j;
    double *p;

    if (!isMatrix(mat) || i < 0 || m <= i)
        error ("bad argument passed to dmatrow1");

    for (j = 0, p = REAL(mat) + i; j < n; j++, p += m, ans++)
        *ans = p;
}



void cmatrow1(const SEXP mat, const int i, char *ans)
/*******************************************************************
 *
 *  Description: Puts the values of the i^th row of *character* mat
 *    into ans[0], ..., ans[ncols(mat) - 1].  Row indexing starts
 *    at 0.
 *
 *******************************************************************/
{
    const int m = nrows(mat);
    const int n = ncols(mat);
    int j;

    if (!isMatrix(mat) || i < 0 || m <= i)
        error ("bad argument passed to cmatrow1");

    for (j = 0; j < n; j++)
        *ans++ = *CHAR(STRING_ELT(mat, m * j));
}



void cmatrow2(const SEXP mat, const int i, char **ans)

  /*******************************************************************
   *
   *  Description: Puts the values of the i^th row of *character*
   *    matrix mat into 
   *    *ans[ 0 ], ..., *ans[ ncols( mat ) - 1 ].  
   *    Row indexing starts at 0.
   *
   *******************************************************************/

{
    const int m = nrows(mat);
    const int n = ncols(mat);
    int j;

    if( !isMatrix( mat ) || i < 0 || m <= i )
        error ( "bad argument passed to cmatrow1" );

    for( j = 0; j < n; ++j )
        *ans++ = CHAR( STRING_ELT( mat, m * j ));
}



void imatrow2(SEXP mat, const long i, int **ans)

  /*******************************************************************
   *
   *  Description: Makes ans[j] point to the jth element of the
   *    ith row of mat.
   *
   *******************************************************************/

{

  const long m = nrows(mat);
  const long n = ncols(mat);
  long j;
  int *p;

  if (!isMatrix(mat) || i < 0 || m <= i)
    error ("bad argument passed to dmatrow1");

  for (j = 0, p = INTEGER(mat) + i; j < n; j++, p += m, ans++)
    *ans = p;

}



SEXP getListElt(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  long i;

  if (!isNewList(list))
    error ("bad args getListElt");

  for (i = 0; i < length(list); i++)
    if (!strcmp(CHAR(STRING_ELT(names, i)), str)) {

      elmt = GET_ELT(list, i);
      break;

    }

  return elmt;

}



SEXP feval(const SEXP fn, const SEXP par, const SEXP env)
{
  SEXP R_fcall;

  if (!isFunction(fn) || !isNumeric(par) || !isEnvironment(env))
    error("C error: feval");

  PROTECT(R_fcall = lang2(fn, R_NilValue));

  SETCADR(R_fcall, par);

  UNPROTECT(1);

  return eval(R_fcall, env);
}


int imax(const int *x, const long n)

  /********************************************************************
  *
  *   Description: Returns the maximum of x[0], ..., x[n-1].
  *
  ********************************************************************/

{
  register int ans = *x;
  int const *end = x + n;

  if (n <= 0)
    error ("n <= 0 passed to int_max");

  while (++x < end)
    ans = ans < *x ? *x : ans;

  return ans;
}



int imin(const int *x, const long n)

  /********************************************************************
  *
  *   Description: Returns the maximum of x[0], ..., x[n-1].
  *
  ********************************************************************/

{
  register int ans = *x;
  int const *end = x + n;

  if (n <= 0)
    error ("n <= 0 passed to int_max");

  while (++x < end)
    ans = ans > *x ? *x : ans;

  return ans;
}



long lmax(const long *x, long len)

  /*******************************************************************
   *
   *  Description: Returns maximum of x[0], ..., x[len - 1].
   *
   *******************************************************************/

{
  long ans = 0;

  if (len < 0)
    error ("negative length passed to lmax");
  
  while (len--) {

    ans = MAX(ans, *x);
    ++x;

  }

  return ans;

}


double double_max(const double *x, const int n)

  /********************************************************************
  *
  *   Description: Returns the maximum of x[0], ..., x[n - 1].
  *
  ********************************************************************/

{
  register double ans = *x;
  double const    *xn = x + n;

  if (n <= 0)
    error ("n <= 0 passed to double_max");

  while (++x < xn)
    ans = ans < *x ? *x : ans;

  return ans;
}


double double_min(const double *x, const int n)
{
  register double ans = *x;
  double const    *xn = x + n;

  if (n <= 0)
    error ("n <= 0 passed to double_min");



  while (++x < xn)
    ans = ans > *x ? *x : ans;

  return ans;
}



void m_avg(double *y, int *order, int *n, double *ans)

  /********************************************************************
  *
  *   Description: Computes the moving average of order (*order) of
  *     y[0],...,y[n-1].  Puts result in double *ans.
  *
  ********************************************************************/

{
  int s, t;
  double *yPtr;
  const int new_order = *order - 1; 
 
  for (t = 0; t < new_order; t++)
    *ans++ = NA_REAL;

  y += new_order;
  for (t = new_order; t < (*n); t++) {
    yPtr = y++;
    *ans = *yPtr;

    for (s = new_order; s > 0; s--)
      *ans += *(--yPtr);

    *ans++ /= (new_order + 1);
  }
}


void m_max(double *y, int *order, int *n, double *ans)

  /********************************************************************
  *
  *   Description: Computes moving max of order (*order) of y[0],...,
  *     y[n-1] and puts result in double *ans.
  *
  ********************************************************************/


{
  int s, t;
  double *yPtr;
  const int new_order = *order - 1; 
 
  for (t = 0; t < new_order; t++)
    *ans++ = NA_REAL;

  y += new_order;

  for (t = new_order; t < (*n); t++) {
    yPtr = y++;
    *ans = *yPtr--;

    for (s = new_order; s > 0; s--) {
      *ans = MAX(*yPtr, *ans);
      --yPtr;
    }

    ++ans;
  }
}



void m_min(double *y, int *order, int *n, double *ans)
{
  int s, t;
  double *yPtr;
  const int new_order = *order - 1; 
 
  for (t = 0; t < new_order; t++)
    *ans++ = NA_REAL;

  y += new_order;

  for (t = new_order; t < (*n); t++) {
    yPtr = y++;
    *ans = *yPtr--;

    for (s = new_order; s > 0; s--) {
      *ans = MIN(*yPtr, *ans);
      --yPtr;
    }

    ++ans;
  }
}







long placeIndex(long n, double *x, long xlen)

  /*******************************************************************
   *
   *  Description: Returns i such that x[i] is in the n-th
   *    lowest value among x[0], ..., x[xlen - 1].
   *
   *******************************************************************/

{

  long i = 0, ans = NA_INTEGER;
  int *ind = (int *) MALLOC(xlen, sizeof(int));
  double *xCopy = (double *) MALLOC(xlen, sizeof(double));

  if (n > xlen - 1 || n < 0)
    error("bad inputs to C function which_index");
	
  memcpy(xCopy, x, xlen * sizeof(double));

  for(i = 0; i < xlen; i++) 
    ind[i] = i;

  rsort_with_index(xCopy, ind, xlen);

  ans = ind[n];

  FREE(ind);
  FREE(xCopy);

  return ans;
}



double placeVal(long n, double *x, long xlen)
  /* Last modified 26-Aug-2002.
     Returns the value x[i] 
     such that x[i] is in the mth place (starting at 0!)
     after sorting i.e. x[i] is the mth-lowest value in x.
     xlen is the length of [x]. */
{
  double ans = NA_REAL, *xCopy = (double *) MALLOC(xlen, sizeof(double));

  if (n > xlen - 1 || n < 0)
    error("bad inputs to C function which1");
	
  memcpy(xCopy, x, xlen*sizeof(double));

  R_rsort(xCopy, xlen);

  ans = xCopy[n];
  FREE(xCopy);

  return ans;
}




double nearestVal(double y, long n, double *x, long xlen)
  /* Finds the value x[i] such that x[i] is closest to y.
     long n is the length of x */
{
  double ans = NA_REAL, *dist = (double*) MALLOC(xlen, sizeof(double));
  long i = 0;

  for (i = 0; i < xlen; i++) 
    dist[i] = ABS(x[i] - y);

  ans = x[placeIndex(n, dist, xlen)];

  FREE(dist);

  return ans;
}


long nearestIndex(double y, long n, double *x, long xlen)
  /* Returns the index i such that x[i] is closest to y.
     long n is the length of x */
{

  double *dist = (double *) MALLOC(xlen, sizeof(double));
  long i = 0, ans = NA_INTEGER;

  for (i = 0; i < xlen; i++) 
    dist[i] = ABS(x[i] - y);

  ans = placeIndex(n, dist, xlen);

  FREE(dist);

  return ans;

}




double near(double refval, long n, double *refptr,
	    double *ans, long ptr_length)
  /* finds i such that refptr[i] is n-th closest to refval, then
     returns ans[i] */
{
  return ans[nearestIndex(refval, n, refptr, ptr_length)];
}



void double_swap(double **ptr1, double **ptr2)
{
  double *tmp;

  tmp = *ptr1;
  *ptr1 = *ptr2;
  *ptr2 = tmp;
}






double sum_product(const double *x, const double *y, int n)
{
  register double sum = 0;

  while (n-- > 0) 
    sum += (*x++) * (*y++);

  return sum;
}




double sum_array(const double *x, const int n)
{
  register double sum = *x;
  double const    *xn = x + n;

  while (++x < xn) 
    sum += *x;

  return sum;
}


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Stuff for random number generation follows, some of which is
**  adapted from NR in C.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



static unsigned long idum_uriel_scott = 1;

void set_seed_uri(const unsigned long seed)
{
  idum_uriel_scott = seed;
}


#define A1 1664525
#define A0 1013904223
long get_rand_uri(const long max)

  /*******************************************************************
   *
   *  Description: Returns a long between 0 and max - 1.
   *
   *******************************************************************/

{
  idum_uriel_scott = A1 * idum_uriel_scott + A0;

  return (long) (idum_uriel_scott % max);
}
#undef A1
#undef A0




void dperm(double *x, const long n)

  /********************************************************************
  *
  * Description: Randomly permutes x[0],...,x[n-1].
  *
  ********************************************************************/

{
  register long     j = 0;
  register double tmp = 0;
  double         *x0  = x;
  double const   *xn  = x + n;

  set_seed_uri((unsigned long) clock());

  for ( ; x < xn; x++) {

    j = get_rand_uri(n);

    tmp   = x0[j];
    x0[j] = *x;
    *x    = tmp;

  }
}




void unif_rand_uri(double *u, const long n, int jiggle, int perm, int rep)
{
  register int   i  = 0;
  const double  ni  = 1.0 / (double)n;
  double *v = NULL;


  GetRNGstate();

  for (i = 0; i < n; i++)
    *u++ = ni * ((double)i + (jiggle ? unif_rand() : 0.5));

  u -= n;

  PutRNGstate();

  if (perm) {
    
    if (!rep)
      dperm(u, n);
    else {

      v = (double *) MALLOC(n,sizeof(double));
      
      if (!v) error ("malloc failed in unif_rand_uri");

      get_random_sample(u, v, n, n);

      memcpy(u, v, n * sizeof(double));
      
      FREE(v);

    }

  }

}



void norm_rand_uri(double *z, const long n, int jiggle, int perm, int rep)
{
  double const *zn = z + n;

  unif_rand_uri(z, n, jiggle, perm, rep);

  for ( ; z < zn; z++)
    *z = qnorm(*z,/*mu*/0.0,/*sigma*/1.0,/*lower*/1,/*log*/ 0);


}





void get_random_sample(const double *src, double *out, 
                       const long srclen, long outlen)

  /*******************************************************************
   *
   *  Description: Uniformly sample out[0], ..., out[outlen - 1]
   *    from src[0], ..., src[srclen - 1].  
   *
   *******************************************************************/

{
  if (out && (!src || srclen <= 0))
    error ("get_random_sample:  Cannot sample from NULL or length 0");

  while (outlen-- > 0)
    *out++ = src[get_rand_uri(srclen)];
}



void get_random_sample2(const double *src, double *out, 
                        const long srclen, long outlen)

  /*******************************************************************
   *
   *  Description: Uniformly sample out[0], ..., out[outlen - 1]
   *    from src[0], ..., src[srclen - 1] unless src == NULL or
   *    0 == srclen, in which case sampling occurs from
   *    standard normal distribution.
   *
   *******************************************************************/

{
  if (!src || srclen <= 0) {

    GetRNGstate();

    while (outlen-- > 0)
      *out++ = norm_rand();

    PutRNGstate();

  }
  else
    get_random_sample(src, out, srclen, outlen);


}






SEXP vector_norm(SEXP x)
  // Last modified 7/22/02.
{
  int n;
  double *xPtr;
  SEXP norm;

  if(!isNumeric(x))
    error("non-numeric vector input to C function vector_norm");

  PROTECT(x = AS_NUMERIC(x));
  PROTECT(norm = NEW_NUMERIC(1));

  *REAL(norm) = 0.0;
  xPtr = REAL(x);
  n = length(x);

  while (n-- && R_FINITE(*REAL(norm))) {
    *REAL(norm) += (*xPtr) * (*xPtr);
    ++xPtr;
  }

  if(R_FINITE(*REAL(norm))) *REAL(norm) = sqrt(*REAL(norm));

  UNPROTECT(2);
  return(norm);
}

SEXP unit_vector(SEXP x) 
  // Last modified 7/22/02.
{
  int n;
  double *xPtr, *ansPtr;
  double norm;
  SEXP ans;

  if(!isNumeric(x)) 
    error("non-numeric vector input to C function unit_vector");

  n = length(x);
  PROTECT(x = AS_NUMERIC(x));
  PROTECT(ans = NEW_NUMERIC(n));

  norm = *REAL(vector_norm(x));
  xPtr = REAL(x);
  ansPtr = REAL(ans);

  if(R_FINITE(norm))
    while (n--) *ansPtr++ = *xPtr++ / norm;
  else 
    while (n--) *ansPtr++ = NA_REAL;

  UNPROTECT(2);
  return(ans);
}








SEXP num_deriv( SEXP fn, SEXP par, SEXP env )

  /*******************************************************************
   *
   *  Description: Computes numeric derivative of function fn
   *    at the value par.
   *
   *******************************************************************/
{
    const int 
        n = length(par);
    int 
        numprot = 0; /* Counts number of PROTECTs called */
    int 
        i, ok;
    double 
        eps = sqrt( DOUBLE_EPS );
    double 
        val, valnew;
    double 
        *par_ptr, *par_new_ptr, *grad_ptr;

    SEXP 
        gradient = R_NilValue, R_fcall, par_new;

    if( !isFunction( fn )) 
        error( "Bad function input in feval" );
    if( !isEnvironment( env )) 
        error( "Bad environment input in feval" );
    if( !isNumeric( par )) 
        error( "Bad parameters input in feval" );

    PROT2( par     = AS_NUMERIC( par ),       numprot );
    PROT2( R_fcall = lang2( fn, R_NilValue ), numprot );

    SETCADR( R_fcall, par );

    val = *REAL( eval( R_fcall, env ));
  
    ok = R_FINITE( val );

    if ( ok ) 
    {
        PROT2( gradient = NEW_NUMERIC( n ), numprot);
        PROT2( par_new  = NEW_NUMERIC( n ), numprot);

        par_ptr = REAL(par);
        par_new_ptr = REAL(par_new);
        grad_ptr = REAL(gradient);

        memcpy( par_new_ptr, par_ptr, n * sizeof( double ));

        for( i = 0; i < n; i++, ++grad_ptr, ++par_new_ptr)
        {
            *par_new_ptr += eps;

            SETCADR( R_fcall, par_new );
    
            valnew = *REAL( eval(R_fcall, env ));

            ok = R_FINITE( valnew );

            *grad_ptr = ok ? (valnew - val) / eps : NA_REAL;

            *par_new_ptr -= eps;
        }

    }

    UNPROTECT( numprot );

    return( gradient );
}





