#ifndef URI_H_
#define URI_H_

#ifdef I_AM_MAIN
#define GLOBAL
#else
#define GLOBAL extern
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "nrutil2.h"

typedef double (*OptimFunc)(double *, long, void *);
typedef double (*optimFun)(SEXP, SEXP);

static const size_t DSIZE = sizeof(double);

#define SET_ELT(x, i, val) SET_ELEMENT(x, i, val)
#define GET_ELT(x, i)      VECTOR_ELT(x, i)
#define IS_MATRIX(x)       isMatrix(x)

#define isPosInf(x) (R_PosInf == (x))
#define isNegInf(x) (R_NegInf == (x))
#define IsPosInf    isPosInf
#define IsNegInf    isNegInf
#define IS_POS_INF  isPosInf
#define IS_NEG_INF  isNegInf
#define CHAR_PTR(x) CHARACTER_POINTER(x)
#define AS_CHAR(x)  AS_CHARACTER(x)
#define SQR(x)      ((x) * (x))
#define ABS(x)      ((x) < 0 ? -(x) : (x))
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#define SIGN_OF(x)  ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
GLOBAL long lmaxarg1, lmaxarg2;
#define LMAX(a, b) \
    ((lmaxarg1 = (long)(a)) > (lmaxarg2 = (long)(b)) ? lmaxarg1 : lmaxarg2)
GLOBAL long lminarg1, lminarg2;
#define LMIN(a, b) \
    ((lminarg1 = (long)(a)) < (lminarg2 = (long)(b)) ? lminarg1 : lminarg2)
GLOBAL int iminarg1, iminarg2;
#define IMIN(a, b) \
    ((iminarg1 = (int)(a)) < (iminarg2 = (int)(b)) ? iminarg1 : iminarg2)
GLOBAL double dmaxarg1, dmaxarg2;
#define DMAX(a, b) \
    ((dmaxarg1 = (double)(a)) > (dmaxarg2 = (double)(b)) ? dmaxarg1 : dmaxarg2)
GLOBAL double dsqrarg;
#define DSQR(a) ((dsqrarg = (double)(a)) == 0.0 ? 0.0 : dsqrarg * dsqrarg)
GLOBAL double dabsarg;
#define DABS(x) ((dabsarg = (double)(x)) >= 0.0 ? dabsarg : -dabsarg)
GLOBAL double dsignarg;
#define DSIGN(x) \
    ((dsignarg = (double)(x)) > 0.0 ? 1.0 : (dsignarg < 0.0 ? -1.0 : 0.0))
GLOBAL double lsignarg;
#define LSIGN(x) ((lsignarg = (long)(x)) > 0 ? 1 : (lsignarg < 0 ? -1 : 0))

#define NORM_DIST(x) pnorm(x, 0, 1, 1, 0)
#define NORM_DENS(x) dnorm(x, 0, 1, 0)

GLOBAL SEXP *ensSxpPtr;
GLOBAL long  ensLenArg;
#define PROT2(x, i) (((i) += 1), PROTECT(x))
#define UNPROT2     UNPROTECT(numprot)
#define ENSURE_NUMERIC(x, i)                                     \
    ((ensSxpPtr = &(x))                                          \
         ? (IS_NUMERIC(*ensSxpPtr)                               \
                ? *ensSxpPtr                                     \
                : PROT2(*ensSxpPtr = AS_NUMERIC(*ensSxpPtr), i)) \
         : R_NilValue)
#define ENSURE_CHAR(x, i)                                                    \
    ((ensSxpPtr = &(x)) ? (IS_CHARACTER(*ensSxpPtr)                          \
                               ? *ensSxpPtr                                  \
                               : PROT2(*ensSxpPtr = AS_CHAR(*ensSxpPtr), i)) \
                        : R_NilValue)
#define ENSURE_INTEGER(x, i)                                     \
    ((ensSxpPtr = &(x))                                          \
         ? (IS_INTEGER(*ensSxpPtr)                               \
                ? *ensSxpPtr                                     \
                : PROT2(*ensSxpPtr = AS_INTEGER(*ensSxpPtr), i)) \
         : R_NilValue)
#define ENSURE_LENGTH(x, n, i)                                                \
    ((ensSxpPtr = &(x)) ? (length(*ensSxpPtr) == (ensLenArg = (long)(n))      \
                               ? *ensSxpPtr                                   \
                               : PROT2(SET_LENGTH(*ensSxpPtr, ensLenArg), i)) \
                        : R_NilValue)

#ifndef HUGE
#define HUGE (1e100)
#endif

#define HUGE_NEGATIVE   (-HUGE)
#define BS_PTR_LEN      3
#define DIGITAL_PTR_LEN BS_PTR_LEN
enum
{ /* Used to be #define's */
  MARGRABE0_PRICE_INDEX,
  MARGRABE0_D1_INDEX,
  MARGRABE0_D2_INDEX,
  MARGRABE0_G1_INDEX,
  MARGRABE0_G2_INDEX,
  MARGRABE0_V1_INDEX,
  MARGRABE0_V2_INDEX,
  MARGRABE0_ETA_INDEX,
  MARGRABE0_PTR_LEN
};
#define MARGRABE0_RET_NAMES   {"p", "d1", "d2", "g1", "g2", "v1", "v2", "eta"}
#define MARGRABE2_0_PTR_LEN   3
#define MARGRABE2_0_RET_NAMES {"p", "d1", "d2"}
#define PEARSON_VEC_LIST_LEN  8
enum
{ /* Used to be #define's */
  PEARSON0_PRICE_INDEX,
  PEARSON0_D1_INDEX,
  PEARSON0_D2_INDEX,
  PEARSON0_PTR_LEN
};
#define PEARSON0_RET_NAMES {"p", "d1", "d2"}

#define QUIT_IF(x, y, z) \
    {                    \
        if (x)           \
        {                \
            (y) = (z);   \
            return;      \
        }                \
    }
#define SET_NA_TO_ZERO(x) (ISNA((double)(x)) ? (x) = 0, 1 : 0)

#define MALLOC(n, size) R_alloc(n, size)
#define FREE(x) \
    {           \
        ;       \
    }
#define MALLOC_CHECK() \
    {                  \
        ;              \
    }

#define CALL_PAYOFF(S, K) DMAX((S) - (K), 0.0)
#define PUT_PAYOFF(S, K)  DMAX((K) - (S), 0.0)
#define STRD_PAYOFF(S, K) DABS((S) - (K))

GLOBAL char optionPayoffArg;
#define OPTION_PAYOFF(S, K, type)    \
    (optionPayoffArg = (char)(type), \
     ('c' == optionPayoffArg         \
          ? CALL_PAYOFF(S, K)        \
          : ('p' == optionPayoffArg  \
                 ? PUT_PAYOFF(S, K)  \
                 : ('s' == optionPayoffArg ? STRD_PAYOFF(S, K) : NA_REAL))))

GLOBAL double stepFunArg;
#define STEP_FUN(x)              \
    ((stepFunArg = (double)(x)), \
     (stepFunArg > 0.0 ? 1.0 : (0.0 == stepFunArg ? NA_REAL : 0.0)))

#define CALL_EXP_DELTA(S, K) STEP_FUN((S) - (K))
#define PUT_EXP_DELTA(S, K)  (-STEP_FUN((K) - (S)))
#define STRD_EXP_DELTA(S, K) DSIGN((S) - (K))

GLOBAL char expDeltaArg;
#define EXP_DELTA(S, K, type)          \
    (expDeltaArg = (char)(type),       \
     ('c' == expDeltaArg               \
          ? CALL_EXP_DELTA(S, K)       \
          : ('p' == expDeltaArg        \
                 ? PUT_EXP_DELTA(S, K) \
                 : ('s' == expDeltaArg ? STRD_EXP_DELTA(S, K) : NA_REAL))))

#define MAX_DIM_LENGTH 4

#define VECTOR(x)     ((x).vec)
#define MATRIX(x)     ((x).mat)
#define ARRAY1(x)     ((x).vec)
#define ARRAY2(x)     ((x).mat)
#define ARRAY3(x)     ((x).arr3)
#define ARRAY4(x)     ((x).arr4)
#define DIM(x)        ((x).dim)
#define NROW(x)       ((x).dim[0])
#define NCOL(x)       ((x).dim[1])
#define DIM_LENGTH(x) ((x).ndim)

typedef struct Carray
{
    double    *vec;
    double   **mat;
    double  ***arr3;
    double ****arr4;
    int        dim[MAX_DIM_LENGTH];
    int        ndim;
} carray;

typedef struct uritm
{
    int year;
    int month;
    int day;
} uritm;

int     test_array_conform(const carray a1, const carray a2);
int     is_matrix(const carray a);
carray  make_array(double *vec, int *dim, int ndim);
carray  make_zero_array(int *dim, int ndim);
carray  make_matrix(double *vec, int nrow, int ncol);
carray  make_zero_matrix(int nrow, int ncol);
carray  make_identity_matrix(int n);
carray  make_vec(double *vec, const int len);
carray  make_zero_vec(const int len);
carray  subarray(carray a, int index);
double *column(carray a, int j);
long    vector_length(carray a);
void    set_array_to_zero(carray arr);
void    copy_array(carray orig, carray ans);
void    array_op(carray arr1, carray arr2, char op, carray ans);
carray  array_op_b(carray arr1, carray arr2, char op);
void    scalar_op(carray arr, double s, char op, carray ans);
void    transpose_matrix(carray mat, carray ans);
void    transpose_matrix_ns(carray mat, carray ans);
carray  transpose_matrix_b(carray);
void    transpose_matrix_ow(carray a);
void matrix_prod(carray mat1, carray mat2, int trans1, int trans2, carray ans);
void matrix_prod_ns(carray mat1, carray mat2, int trans1, int trans2,
                    carray ans);
carray matrix_prod_b(carray mat1, carray mat2, int trans1, int trans2);
void   matrix_prod_3(carray mat1, carray mat2, carray mat3, int trans1,
                     int trans2, int trans3, carray ans);
carray matrix_prod_3b(carray mat1, carray mat2, carray mat3, int trans1,
                      int trans2, int trans3);
void   matrix_prod_4(carray mat1, carray mat2, carray mat3, carray mat4,
                     int trans1, int trans2, int trans3, int trans4, carray ans);
void   matrix_prod_n(carray *mat, int *trans_ptr, int length, carray ans);
void   matrix_ptr_prod(carray *mat, int *trans_ptr, int length, carray ans);
int    matrix_inverse(carray a, carray ainverse);
carray matrix_inverse_b(carray, int *invertible);
int    matrix_inverse_2x2(carray, carray);
carray vec(carray);
void   kronecker(carray, carray, carray);
double sumsq(carray);
double det(carray);
carray sexp_to_carray(SEXP, int);
SEXP   carray_to_sexp(carray);
void   delta_to_strike(const int *n, const double *T, const double *S,
                       const double *vol, const double *del, const char **opt,
                       const double *r, double *q, double *ans);
SEXP   bvt_normuri(SEXP n, SEXP rho, SEXP jiggle);
SEXP   simulate_garch1(SEXP coefs, SEXP days, SEXP r0, SEXP h0);
SEXP   roll_apply(SEXP x, const SEXP fn, const SEXP order, const SEXP env);
void bs_value(const double T, const double S, const double K, const double vol,
              const double r, const double q, char opt, char ret, double *ans);
SEXP deltaToStrike(SEXP T, SEXP S, SEXP vol, SEXP del, SEXP optionType, SEXP r,
                   SEXP q);
SEXP bs(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP r, SEXP q, SEXP opt, SEXP ret);
int  bs_price_chkargs(const double T, const double S, const double K,
                      const double vol, const double r, const double q);
double bs_price(double T, double S, double K, double vol, double r, double q,
                char o);
double bs_delta(double T, double S, double K, double vol, double r, double q,
                char o);
double bs_gamma(double T, double S, double K, double vol, double r, double q,
                char o);
double bs_vega(double T, double S, double K, double vol, double r, double q,
               char o);
double bs_theta(double T, double S, double K, double vol, double r, double q,
                char o);
SEXP   maxDrawdown(SEXP x, SEXP incremental);
SEXP bsDeltaHedge1(SEXP T, SEXP S, SEXP K, SEXP vol, SEXP posn, SEXP r, SEXP q,
                   SEXP hedgePeriod, SEXP optionTypes, SEXP transCosts,
                   SEXP reltc);
SEXP bsDeltaHedge2(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP posn, SEXP r,
                   SEXP q, SEXP hedgePeriod, SEXP optionTypes, SEXP transCosts,
                   SEXP reltc, SEXP returnDaily);
SEXP bsDeltaHedge3(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP posn, SEXP r,
                   SEXP q, SEXP hedgePeriod, SEXP optionTypes, SEXP transCosts,
                   SEXP reltc, SEXP returnDaily);
SEXP dvHedge2(SEXP posn, SEXP posnStrikes, SEXP optionTypes, SEXP hedgePeriod,
              SEXP T, SEXP S, SEXP vols, SEXP volStrikes, SEXP r, SEXP q,
              SEXP undTransCosts, SEXP optTransCosts, SEXP isRelUndCosts,
              SEXP isRelOptCosts);
SEXP dvHedge(SEXP posn, SEXP posnStrikes, SEXP optionTypes, SEXP hedgePeriod,
             SEXP T, SEXP S, SEXP vols, SEXP volStrikes, SEXP r, SEXP q,
             SEXP undTransCosts, SEXP optTransCosts, SEXP isRelUndCosts,
             SEXP isRelOptCosts);
SEXP sprdOptDh1(SEXP T, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2,
                SEXP rho, SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP grdPts);
void gls_sumsq1(int T, int n, int p, double *y, double **X, double *beta,
                double *phi, double *sumsq, double *res);
SEXP gls_sumsq2(SEXP y, SEXP X, SEXP beta, SEXP phi);
void margrabe0(const double T, const double S1, const double S2,
               const double vol1, const double vol2, const double rho,
               const double r, const double q1, const double q2, const char opt,
               double *ans);
SEXP margrabe(SEXP tau, SEXP S1, SEXP S2, SEXP vol1, SEXP vol2, SEXP rho,
              SEXP r, SEXP q1, SEXP q2, SEXP opt);
SEXP margrabe2(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2,
               SEXP rho, SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP npaths,
               SEXP eps);
void pearson0(double tau, double S1, double S2, double K, double vol1,
              double vol2, double rho, double r, double q1, double q2,
              char optionType, long nGrdPts, int calcDeltas, double *valG);
SEXP pearson(SEXP tau, SEXP S1, SEXP S2, SEXP K, SEXP vol1, SEXP vol2, SEXP rho,
             SEXP r, SEXP q1, SEXP q2, SEXP opt, SEXP calcDeltas, SEXP grdPts);
SEXP optimGradient3(optimFun f, SEXP initPar, SEXP controlPar, const double tol,
                    const int relTol, const int minimize, const long maxit);
SEXP optimGa3(optimFun f, SEXP initPar, SEXP controlPar, SEXP parType,
              const int parLen, const int basePopSize, const long minit,
              const long stopLag, const long maxit, const int minimize,
              const double tol, const int relTol);
SEXP optim_gradient1(SEXP, SEXP, SEXP, SEXP, SEXP maxit, SEXP env);
SEXP optim_ga1(SEXP par, SEXP fn, SEXP gens, SEXP minimize, SEXP tol, SEXP env);
SEXP optimGradient2(optimFun f, SEXP initPar, SEXP controlPar, const double tol,
                    const int relTol, const int minimize, const long maxit);
SEXP optimGa2(optimFun f, SEXP initPar, SEXP controlPar, const int parLen,
              const int basePopSize, const long stopLags, const long minit,
              const long maxit, const int minimize, const double tol,
              const int relTol);
SEXP asRealMatrix(SEXP x);
double dmean(const double *x, long n);
double dvar(const double *x, long n);
SEXP   dSubMatrix(SEXP x, long row1, long row2, long col1, long col2);
SEXP   dSubVector(SEXP x, long ind1, long ind2);
void   dswap(double *x, double *y);
SEXP   numrep(SEXP x, const long n);
SEXP   intrep(SEXP x, const long n);
SEXP   charrep(SEXP x, const long n);
SEXP   listrep(SEXP x, const long n);
SEXP   rep(SEXP x, SEXPTYPE type, const long n);
SEXP   tolist(SEXP x);
void   cvect_to_string(SEXP src, char *dest);
void   dnormalize(double *x, long len);
void   dsetval(double *x, long len, double val);
#define dsetzero(x, len) dsetval(x, len, 0.0)
#define dsetna(x, len)   dsetval(x, len, NA_REAL)
SEXP na_omit(SEXP);
SEXP make_sexp_vec(const double *x, const int len);
void set_names(SEXP dest, char **names);
#define setNames set_names
double *matcol1(const SEXP mat, const int j);
SEXP    matcol2(const SEXP mat, const int j);
void    dmatrow1(SEXP mat, const int i, double *ans);
void    dmatrow2(SEXP mat, const long i, double **ans);
void    cmatrow1(const SEXP mat, const int i, char *ans);
void    cmatrow2(const SEXP mat, const int i, char **ans);
void    imatrow2(SEXP mat, const long i, int **ans);
void    double_swap(double **ptr1, double **ptr2);
double  sum_product(const double *x, const double *y, int n);
double  sum_array(const double *x, const int n);
int     imax(const int *x, const long len);
#define int_max imax
int    imin(const int *x, const long len);
long   lmax(const long *x, long len);
double double_max(const double *x, const int len);
#define dmax double_max
double double_min(const double *x, const int len);
#define dmin double_min
long   placeIndex(long place, double *x, long xlen);
double placeVal(long place, double *x, long xlen);
double nearestVal(double y, long place, double *x, long xlen);
double near(double refval, long place, double *refptr, double *outptr,
            long ptrlen);
long   nearestIndex(double y, long place, double *x, long xlen);
void   m_avg(double *y, int *order, int *n, double *ans);
void   m_max(double *y, int *order, int *n, double *ans);
void   m_min(double *y, int *order, int *n, double *ans);
void   get_random_sample(const double *in, double *out, const long inlen,
                         long outlen);
void   get_random_sample2(const double *src, double *out, const long srclen,
                          long outlen);
void   set_seed_uri(const unsigned long seed);
long   get_rand_uri(const long max);
void   dperm(double *x, const long n);
void   unif_rand_uri(double *u, const long n, int jiggle, int perm, int rep);
void   norm_rand_uri(double *z, const long n, int jiggle, int perm, int rep);
SEXP   num_deriv(SEXP par, SEXP fn, SEXP env);
SEXP   vector_norm(SEXP x);
SEXP   unit_vector(SEXP x);
SEXP   feval(SEXP par, SEXP fn, SEXP env);
SEXP   getListElt(SEXP list, const char *str);
SEXP   kalman1(SEXP yr, SEXP Hr_tr, SEXP Fr, SEXP Qr, SEXP Rr, SEXP xr,
               SEXP Ar_tr);
SEXP   sortedUniqueInteger(SEXP t);
int    iEltOf(long val, SEXP t, long *index);
void   sortWithIndex2(double *x, double *y, long n);
#define sort2 sortWithIndex2
int strptime_uri(const char *src, uritm *dest);
int uritm_lessThan(uritm *tm1, uritm *tm2);
int uritm_greaterThan(uritm *tm1, uritm *tm2);
int uritm_equals(uritm *tm1, uritm *tm2);
#define uritm_lessOrEqual(a, b) (uritm_lessThan(a, b) || uritm_equals(a, b))
#define uritm_greaterOrEqual(a, b) \
    (uritm_greaterThan(a, b) || uritm_equals(a, b))
long uritm_difftime(uritm *tm1, uritm *tm2);
SEXP cum_cov(SEXP x, SEXP y, SEXP reverse);
SEXP cum_var(SEXP x, SEXP reverse);
void drev(long len, const double *in, double *out);
void dmult(long len, const double *in1, const double *in2, double *out);
SEXP cum_corr(SEXP x, SEXP y, SEXP reverse);

#endif /* URI_H */
