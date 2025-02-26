#include <stdio.h>

#include "uri.h"

#define SWAP(a, b)   \
    {                \
        temp = (a);  \
        (a)  = (b);  \
        (b)  = temp; \
    }
#define TINY       1.0e-20
#define USE_MALLOC 1

static void assert2(int condition)
{
    if (!condition) error("assert2 failed in carray.c");
}

static carray init_array()
{
    int    i;
    int   *dim_ptr;
    carray a;

    /* Initialize everything to zero.  Useful for debugging */
    ARRAY1(a) = (double *)'\0';
    ARRAY2(a) = (double **)'\0';
    ARRAY3(a) = (double ***)'\0';
    ARRAY4(a) = (double ****)'\0';

    dim_ptr = DIM(a);

    for (i = 0; i < MAX_DIM_LENGTH; i++) *dim_ptr++ = 0;

    DIM_LENGTH(a) = 0;

    return a;
}

// Gauss-Jordan elimination
// The arguments a, b are destroyed with the matrix double **a being replaced by
// its inverse.

static int gaussj2(double **a, int n, double **b, int m)
{
    int     *ixc, *ixr, *ipiv, *ixc2, *ixr2, *ipiv2, *ipiv3;
    int      i, ic = 0, ir = 0, j, k, l, ll, invertible;
    double   big, dum, pivinv, temp;
    double  *ai1, *bi1, *ai2, *bi2, *aij;
    double **a2;

    if (USE_MALLOC)
    {
        ixc2 = ixc = (int *)malloc(n * sizeof(int));
        ixr2 = ixr = (int *)malloc(n * sizeof(int));
        ipiv2 = ipiv3 = ipiv = (int *)malloc(n * sizeof(int));

        if (!ixc || !ixr || !ipiv) error("Malloc failed.");
    } else
    {
        ixc2 = ixc = (int *)R_alloc(n, sizeof(int));
        ixr2 = ixr = (int *)R_alloc(n, sizeof(int));
        ipiv2 = ipiv3 = ipiv = (int *)R_alloc(n, sizeof(int));
    }

    for (j = 0; j < n; j++) *ipiv2++ = 0;

    for (i = 0, invertible = 1; i < n && invertible; i++)
    {
        ipiv2 = ipiv;
        big   = 0.0;
        a2    = a;

        for (j = 0; j < n; j++, a2++) /* Note assignment: */
            if (*ipiv2++ < 1 && (ipiv3 = ipiv, ai1 = *a2))
            {
                for (k = 0; k < n; k++)
                {
                    if (!(*ipiv3++))
                        if (ABS(*ai1) >= big)
                        {
                            ir  = j;
                            ic  = k;
                            big = ABS(*ai1);
                        }
                }
            }

        ++(ipiv[ic]);

        if (ir != ic)
        {
            ai1 = a[ir], ai2 = a[ic];
            bi1 = b[ir], bi2 = b[ic];
            for (l = 0; l < n; l++, ai1++, ai2++) SWAP(*ai1, *ai2);
            for (l = 0; l < m; l++, bi1++, bi2++) SWAP(*bi1, *bi2);
        }

        *ixr2++ = ir;
        *ixc2++ = ic;

        ai1 = a[ic];
        bi1 = b[ic];

        aij = ai1 + ic;

        invertible = !(0.0 == *aij);

        pivinv = 1.0 / (*aij);
        *aij   = 1.0;

        for (l = 0; l < n; l++) *ai1++ *= pivinv;
        for (l = 0; l < m; l++) *bi1++ *= pivinv;

        for (ll = 0; ll < n; ll++)
            if (ll != ic)
            {
                ai1 = a[ll], ai2 = a[ic];
                bi1 = b[ll], bi2 = b[ic];
                dum     = ai1[ic];
                ai1[ic] = 0.0;

                for (l = 0; l < n; l++) *ai1++ -= *ai2++ * dum;
                for (l = 0; l < m; l++) *bi1++ -= *bi2++ * dum;
            }
    }

    ixr2 = ixr + n - 1;
    ixc2 = ixc + n - 1;

    if (invertible)
        for (l = 0; l < n; l++, ixr2--, ixc2--)
        {
            if (*ixr2 != *ixc2 && (a2 = a))
            { /* Note assignment */
                for (k = 0; k < n; k++, a2++) SWAP((*a2)[*ixr2], (*a2)[*ixc2]);
            }
        }

    if (USE_MALLOC)
    {
        free(ipiv);
        free(ixr);
        free(ixc);
    }

    return invertible;
}

static int ludcmp2(double **a, int n, int *indx, double *d)

/*******************************************************************
 *
 *  Description: LU Decomposition from based on NR in C.  The matrix a
 *    is destroyed in the process and replaced by its LU
 *    decomposition.
 *
 *******************************************************************/

{
    int     i, imax = 0, j, k, invertible;
    double  big, dum, sum, temp;
    double *vv;

    vv = dvector(0, n - 1);
    *d = 1.0;

    for (i = 0, invertible = 1; i < n && invertible; i++)
    {
        big = 0.0;

        for (j = 0; j < n; j++)
            if ((temp = ABS(a[i][j])) > big) big = temp;

        invertible = !(0.0 == big);
        vv[i]      = 1.0 / big;
    }

    for (j = 0; j < n && invertible; j++)
    {
        for (i = 0; i < j; i++)
        {
            sum = a[i][j];
            for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }

        big = 0.0;
        for (i = j; i < n; i++)
        {
            sum = a[i][j];
            for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];

            a[i][j] = sum;

            if ((dum = vv[i] * ABS(sum)) >= big)
            {
                big  = dum;
                imax = i;
            }
        }

        if (j != imax)
        {
            for (k = 0; k < n; k++)
            {
                dum        = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k]    = dum;
            }

            *d       = -(*d);
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (0.0 == a[j][j]) a[j][j] = TINY;

        if (j < n - 1)
        {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i < n; i++) a[i][j] *= dum;
        }
    }

    free_dvector(vv, 0, n - 1);

    return invertible;
}

long vector_length(carray a)
{
    int       *dimp  = DIM(a);
    int const *dimpE = dimp + DIM_LENGTH(a);
    long       len;

    for (len = 1; dimp < dimpE; dimp++) len *= *dimp;

    return len;
}

carray make_array(double *vec, int *dim, int ndim)
{
    int    d, i, j, *dimp;
    long   len[MAX_DIM_LENGTH + 1];
    carray a;

    if (ndim > MAX_DIM_LENGTH) error("Too many dimensions in make_array");

    a = init_array();

    len[ndim] = 1;
    dimp      = dim;
    for (d = ndim; d >= 1; d--) len[d - 1] = len[d] * (*dimp++);

    for (d = 1; d <= ndim; d++)
    {
        switch (d)
        {
            case 1:
                VECTOR(a) = vec;
                break;
            case 2:
                ARRAY2(a) = (double **)R_alloc(len[2 - 1], sizeof(double *));
                for (i = 0, j = 0; i < len[2 - 1]; i++, j += dim[ndim - 2 + 1])
                {
                    ARRAY2(a)[i] = ARRAY1(a) + j;
                }
                break;
            case 3:
                ARRAY3(a) = (double ***)R_alloc(len[3 - 1], sizeof(double **));
                for (i = 0, j = 0; i < len[3 - 1]; i++, j += dim[ndim - 3 + 1])
                {
                    ARRAY3(a)[i] = ARRAY2(a) + j;
                }
                break;
            case 4:
                ARRAY4(a) =
                    (double ****)R_alloc(len[4 - 1], sizeof(double ***));
                for (i = 0, j = 0; i < len[4 - 1]; i++, j += dim[ndim - 4 + 1])
                {
                    ARRAY4(a)[i] = ARRAY3(a) + j;
                }
                break;
            default:
                break;
        }
    }

    for (d = 0; d < ndim; d++) DIM(a)[d] = dim[d];

    DIM_LENGTH(a) = ndim;

    return a;
}

carray make_zero_array(int *dim, int ndim)
{
    int     n, *dim2;
    long    len;
    double *vec;

    len  = 1;
    n    = ndim;
    dim2 = dim;

    while (n--) len *= (*dim2++);

    vec = (double *)R_alloc(len, sizeof(double));

    for (n = 0; n < len; n++) *vec++ = 0.0;

    return make_array(vec - len, dim, ndim);
}

carray make_matrix(double *vec, int nrow, int ncol)
{
    int dim[2];

    dim[0] = nrow;
    dim[1] = ncol;
    return make_array(vec, dim, 2);
}

carray make_zero_matrix(int nrow, int ncol)
{
    int dim[2] = {nrow, ncol};

    return make_zero_array(dim, 2);
}

carray make_vec(double *vec, const int len)
{
    int dim[] = {len};

    return make_array(vec, dim, 1);
}

carray make_zero_vec(const int len)
{
    int dim[] = {len};

    return make_zero_array(dim, 1);
}

void set_array_to_zero(carray arr)
{
    int     n;
    double *arr_ptr;

    n       = vector_length(arr);
    arr_ptr = VECTOR(arr);

    while (n--) *arr_ptr++ = 0.0;
}

carray make_identity_matrix(int n)
{
    int     i;
    carray  a;
    double *a_ptr;

    a     = make_zero_matrix(n, n);
    a_ptr = VECTOR(a);

    for (i = 0; i < n; i++)
    {
        *a_ptr  = 1.0;
        a_ptr  += (n + 1);
    }

    return a;
}

carray subarray(carray a, int index)
/* Return subarray of array a in the form of an carray
structure so it can be manipulated by other functions
NB The data are not copied, so any changes made to the
subarray will affect the original array.
*/
{
    int    i;
    long   offset;
    carray b;

    b = init_array();

    /* is index in range? */
    if (index < 0 || index >= *DIM(a)) error("Index out of range in subarray");

    offset = index;
    switch (DIM_LENGTH(a))
    {
            /* NB Falling through here */
        case 4:
            offset    *= DIM(a)[DIM_LENGTH(a) - 4 + 1];
            ARRAY3(b)  = ARRAY3(a) + offset;
        case 3:
            offset    *= DIM(a)[DIM_LENGTH(a) - 3 + 1];
            ARRAY2(b)  = ARRAY2(a) + offset;
        case 2:
            offset    *= DIM(a)[DIM_LENGTH(a) - 2 + 1];
            ARRAY1(b)  = ARRAY1(a) + offset;
            break;
        default:
            break;
    }

    DIM_LENGTH(b) = DIM_LENGTH(a) - 1;

    for (i = 0; i < DIM_LENGTH(b); i++) DIM(b)[i] = DIM(a)[i + 1];

    return b;
}

double *column(carray a, int j)

/*******************************************************************
 *
 *  Description: Returns a _copy_ of the jth column of a.
 *
 *******************************************************************/

{
    int     m;
    double *ans;
    double *a_ptr, *ans_ptr;

    if (2 != DIM_LENGTH(a) || j < 0 || j >= NCOL(a))
        error("Bad arguments for carray column function");

    m       = NROW(a);
    ans     = (double *)R_alloc(m, sizeof(double));
    ans_ptr = ans;
    a_ptr   = ARRAY1(a) + j;

    while (m--)
    {
        *ans_ptr++  = *a_ptr;
        a_ptr      += NCOL(a);
    }

    return ans;
}

int test_array_conform(const carray a1, const carray a2)
{
    int  d, ans;
    int *dim1, *dim2;

    ans = DIM_LENGTH(a1) == DIM_LENGTH(a2);

    if (ans)
    {
        d    = DIM_LENGTH(a1);
        dim1 = (int *)DIM(a1);
        dim2 = (int *)DIM(a2);

        while (ans && d--) ans = *dim1++ == *dim2++;
    }

    return ans;
}

int is_matrix(const carray a) { return 2 == DIM_LENGTH(a); }

void copy_array(carray orig, carray dest)

/*******************************************************************
 *
 *  Description: Copies orig to dest.
 *
 *******************************************************************/

{
    int       t;
    double   *o = ARRAY1(orig);
    double   *d = ARRAY1(dest);
    const int T = vector_length(dest);

    if (!test_array_conform(orig, dest))
        error("Non-conformable arrays in copy_array");

    for (t = 0; t < T; t++, d++, o++) *d = *o;
}

void transpose_matrix(carray mat, carray ans)
{
    int    i, j;
    char  *vmax;
    carray tmp;

    tmp = init_array(); /* Use swap matrix in case mat == ans */

    if (2 != DIM_LENGTH(mat) || 2 != DIM_LENGTH(ans))
        error("Non-matrices in transpose_matrix");
    if (NCOL(mat) != NROW(ans) || NROW(mat) != NCOL(ans))
        error("Non-conformable matrices in transpose_matrix");

    vmax = vmaxget();

    tmp = make_zero_matrix(NROW(ans), NCOL(ans));

    for (i = 0; i < NROW(ans); i++)
        for (j = 0; j < NCOL(ans); j++) ARRAY2(tmp)[i][j] = ARRAY2(mat)[j][i];

    copy_array(tmp, ans);

    vmaxset(vmax);
}

void transpose_matrix_ns(carray mat, carray ans)

/*******************************************************************
 *
 *  Description: The "_ns" suffix stands for "no swap".
 *    Use if you are sure that mat != ans.
 *
 *******************************************************************/

{
    int i, j;

    if (ARRAY1(mat) == ARRAY1(ans))
        error("ARRAY1(mat) == ARRAY1(ans) tranpose_matrix_ns");

    if (2 != DIM_LENGTH(mat) || 2 != DIM_LENGTH(ans))
        error("Non-matrices in transpose_matrix");
    if (NCOL(mat) != NROW(ans) || NROW(mat) != NCOL(ans))
        error("Non-conformable matrices in transpose_matrix");

    for (i = 0; i < NROW(ans); i++)
        for (j = 0; j < NCOL(ans); j++) ARRAY2(ans)[i][j] = ARRAY2(mat)[j][i];
}

carray transpose_matrix_b(carray mat)
{
    carray ans;

    ans = make_zero_matrix(NCOL(mat), NROW(mat));

    transpose_matrix_ns(mat, ans);

    return ans;
}

void transpose_matrix_ow(carray a)

/********************************************************************
 *
 *   Description: Replaces a by its transpose.
 *
 ********************************************************************/

{
    carray b;

    b = make_zero_matrix(NCOL(a), NROW(a));

    transpose_matrix_ns(a, b);

    a = b;
}

void array_op(carray arr1, carray arr2, char op, carray ans)
/* Element-wise array operations */
{
    int     n;
    double *ans_ptr, *ptr1, *ptr2;

    if (!test_array_conform(arr1, arr2) || !test_array_conform(arr2, ans))
        error("Non-conformable matrices in array_op");

    n       = vector_length(ans);
    ans_ptr = VECTOR(ans);
    ptr1    = VECTOR(arr1);
    ptr2    = VECTOR(arr2);

    switch (op)
    {
        case '*':
            while (n--) *ans_ptr++ = (*ptr1++) * (*ptr2++);
            break;
        case '+':
            while (n--) *ans_ptr++ = (*ptr1++) + (*ptr2++);
            break;
        case '/':
            while (n--) *ans_ptr++ = (*ptr1++) / (*ptr2++);
            break;
        case '-':
            while (n--) *ans_ptr++ = (*ptr1++) - (*ptr2++);
            break;
        default:
            printf("Unknown op in array_op");
    }
}

carray array_op_b(carray arr1, carray arr2, char op)
{
    carray ans;

    ans = make_zero_matrix(NROW(arr1), NCOL(arr1));

    array_op(arr1, arr2, op, ans);

    return ans;
}

void scalar_op(carray arr, double s, char op, carray ans)
/* Elementwise scalar operations */
{
    int     n;
    double *arr_ptr, *ans_ptr;

    assert2(test_array_conform(arr, ans));

    n       = vector_length(ans);
    ans_ptr = VECTOR(ans);
    arr_ptr = VECTOR(arr);

    switch (op)
    {
        case '*':
            while (n--) *ans_ptr++ = (*arr_ptr++) * s;
            break;
        case '+':
            while (n--) *ans_ptr++ = (*arr_ptr++) + s;
            break;
        case '/':
            while (n--) *ans_ptr++ = (*arr_ptr++) / s;
            break;
        case '-':
            while (n--) *ans_ptr++ = (*arr_ptr++) - s;
            break;
        default:
            printf("Unknown op in array_op");
    }
}

carray scalar_op_b(carray arr, double s, char op)
{
    carray ans;
    ans = make_zero_matrix(NROW(arr), NCOL(arr));

    scalar_op(arr, s, op, ans);

    return ans;
}

void matrix_prod(carray mat1, carray mat2, int trans1, int trans2, carray ans)

/*******************************************************************
 *
 *  Description:
 *    General matrix product between mat1 and mat2.
 *    Puts answer in ans.
 *    Parameters trans1 and trans2 are logical flags which
 *    indicate whether mat1 or mat2 respectively are
 *    to be transposed.
 *    Normal matrix multiplication has trans1 = trans2 = 0.
 *
 *******************************************************************/

{
    int             i, j, k, K1, K2, ok;
    char           *vmax;
    register double m1, m2;
    carray          tmp;

    /* Test whether everything is a matrix */
    if (!is_matrix(mat1) || !is_matrix(mat2) || !is_matrix(ans))
        error("Non-matrix passed to matrix_prod");

    /* Test whether matrices conform. K is the dimension that is
       lost by multiplication */

    ok = trans1 ? (K1 = NROW(mat1), NCOL(mat1) == NROW(ans))
                : (K1 = NCOL(mat1), NROW(mat1) == NROW(ans));

    ok *= trans2 ? (K2 = NCOL(mat2), NROW(mat2) == NCOL(ans))
                 : (K2 = NROW(mat2), NCOL(mat2) == NCOL(ans));

    ok *= K1 == K2;

    if (!ok) error("Non-conformable matrices in matrix_prod");

    /*
    **    In case ans is the same as mat1 or mat2, we create a temporary
    **    matrix tmp to hold the answer, then copy it to ans
    */

    tmp = init_array();

    vmax = vmaxget();

    tmp = make_zero_matrix(NROW(ans), NCOL(ans));

    for (i = 0; i < NROW(tmp); i++)
    {
        for (j = 0; j < NCOL(tmp); j++)
        {
            for (k = 0; k < K1; k++)
            {
                m1 = trans1 ? MATRIX(mat1)[k][i] : MATRIX(mat1)[i][k];
                m2 = trans2 ? MATRIX(mat2)[j][k] : MATRIX(mat2)[k][j];

                MATRIX(tmp)[i][j] += m1 * m2;
            }
        }
    }

    copy_array(tmp, ans);

    vmaxset(vmax);
}

void matrix_prod_ns(carray mat1, carray mat2, int trans1, int trans2,
                    carray ans)

/*******************************************************************
 *
 *  Description:
 *    General matrix product between mat1 and mat2. Put answer in ans.
 *    Only use if you know that ans != mat1 and ans != mat2.
 *    trans1 and trans2 are logical flags which indicate if the matrix is
 *    to be transposed.
 *    Normal matrix multiplication has trans1 = trans2 = 0.
 *
 *******************************************************************/

{
    int             i, j, k, K1, K2, ok;
    register double m1, m2;
    carray          tmp;

    /* Test whether ans == mat1 or ans == mat2 */
    if (ARRAY1(ans) == ARRAY1(mat1) || ARRAY1(ans) == ARRAY1(mat2))
        error("ARRAY1(ans) == ARRAY1(mat1) || ARRAY1(ans) == ARRAY1(mat2)");

    /* Test whether everything is a matrix */
    if (!is_matrix(mat1) || !is_matrix(mat2) || !is_matrix(ans))
        error("Non-matrix passed to matrix_prod");

    /* Test whether matrices conform. K is the dimension that is
       lost by multiplication */

    ok = trans1 ? (K1 = NROW(mat1), NCOL(mat1) == NROW(ans))
                : (K1 = NCOL(mat1), NROW(mat1) == NROW(ans));

    ok *= trans2 ? (K2 = NCOL(mat2), NROW(mat2) == NCOL(ans))
                 : (K2 = NROW(mat2), NCOL(mat2) == NCOL(ans));

    ok *= (K1 == K2);

    if (!ok) error("Non-conformable matrices in matrix_prod");

    tmp = ans; /* Here tmp is just basically an alias for ans */

    for (i = 0; i < NROW(tmp); ++i)
        for (j = 0; j < NCOL(tmp); ++j)
        {
            MATRIX(tmp)[i][j] = 0.0;

            for (k = 0; k < K1; ++k)
            {
                m1 = trans1 ? MATRIX(mat1)[k][i] : MATRIX(mat1)[i][k];

                m2 = trans2 ? MATRIX(mat2)[j][k] : MATRIX(mat2)[k][j];

                MATRIX(tmp)[i][j] += m1 * m2;
            }
        }
}

carray matrix_prod_b(carray mat1, carray mat2, int trans1, int trans2)
{
    carray ans;
    int    m, n;

    m = trans1 ? NCOL(mat1) : NROW(mat1);
    n = trans2 ? NROW(mat2) : NCOL(mat2);

    ans = make_zero_matrix(m, n);

    matrix_prod_ns(mat1, mat2, trans1, trans2, ans);

    return ans;
}

void matrix_prod_3(carray mat1, carray mat2, carray mat3, int trans1,
                   int trans2, int trans3, carray ans)
{
    int    dim1, dim2;
    carray a;

    dim1 = trans1 ? NCOL(mat1) : NROW(mat1);
    dim2 = trans2 ? NROW(mat2) : NCOL(mat2);

    a = make_zero_matrix(dim1, dim2);

    matrix_prod_ns(mat1, mat2, trans1, trans2, a);

    matrix_prod(a, mat3, 0, trans3, ans);
}

carray matrix_prod_3b(carray mat1, carray mat2, carray mat3, int trans1,
                      int trans2, int trans3)
{
    int    dim1, dim2;
    carray a;

    dim1 = (trans1 ? NCOL(mat1) : NROW(mat1));
    dim2 = (trans3 ? NROW(mat3) : NCOL(mat3));

    a = make_zero_matrix(dim1, dim2);

    matrix_prod_3(mat1, mat2, mat3, trans1, trans2, trans3, a);

    return a;
}

void matrix_prod_4(carray mat1, carray mat2, carray mat3, carray mat4,
                   int trans1, int trans2, int trans3, int trans4, carray ans)
{
    int    dim1, dim2, dim3, dim4;
    carray tmp1, tmp2;

    dim1 = trans1 ? NCOL(mat1) : NROW(mat1);
    dim2 = trans2 ? NROW(mat2) : NCOL(mat2);
    dim3 = trans3 ? NCOL(mat3) : NROW(mat3);
    dim4 = trans4 ? NROW(mat4) : NCOL(mat4);

    tmp1 = make_zero_matrix(dim1, dim2);
    tmp2 = make_zero_matrix(dim3, dim4);

    matrix_prod(mat1, mat2, trans1, trans2, tmp1);
    matrix_prod(mat3, mat4, trans3, trans4, tmp2);

    matrix_prod(tmp1, tmp2, 0, 0, ans);
}

void matrix_prod_n(carray *mat, int *trans, int n, carray ans)
{
    carray tmp1, tmp2;

    if (n <= 0)
        error("Bad value for number of matrices in matrix_prod_n");
    else if (1 == n)
        if (*trans)
            transpose_matrix(*mat, ans);
        else
            copy_array(*mat, ans);
    else
    {
        tmp1 = make_zero_matrix(*trans ? NCOL(*mat) : NROW(*mat),
                                *trans ? NROW(*mat) : NCOL(*mat));

        if (*trans)
            transpose_matrix(*mat++, tmp1);
        else
            copy_array(*mat++, tmp1);
        ++trans;

        while (--n > 0)
        {
            tmp2 =
                make_zero_matrix(NROW(tmp1), *trans ? NROW(*mat) : NCOL(*mat));
            matrix_prod(tmp1, *mat++, 0, *trans++, tmp2);
            tmp1 = tmp2;
        }

        copy_array(tmp1, ans);
    }
}

void matrix_ptr_prod(carray *mat, int *trans, int n, carray ans)
{
    int    i;
    carray tmp1, tmp2;

    if (n <= 0)
        error("Bad value for number of matrices in matrix_ptr_prod");
    else if (1 == n)
    {
        if (*trans)
            transpose_matrix(*mat, ans);
        else
            copy_array(*mat, ans);

    } else
    {
        tmp1 = make_zero_matrix(*trans ? NCOL(*mat) : NROW(*mat),
                                *trans ? NROW(*mat) : NCOL(*mat));

        if (*trans++)
            transpose_matrix(*mat++, tmp1);
        else
            copy_array(*mat++, tmp1);

        for (i = 1; i < n; i++, mat++, trans++)
        {
            tmp2 =
                make_zero_matrix(NROW(tmp1), *trans ? NROW(*mat) : NCOL(*mat));

            matrix_prod_ns(tmp1, *mat, 0, *trans, tmp2);

            tmp1 = tmp2;
        }

        copy_array(tmp1, ans);
    }
}

int matrix_inverse(carray x, carray xi)

/*******************************************************************
 *
 *  Description: Attempts to invert matrix a, putting the inverse
 *    into ai.  Returns 1 if inversion succeeds and 0 if it fails.
 *
 *******************************************************************/

{
    carray b;
    int    n;

    n = NROW(x);

    if (n != NCOL(x) || 2 != DIM_LENGTH(x))
        error("Attempt to invert non-matrix or non-square matrix");

    if (!test_array_conform(x, xi))
        error("Non-conformable arrays in matrix_inverse");

    copy_array(x, xi);
    b = make_identity_matrix(n);

    return gaussj2(MATRIX(xi), n, MATRIX(b), n);
}

carray matrix_inverse_b(carray a, int *invertible)
{
    carray ai;

    ai = make_zero_matrix(NROW(a), NCOL(a));

    *invertible = matrix_inverse(a, ai);

    return ai;
}

int matrix_inverse_2x2(carray a, carray ai)

/*******************************************************************
 *
 *  Description: Basically just like matrix_inverse but for 2x2
 *    matrices only.  Returns 0 iff matrix a is not invertible.
 *
 *******************************************************************/

{
    double d, tmp[4] = {0};

    if (2 != DIM_LENGTH(a) || 2 != DIM_LENGTH(ai) || 2 != NROW(a) ||
        2 != NCOL(a) || 2 != NROW(ai) || 2 != NCOL(ai))
        error("Bad arguments passed to matrix_inverse_2x2");

    tmp[0] = ARRAY1(a)[0];
    tmp[1] = ARRAY1(a)[1];
    tmp[2] = ARRAY1(a)[2];
    tmp[3] = ARRAY1(a)[3];

    d = tmp[0] * tmp[3] - tmp[1] * tmp[2];

    ARRAY1(ai)[0] = tmp[3] / d;
    ARRAY1(ai)[1] = -tmp[1] / d;
    ARRAY1(ai)[2] = -tmp[2] / d;
    ARRAY1(ai)[3] = tmp[0] / d;

    return (0.0 != d && R_FINITE(d));
}

carray vec(carray a)
/* Returns a vector_length(a) x 1 matrix whose elements
are the columns of a stacked on top of each other. */
{
    return make_matrix(ARRAY1(transpose_matrix_b(a)), vector_length(a), 1);
}

void kronecker(carray a, carray b, carray ans)
{
    int i, j, k, h;
    int m, n, p, q;
    /* double *a_i, *b_i, *ans_i; */
    /* double **a_c, **b_c, **ans_c; */

    if (2 != DIM_LENGTH(a) || 2 != DIM_LENGTH(b))
        error("Non-matrices passed to kronecker");

    if (NROW(a) * NROW(b) != NROW(ans) || NCOL(a) * NCOL(b) != NCOL(ans))
        error("Non-conformable matrices in kronecker");

    m = NROW(a);
    n = NCOL(a);
    p = NROW(b);
    q = NCOL(b);

    for (i = 0; i < n; i++)
        for (j = 0; j < q; j++)
            for (k = 0; k < m; k++)
                for (h = 0; h < p; h++) ARRAY2(ans)
    [h + p * k][j + q * i] = ARRAY2(a)[k][i] * ARRAY2(b)[h][j];
}

double sumsq(carray a)
/* Returns the sum of squares of entries of a */
{
    register double ans  = 0.0;
    double         *aPtr = ARRAY1(a);
    double const   *end  = aPtr + vector_length(a);

    for (; aPtr < end; aPtr++) ans += SQR(*aPtr);

    return ans;
}

double det(carray a)

/*******************************************************************
 *
 *  Description: Returns the determinant the matrix argument.
 *
 *******************************************************************/

{
    int    j, n = NROW(a);
    int   *ix = NULL;
    double d  = 0.0;
    carray b;

    if (2 != DIM_LENGTH(a) || n != NCOL(a))
        error("C error: det.  Expected square matrix argument.");

    ix = (int *)malloc(n * sizeof(int));

    if (!ix) error("C error: det.  Malloc failed to allocate memory.");

    b = make_zero_matrix(n, n);
    copy_array(a, b);

    if (ludcmp2(ARRAY2(b), n, ix, &d))
    {
        for (j = 0; j < n; j++) d *= ARRAY2(b)[j][j];
    } else
        d = 0.0;

    free(ix);

    return d;
}

carray sexp_to_carray(SEXP A, int dup)

/********************************************************************
 *
 *   Description: Converts an R matrix to a Carray and returns
 *     result.  Note that the dimensions have to be switched as
 *     R matrix representation is column-major whereas carray
 *     matrices are row-major.
 *
 *   Caution: Use only dup = 1 for now.
 *
 ********************************************************************/

{
    int    ndim, dims[1];
    SEXP   Adim, A1;
    carray ans, tmp;

    dup = 1; /* Temporary until we debug the dup = 0 case */

    PROTECT(Adim = getAttrib(A, R_DimSymbol));

    PROTECT(A1 = dup ? AS_NUMERIC(duplicate(A)) : A);
    setAttrib(A1, R_DimSymbol, R_NilValue);

    ndim    = isNull(Adim) ? 1 : length(Adim);
    dims[0] = length(A1);

    switch (ndim)
    {
        case 1:

            ans = make_array(REAL(A1), dims, 1);
            break;

        case 2:

            /* Note reversal of dimensions: */
            tmp = make_matrix(REAL(A1), INTEGER(Adim)[1], INTEGER(Adim)[0]);
            ans = make_zero_matrix(NCOL(tmp), NROW(tmp));
            transpose_matrix(tmp, ans);
            break;

        default:

            error("SEXP has to many dimensions");
            break;
    }

    UNPROTECT(2);

    return ans;
}

SEXP carray_to_sexp(carray A)
/* As remarked in the previous function the dimensions need
   to be reversed in the transformation. */
{
    const int ndim = DIM_LENGTH(A);
    carray    B;
    SEXP      ans = R_NilValue, dim = R_NilValue;

    if (2 < DIM_LENGTH(A)) error("Too many dimensions in carray_to_sexp");

    switch (ndim)
    {
        case 1:

            PROTECT(ans = NEW_NUMERIC(vector_length(A)));
            memcpy(REAL(ans), ARRAY1(A), vector_length(A) * sizeof(double));
            break;

        case 2:

            B = make_zero_matrix(NCOL(A), NROW(A));
            transpose_matrix(A, B);
            PROTECT(ans = NEW_NUMERIC(vector_length(A)));
            PROTECT(dim = NEW_INTEGER(2));
            memcpy(INTEGER(dim), DIM(A), 2 * sizeof(int));
            memcpy(REAL(ans), ARRAY1(B), vector_length(A) * sizeof(double));
            setAttrib(ans, R_DimSymbol, dim);
            break;

        default:
            break;
    }

    UNPROTECT(ndim);

    return ans;
}

#undef SWAP
#undef NRANSI
#undef TINY
#undef USE_MALLOC
