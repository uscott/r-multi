

#ifdef NOT_NOW

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#ifndef URI_H_
#include "multi.h"
#endif

SEXP test1_10(SEXP x)
{
    carray y;

    y = make_zero_matrix(ncols(x), nrows(x));

    transpose_matrix_ns(sexp_to_carray(x, 1), y);

    return carray_to_sexp(y);
}

SEXP test1_9(SEXP x)
{
    double vec[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int    dim[2] = {3, 2};
    int    ndim   = 2;

    return carray_to_sexp(make_array(vec, dim, ndim));
}

SEXP test1_8(SEXP x) { return carray_to_sexp(make_zero_matrix(3, 3)); }

SEXP test1_7(SEXP x) { return (AS_NUMERIC(x)); }

SEXP test1_6(SEXP x)
{
    int   i;
    char *names[] = {"xinverse", "invertible"};
    SEXP  ans;

    PROTECT(ans = NEW_LIST(2));
    set_names(ans, names);

    SET_ELT(ans, 0, carray_to_sexp(matrix_inverse_b(sexp_to_carray(x, 1), &i)));
    SET_ELT(ans, 1, ScalarInteger(i));

    UNPROTECT(1);

    return ans;
}

SEXP test1_5(SEXP x)
{
    return ScalarReal(det(sexp_to_carray(x, /* dup = */ 1)));
}

SEXP test1_4(SEXP x)
{
    return carray_to_sexp(
        matrix_prod_b(sexp_to_carray(x, 0), sexp_to_carray(x, 0), 1, 1));
}

SEXP test1_3(SEXP x)
{
    carray y;

    if (nrows(x) != ncols(x)) error("!!");

    y = make_zero_matrix(nrows(x), ncols(x));

    matrix_prod(sexp_to_carray(x, 1), sexp_to_carray(x, 1), 0, 1, y);

    return carray_to_sexp(y);
}

SEXP test1_2(SEXP x)
{
    carray X;
    X = make_zero_matrix(nrows(x), ncols(x));

    matrix_inverse(sexp_to_carray(x, 1), X);

    return carray_to_sexp(X);
}

SEXP test1_1(SEXP x)
{
    return carray_to_sexp(transpose_matrix_b(sexp_to_carray(x, 1)));
}

SEXP test1_0(SEXP x)
{
    carray y;

    y = make_zero_matrix(nrows(x), ncols(x));

    copy_array(sexp_to_carray(x, 1), y);

    return carray_to_sexp(y);
}

void test2_8(double *x, int *n)
{
    int i;

    for (i = 0; i < *n; i++) x[i] = SQR(x[i]);
}

SEXP test2_6(SEXP x, SEXP i)
{
    SEXP y;
    int  n;

    if (!isMatrix(x)) error("asdl;fja");

    n = ncols(x);

    y = PROTECT(NEW_NUMERIC(n));

    dmatrow1(x, asInteger(i), REAL(y));

    UNPROTECT(1);

    return y;
}

SEXP test2_5(SEXP x, SEXP i)
{
    char     *a;
    const int n = ncols(x);
    SEXP      ans;

    a = (char *)R_alloc(n, sizeof(char));

    cmatrow1(x, asInteger(i), a);

    PROTECT(ans = NEW_STRING(1));

    STRING_PTR(ans)[0] = mkString(a);

    UNPROTECT(1);

    return ans;
}

SEXP test2_4(SEXP x, SEXP y)

/*******************************************************************
 *
 *  Description: Tests Kronecker product code in carray.c
 *
 *******************************************************************/

{
    carray z;

    z = make_zero_matrix(nrows(x) * nrows(y), ncols(x) * ncols(y));

    kronecker(sexp_to_carray(x, 1), sexp_to_carray(y, 1), z);

    return carray_to_sexp(z);
}

SEXP test2_3(SEXP x, SEXP i)
{
    return ScalarInteger('a' == *CHAR(STRING_ELT(x, asInteger(i))));
}

SEXP test2_2(SEXP x, SEXP t)
{
    int      i, n;
    carray  *mat1, *mat2, a, b;
    carray **mat1_1, **mat2_1;
    SEXP     ans;

    if (!isNewList(x) || length(x) != length(t)) error("!");

    n = length(x);
    PROTECT(t = AS_INTEGER(t));

    mat2 = mat1 = (carray *)R_alloc(n, sizeof(carray));
    mat2_1 = mat1_1 = (carray **)R_alloc(n, sizeof(carray *));

    for (i = 0; i < n; i++)
    {
        *mat2     = sexp_to_carray(AS_NUMERIC(GET_ELT(x, i)), 1);
        *mat2_1++ = mat2++;
    }

    a = make_zero_matrix(NROW(*mat1), NCOL(mat1[n - 1]));
    b = make_zero_matrix(NROW(*mat1), NCOL(mat1[n - 1]));

    matrix_prod_n(mat1, INTEGER(t), n, a);
    matrix_ptr_prod(mat1_1, INTEGER(t), n, b);

    PROTECT(ans = NEW_LIST(2));
    SET_ELT(ans, 0, carray_to_sexp(a));
    SET_ELT(ans, 1, carray_to_sexp(b));

    UNPROTECT(2);

    return ans;
}

SEXP test2_(SEXP x, SEXP y)
{
    carray x1, y1, z1;

    if (!isMatrix(x) || !isMatrix(y)) error("!");

    x1 = sexp_to_carray(x, 1);
    y1 = sexp_to_carray(y, 1);

    z1 = make_zero_matrix(NROW(x1), NCOL(y1));
    matrix_prod(x1, y1, 0, 0, z1);

    return carray_to_sexp(z1);
}

SEXP test3_1(SEXP x, SEXP i, SEXP j)
{
    int m;

    if (!isMatrix(x)) error("fuck off");

    m = nrows(x);

    return mkString(CHAR(STRING_ELT(x, asInteger(i) + m * asInteger(j))));
}

#endif
