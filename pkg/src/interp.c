#include "multi.h"

void interp_crv(long len, double *x_in, double *y_in, double x0, double *y_out)
{
    const double TOL = 1e-10;
    double max_sprd  = -1, tmp_sprd, sprd_lower, sprd_upper, x1, x2, y1, y2, w1,
           w2;
    int  found_lower = 0, found_upper = 0;
    long i;

    /* Get max sprd from x0 */
    for (i = 0; i < len; ++i)
        if ((tmp_sprd = ABS(x_in[i] - x0)) > max_sprd) max_sprd = tmp_sprd;

    sprd_lower = sprd_upper = max_sprd + 1;

    for (i = 0; i < len; ++i)
    {
        tmp_sprd = ABS(x_in[i] - x0);

        if (x_in[i] < x0 && tmp_sprd < sprd_lower)
        {
            sprd_lower  = tmp_sprd;
            x1          = x_in[i];
            y1          = y_in[i];
            found_lower = 1;
        }

        if (x_in[i] >= x0 && tmp_sprd < sprd_upper)
        {
            sprd_upper  = tmp_sprd;
            x2          = x_in[i];
            y2          = y_in[i];
            found_upper = 1;
        }
    }

    if (!found_lower || ABS(x2 - x0) < TOL)
        *y_out = y2;
    else if (!found_upper)
        *y_out = y1;
    else
    {
        w1 = (x2 - x0) / (x2 - x1);
        w2 = (x0 - x1) / (x2 - x1);

        *y_out = w1 * y1 + w2 * y2;
    }
}

SEXP R_interp_crv(SEXP x_in, SEXP y_in, SEXP x0)
{
    long x0_len = length(x0), i;
    SEXP y_out  = NEW_NUMERIC(x0_len);
    PROTECT(y_out);

    for (i = 0; i < x0_len; ++i)
        interp_crv(length(x_in), REAL(x_in), REAL(y_in), REAL(x0)[i],
                   REAL(y_out) + i);

    UNPROTECT(1);

    return y_out;
}

void interp_sfc(long len, double *t, double *m, double *v, double t0, double m0,
                double *v_out, double *m_buf, double *v_buf)
{
    const double TOL     = 1e-10;
    const int    BUF_LEN = 64;
    double sprd_max = -1, tmp, sprd_lower = -1, sprd_upper = -1, t1, t2, v1, v2,
           w1, w2, fwd_var, buf1[BUF_LEN], buf2[BUF_LEN];
    int  found_lower = 0, found_upper = 0, used_malloc_m = 0, used_malloc_v = 0;
    long i, len_loc;

    /* Get max sprd from t0 */
    for (i = 0; i < len; ++i)
        if ((tmp = ABS(t[i] - t0)) > sprd_max) sprd_max = tmp;

    sprd_lower = sprd_upper = sprd_max + 1;

    for (i = 0; i < len; ++i)
    {
        tmp = ABS(t[i] - t0);

        if (t[i] < t0 && tmp < sprd_lower)
        {
            t1          = t[i];
            sprd_lower  = tmp;
            found_lower = 1;
        }

        if (t[i] >= t0 && tmp < sprd_upper)
        {
            t2          = t[i];
            sprd_upper  = tmp;
            found_upper = 1;
        }
    }

    if (!m_buf)
    {
        used_malloc_m = len > BUF_LEN;
        m_buf =
            (double *)(len <= BUF_LEN ? buf1 : malloc(len * sizeof(double)));
    }

    if (!v_buf)
    {
        used_malloc_v = len > BUF_LEN;
        v_buf =
            (double *)(len <= BUF_LEN ? buf2 : malloc(len * sizeof(double)));
    }

    if (found_lower)
    {
        len_loc = 0;

        for (i = 0; i < len; ++i)
            if (ABS(t[i] - t1) < TOL)
            {
                m_buf[len_loc] = m[i];
                v_buf[len_loc] = v[i];
                ++len_loc;
            }

        interp_crv(len_loc, m_buf, v_buf, m0, &v1);
    }

    if (found_upper)
    {
        len_loc = 0;

        for (i = 0; i < len; ++i)
            if (ABS(t[i] - t2) < TOL)
            {
                m_buf[len_loc] = m[i];
                v_buf[len_loc] = v[i];
                ++len_loc;
            }

        interp_crv(len_loc, m_buf, v_buf, m0, &v2);
    }

    if (!found_lower || ABS(t2 - t0) < TOL)
        *v_out = v2;
    else if (!found_upper)
        *v_out = v1;
    else
    {
        w1 = (t2 - t0) / (t2 - t1);
        w2 = (t0 - t1) / (t2 - t1);

        *v_out = w1 * v1 + w2 * v2;
    }

    if (used_malloc_m) free(m_buf);

    if (used_malloc_v) free(v_buf);
}

SEXP R_interp_sfc(SEXP t_in, SEXP m_in, SEXP v_in, SEXP t_arg, SEXP m_arg)
{
    long    arg_len = length(m_arg), len_in = length(m_in), i;
    double *m_buf, *v_buf;
    SEXP    v_out = NEW_NUMERIC(arg_len);
    PROTECT(v_out);

    m_buf = (double *)malloc(len_in * sizeof(double));
    v_buf = (double *)malloc(len_in * sizeof(double));

    for (i = 0; i < arg_len; ++i)
        interp_sfc(len_in, REAL(t_in), REAL(m_in), REAL(v_in), REAL(t_arg)[i],
                   REAL(m_arg)[i], REAL(v_out) + i, m_buf, v_buf);

    free(m_buf);
    free(v_buf);

    UNPROTECT(1);

    return v_out;
}
