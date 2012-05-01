#ifdef KISS_MY_ASS
#include "uri.h"

void get_init_params(int T, double *alpha, double *delta,
		     double *sigv, carray m, carray M, double *N,
		     double *lambda, double *h)
{

  int t;
  double *Lh = h, *ptr, rho, sig1, sig2, e1, e2;
  const double *h0 = h;
  carray y, X, Xt, Mi, res, b;

  GetRNGstate();

  /* Draw alpha, delta, sigv from priors */
  *N = 4.0;
  *lambda  = 1.0 / 0.005;
  *sigv = sqrt(1.0 / rgamma(/* shape = */ 0.5 * (*N),
			       /* scale = */ 2.0 / (*lambda)));

  *delta = 2.0 * unif_rand() - 1.0;
  *alpha = norm_rand();

  /* Simulate h[t] */
  *h++ = rnorm(/* mu = */ (*alpha) / (1.0 - (*delta)),
	       /* sigma = */ *sigv / sqrt(1.0 - (*delta) * (*delta)));

  for (t = 1; t < T; t++)
    *h++ = (*alpha) + (*delta) * (*Lh++) + (*sigv) * rnorm(0.0, 1.0);

  y         = make_zero_matrix(T - 1, 1);
  res       = make_zero_matrix(T - 1, 1);
  Xt        = make_zero_matrix(2, T - 1);
  X         = make_zero_matrix(T - 1, 2);
  Mi = make_zero_matrix(2, 2);
  b         = make_zero_matrix(2, 1);

  /*  Calculate X, Xt, y */
  ptr = *ARRAY2(Xt);
  for (t = 1; t < T; t++)
    *ptr++ = 1.0;

  memcpy(ARRAY2(Xt)[1], h0, (T - 1) * sizeof(double));
  memcpy(ARRAY1(y), h0 + 1, (T - 1) * sizeof(double));

  transpose_matrix(Xt, X);
  matrix_prod(Xt, X, 0, 0, Mi);
  matrix_inverse_2x2(Mi, M);
  matrix_prod_3(M, Xt, y, 0, 0, 0, b);
  matrix_prod_3(M, Xt, y, 0, 0, 0, m);
  
  /* Update N and lambda */
  *N = *N + T - 1;
  array_op(y, matrix_prod_b(X, b, 0, 0), '-', res);

  *lambda += sumsq(res);
  
  /* Draw new sigv */
  
  *sigv = sqrt(1.0 / rgamma(0.5 * (*N), 2.0 / (*lambda)));
 

  /* Draw new alpha, delta conditional on new sigma.v */

  sig1 = (*sigv) * sqrt(ARRAY1(M)[0]);
  sig2 = (*sigv) * sqrt(ARRAY1(M)[3]);
  rho  = (*sigv) * (*sigv) * ARRAY1(M)[1] / (sig1 * sig2);

  e1 = rnorm(0, 1);
  e2 = rho * e1 + sqrt(1 - rho * rho) * rnorm(0, 1);
  (*alpha) = ARRAY1(m)[0] + sig1 * e1;
  (*delta) = ARRAY1(m)[1] + sig2 * e2;

  /* Continue until |delta| < 1 for stationarity */
    while (fabs(*delta) >= 1) {

    e1 = rnorm(0, 1);
    e2 = rho * e1 + sqrt(1 - rho * rho) * rnorm(0, 1);
    (*alpha) = ARRAY1(m)[0] + sig1 * e1;
    (*delta) = ARRAY1(m)[1] + sig2 * e2;

  }
  

  /* Transform h -> exp(h) */

    h = (double *) h0;
    for (t = 0; t < T; t++, h++)
      *h = exp(*h);
      

  PutRNGstate();
}




void update_params(int T, double *alpha, double *delta,
		   double *sigv, carray m, carray M, double *N,
		   double *lambda, double *h)
{
  int t;
  double *ptr1, *ptr2, *ptr3, rho, sig1, sig2, e1, e2;
  /* const double *h0 = h; */
  carray y, X, Xt, Mi, res, b, XtXi, tmp_2x2, tmp_1x2, tmp_2x1;

  GetRNGstate();

  y         = make_zero_matrix(T - 1, 1);
  res       = make_zero_matrix(T - 1, 1);
  Xt        = make_zero_matrix(2, T - 1);
  X         = make_zero_matrix(T - 1, 2);
  Mi = make_zero_matrix(2, 2);
  tmp_2x2   = make_zero_matrix(2, 2);
  tmp_1x2   = make_zero_matrix(1, 2);
  tmp_2x1   = make_zero_matrix(2, 1);
  b         = make_zero_matrix(2, 1);
  XtXi      = make_zero_matrix(2, 2);


  ptr1 = ARRAY2(Xt)[0];
  ptr2 = ARRAY2(Xt)[1];
  ptr3 = ARRAY1(y);
  for (t = 1; t < T; t++) {
    *ptr1++ = 1.0;
    *ptr2++ = log(*h++);
    *ptr3++ = log(*h);
  }

  transpose_matrix(Xt, X);
  matrix_inverse_2x2(matrix_prod_b(Xt, X, 0, 0), XtXi);

  matrix_prod_3(XtXi, Xt, y, 0, 0, 0, b);
  
  // Update N and lambda
  *N = (*N) + T - 1.0;

  array_op(y, matrix_prod_b(X, b, 0, 0), '-', res);  

  matrix_inverse_2x2(M, Mi);
  
  matrix_prod(Xt, X, 0, 0, tmp_2x2);
  array_op(tmp_2x2, Mi, '+', tmp_2x2);
  matrix_prod_4(Mi, tmp_2x2, Xt, X, 0, 0, 0, 0, tmp_2x2);
  array_op(b, m, '-', tmp_2x1);
  matrix_prod(tmp_2x1, tmp_2x2, 1, 0, tmp_1x2);
  *lambda += sumsq(res) + *ARRAY1(matrix_prod_b(tmp_1x2, tmp_2x1, 0, 0));

  
  /* Draw new sigv */  
  *sigv = sqrt(1.0 / rgamma(0.5 * (*N), 2.0 / (*lambda)));
  
  // Draw new alpha, delta conditional on new sigma.v
  matrix_prod(Xt, X, 0, 0, tmp_2x2);
  array_op(Mi, tmp_2x2, '+', tmp_2x2);
  matrix_inverse_2x2(tmp_2x2, M);
  matrix_prod(Mi, m, 0, 0, m);
  matrix_prod(Xt, y, 0, 0, tmp_2x1);
  ARRAY1(m)[0] += ARRAY1(tmp_2x1)[0];
  ARRAY1(m)[1] += ARRAY1(tmp_2x1)[1];
  matrix_prod(M, m, 0, 0, m);
  matrix_inverse_2x2(M, Mi);

  /* Draw new alpha, delta conditional on new sigv */

  sig1 = (*sigv) * sqrt(ARRAY1(M)[0]);
  sig2 = (*sigv) * sqrt(ARRAY1(M)[3]);
  rho  = (*sigv) * (*sigv) * ARRAY1(M)[1] / (sig1 * sig2);

  e1 = rnorm(0, 1);
  e2 = rho * e1 + sqrt(1 - rho * rho) * rnorm(0, 1);
  (*alpha) = ARRAY1(m)[0] + sig1 * e1;
  (*delta) = ARRAY1(m)[1] + sig2 * e2;

  /* Continue until |delta| < 1 for stationarity */

    while (fabs(*delta) >= 1) {

    e1 = rnorm(0, 1);
    e2 = rho * e1 + sqrt(1 - rho * rho) * rnorm(0, 1);
    (*alpha) = ARRAY1(m)[0] + sig1 * e1;
    (*delta) = ARRAY1(m)[1] + sig2 * e2;

  }

  PutRNGstate();

}



void update_h(double *alpha, double *delta, double *sigv,
	      double *h, double *r, int *T, double *h_new)
{
  int t;
  double *h_minus, *h_plus;
  double mu, sigma2, phi, phi_LN, theta, theta_1, theta_LN, 
    x_mode, c, p1, p2, q1, q2, num, denom;

  r      += 1;
  h_new  += 1;
  h_minus = ++h - 1;
  h_plus  =   h + 1;

  GetRNGstate();

  for (t = 1; t < (*T) - 1; t++, h++, h_new++, r++) {

    mu  = (*alpha) * (1.0 - (*delta)) 
      + (*delta) * (log(*h_minus++) + log(*h_plus++));
    mu /= (1.0 + (*delta) * (*delta));

    sigma2 = (*sigv) * (*sigv) * (1.0 + (*delta) * (*delta));

    phi_LN = (1.0 - 2.0 * exp(sigma2)) / (1.0 - exp(sigma2));
    phi    = -0.5 + phi_LN - 1.0;

    theta_LN = (phi_LN - 0.5) * exp(mu + 0.5 * sigma2);
    theta_1  = 0.5 * (*r) * (*r);
    theta    = theta_1 + theta_LN;

    x_mode = theta / (phi + 1.0);
    c = 1.1 * R_pow(x_mode, -1.5) * exp(-0.5 * (*r) * (*r) / x_mode)
      * exp(-0.5 * (x_mode - mu) * (x_mode - mu) / sigma2)
      / dgamma(1.0 / x_mode, /* shape = */ phi + 2.0,
	       /* scale = */ 1.0 / theta, /* give_log = */ 0);
	       
    do {

      *h_new = 1.0 / rgamma(/* shape = */ phi + 2.0,
			    /* scale = */ 1.0 / theta);

      p2 = R_pow(*h_new, -1.5) * exp(-0.5 * (*r) * (*r) / (*h_new))
	* exp(-0.5 * (log(*h_new) - mu) * (log(*h_new) - mu) / sigma2);

      q2 = dgamma(1.0 / (*h_new), /* shape = */ phi + 2.0,
		  /* scale = */ 1.0 / theta, /* give_log = */ 0);

    } while (unif_rand() > p2 / (c * q2));
    
      p1 = R_pow(*h, -1.5) * exp(-0.5 * (*r) * (*r) / (*h))
	* exp(-0.5 * (log(*h) - mu) * (log(*h) - mu) / sigma2);

    q1 = dgamma(1.0 / (*h), /* shape = */ phi + 2.0,
		/* scale = */ 1.0 / theta, /* give_log = */ 0);

    num   = p2 / uri_min(p2, c * q2);
    denom = p1 / uri_min(p1, c * q1);

    *h_new = unif_rand() < num / denom ? (*h_new) : (*h);
  }

  /* By now h_new has been incremented (*T) - 1 times */

  h_new[0]
    = exp((*alpha) + (*delta) * log(h_new[-1]) + (*sigv) * norm_rand());

  h_new -= (*T - 1);

  h_new[0]
    = exp((*alpha) + (*delta) * log(h_new[1]) + (*sigv) * norm_rand());

  PutRNGstate();

}







SEXP get_init_params_test(SEXP r)
{
  int T, num_protects = 0;
  SEXP alpha, delta, sigv, lambda, N, ans, h;
  carray M, m;

  T = length(r);

  M = make_zero_matrix(2, 2);
  m = make_zero_matrix(2, 1);

  PROTECT(ans = NEW_LIST(8));
  PROTECT(alpha = NEW_NUMERIC(1));
  PROTECT(delta = NEW_NUMERIC(1));
  PROTECT(sigv = NEW_NUMERIC(1));
  PROTECT(lambda = NEW_NUMERIC(1));
  PROTECT(N = NEW_NUMERIC(1));
  PROTECT(h = NEW_NUMERIC(T));

  num_protects += 7;

  get_init_params(T, REAL(alpha), REAL(delta),
		  REAL(sigv), m, M, REAL(N),
		  REAL(lambda), REAL(h));
  
  SET_ELT(ans, 0, alpha);
  SET_ELT(ans, 1, delta);
  SET_ELT(ans, 2, sigv);
  SET_ELT(ans, 3, carray_to_robj(m));
  SET_ELT(ans, 4, carray_to_robj(M));
  SET_ELT(ans, 5, N);
  SET_ELT(ans, 6, lambda);
  SET_ELT(ans, 7, h);

  UNPROTECT(num_protects);

  return ans;
}


SEXP update_params_test(SEXP h, SEXP m, SEXP M, SEXP lambda, SEXP N)
{
  
  int T, num_protects = 0;
  SEXP alpha, delta, sigv, ans, lambda_out, N_out;
  carray M_out, m_out;

  T = length(h);

  M_out = robj_to_carray(M, /* dup = */ true);
  m_out = robj_to_carray(m, /* dup = */ true);

  PROTECT(ans = NEW_LIST(7));
  PROTECT(alpha = NEW_NUMERIC(1));
  PROTECT(delta = NEW_NUMERIC(1));
  PROTECT(sigv = NEW_NUMERIC(1));
  PROTECT(lambda_out = NEW_NUMERIC(1));
  PROTECT(N_out = NEW_NUMERIC(1));

  num_protects += 6;

  *REAL(N_out) = *REAL(N);
  *REAL(lambda_out) = *REAL(lambda);

  update_params(T, REAL(alpha), REAL(delta),
		REAL(sigv), m_out, M_out, REAL(N_out),
		REAL(lambda_out), REAL(h));

  SET_ELT(ans, 0, alpha);
  SET_ELT(ans, 1, delta);
  SET_ELT(ans, 2, sigv);
  SET_ELT(ans, 3, carray_to_robj(m_out));
  SET_ELT(ans, 4, carray_to_robj(M_out));
  SET_ELT(ans, 5, N_out);
  SET_ELT(ans, 6, lambda_out);
  SET_ELT(ans, 7, h);

  UNPROTECT(num_protects);

  return ans;

}




#define STATE_LEN 9
#define ANSLEN 9
SEXP algo3(SEXP r, SEXP info, SEXP initState)
{

  int numprot = 0;
  int maxit, drawSuccess, diffusePrior;
  long i, it, t, numLoops, initLoops, T, n;
  double *hPtr = NULL, *rPtr = NULL;
  double alpha, delta, sigv, N, lambda, hcand, c, phi, mu, theta, qMode, sigma;
  double rho, sig1, sig2, e1, e2;
  double p, q, pnew, qnew, fnew, f;
  char *ansNames[ANSLEN] = {"alpha", "delta", "sigv", 
                            "m", "M", "N", 
                            "lambda", "h", "hTerm"};
  carray m, M;
  carray y, X, Xt, Mi, Mnew, res, b, Xb, XtX, XtXi, Xty;
  carray tmp_2x2, tmp_1x2, tmp_2x1, tmp_1x1;
  SEXP ans, alphaOut, deltaOut, sigvOut, hOut;
  SEXP newState, h;

  diffusePrior = asInteger(getListElt(info, "diffusePrior"));

  PROT2(initState = getInitState(r, initState, diffusePrior), numprot);

  if (!isNewList(initState) || STATE_LEN != length(initState))
    error ("initState must be a list");

  hUpdate1_chk(&initState, &numprot);
  
  numLoops = asInteger(getListElt(info, "numLoops"));
  initLoops = asInteger(getListElt(info, "initLoops"));
  maxit = asInteger(getListElt(info, "maxit"));
  alpha = asReal(getListElt(initState, "alpha"));
  delta = asReal(getListElt(initState, "delta"));
  sigv  = asReal(getListElt(initState, "sigv"));
  m = sexp_to_carray(getListElt(initState, "m"), 1);
  M = sexp_to_carray(getListElt(initState, "M"), 1);
  N = asReal(getListElt(initState, "N"));
  lambda = asReal(getListElt(initState,"lambda"));
  PROT2(h = getListElt(initState, "h"), numprot);
  PROT2(r = getListElt(initState, "r"), numprot);
  T = length(h);
  hPtr = REAL(h);
  rPtr = REAL(r);

  PROT2(alphaOut = NEW_NUMERIC(numLoops), numprot);
  PROT2(deltaOut = NEW_NUMERIC(numLoops), numprot);
  PROT2(sigvOut = NEW_NUMERIC(numLoops), numprot);
  PROT2(hOut = NEW_NUMERIC(numLoops), numprot);
  PROT2(ans = NEW_LIST(ANSLEN), numprot);
  PROT2(newState = NEW_LIST(STATE_LEN), numprot);

  y       = make_zero_matrix(T - 1, 1);
  res     = make_zero_matrix(T - 1, 1);
  X       = make_zero_matrix(T - 1, 2);
  Xb      = make_zero_matrix(T - 1, 1);
  Xt      = make_zero_matrix(2, T - 1);
  Mi      = make_zero_matrix(2, 2);
  Mnew    = make_zero_matrix(2, 2);
  b       = make_zero_matrix(2, 1);
  XtX     = make_zero_matrix(2, 2);
  XtXi    = make_zero_matrix(2, 2);
  Xty     = make_zero_matrix(2, 1);
  tmp_2x2 = make_zero_matrix(2, 2);
  tmp_1x2 = make_zero_matrix(1, 2);
  tmp_2x1 = make_zero_matrix(2, 1);
  tmp_1x1 = make_zero_matrix(1, 1);

  for (t = 0; t < T-1; t++)
    ARRAY2(Xt)[0][t] = 1;

  GetRNGstate();

  for (n = -initLoops; n < numLoops; n++) {

    /* 
    **  Update h 
    */

    sigma = sigv/sqrt(1 + delta*delta);
    phi = exp(sigv*sigv)/expm1(sigv*sigv) + 1 + 0.5;


    for (i = 1; i < T-1; i++) {
    
      mu = alpha * (1-delta) + delta * (log(hPtr[i-1]) + log(hPtr[i+1]));
      mu /= (1 + delta*delta);
      theta = 0.5*SQR(rPtr[i]) + (phi-1)*exp(mu + 0.5*SQR(sigma));
      qMode = theta / (phi+1);
    
      c = 1.1 * R_pow(qMode, -1.5) * exp(-0.5*SQR(rPtr[i])/qMode);
      c *= exp(-0.5 * SQR(log(qMode) - mu) / SQR(sigma));
      c /= dgamma(1/qMode, /* shape = */ phi + 2, /* scale = */ 1/theta, 0);

      it = drawSuccess = 0;

      do {

        hcand = 1/rgamma(/*shape = */ phi, /* scale = */ 1/theta);
 
        pnew  = R_pow(hcand, -1.5) * exp(-0.5*SQR(rPtr[i])/hcand);
        pnew *= exp(-0.5 * SQR(log(hcand) - mu) / SQR(sigma));
        qnew  = dgamma(1/hcand, /* shape = */ phi + 2, /* scale = */ 1/theta, 0);

        drawSuccess = unif_rand() < pnew / (c*qnew);
        it += 1;

      } while (!drawSuccess && it < maxit);

      if (drawSuccess) {

        p = R_pow(hPtr[i], -1.5) * exp(-0.5*SQR(rPtr[i])/hPtr[i]);
        p *= exp(-0.5 * SQR(log(hPtr[i]) - mu) / SQR(sigma));
        q = dgamma(1/hPtr[i], phi + 2, 1/theta, 0);
    
        fnew = MIN(pnew, c * qnew);
        f = MIN(p, c * q);

        hPtr[i] = unif_rand() < (pnew/fnew)/(p/f) ? hcand : hPtr[i];

      }

    }

    hPtr[0]   = exp(alpha + delta * log(hPtr[1])   + sigv * norm_rand());
    hPtr[T-1] = exp(alpha + delta * log(hPtr[T-2]) + sigv * norm_rand());



    /* 
    **  Update alpha, delta, sigv 
    */

    for (t = 0; t < T-1; t++) {
      ARRAY2(Xt)[1][t] = log(hPtr[t]);
      ARRAY1(y)[t] = log(hPtr[t+1]);
    }

    transpose_matrix(Xt, X);
    matrix_prod_ns(X, X, 1, 0, XtX);
    matrix_inverse_2x2(XtX, XtXi);
    matrix_prod_ns(X, y, 1, 0, Xty);
    matrix_prod_ns(XtXi, Xty, 0, 0, b);
    matrix_prod_ns(X, b, 0, 0, Xb);

    /* Update N and lambda */
    N += T - 1;

    array_op(y, Xb, '-', res);  

    if ( ! matrix_inverse_2x2(M, Mi) )
      set_array_to_zero(Mi);
  
    array_op(XtX, Mi, '+', Mnew); /* Mnew <- XtX + Mi */
    matrix_inverse_2x2(Mnew, Mnew); /* Mnew <- (XtX + Mi)^-1 */
    matrix_prod_ns(Mi, Mnew, 0, 0, tmp_2x2); /* tmp_2x2 <- Mi * Mnew*/
    matrix_prod(tmp_2x2, XtX, 0, 0, tmp_2x2); /* tmp_2x2 <- Mi * Mnew * t(X) * X */
    array_op(b, m, '-', tmp_2x1);
    /* tmp_1x2 <- t(b-m) * Mi * Mnew * t(X) * X : */
    matrix_prod_ns(tmp_2x1, tmp_2x2, 1, 0, tmp_1x2);
    /* tmp_1x1 <- t(b-m) * Mi * Mnew * t(X) * X * (b-m) : */
    matrix_prod_ns(tmp_1x2, tmp_2x1, 0, 0, tmp_1x1);

    lambda += sumsq(res) + ARRAY1(tmp_1x1)[0];

  
    /* Draw new sigv */  
    sigv = sqrt(1.0 / rgamma(N/2, 2/lambda));
  
    /* Draw new alpha, delta conditional on new sigma.v */
    matrix_prod(Mi, m, 0, 0, m);
    array_op(m, Xty, '+', m);
    matrix_prod(Mnew, m, 0, 0, m);
    memcpy(ARRAY1(M), ARRAY1(Mnew), 4 * sizeof(double));
    matrix_inverse_2x2(M, Mi);

    /* Draw new alpha, delta conditional on new sigv */

    sig1 = sigv * sqrt(ARRAY1(M)[0]);
    sig2 = sigv * sqrt(ARRAY1(M)[3]);
    rho  = sigv * sigv * ARRAY1(M)[1] / (sig1 * sig2);

    do {

      e1 = norm_rand();
      e2 = rho * e1 + sqrt(1 - rho * rho) * norm_rand();

      alpha = ARRAY1(m)[0] + sig1 * e1;
      delta = ARRAY1(m)[1] + sig2 * e2;

      /* Continue until |delta| < 1 for stationarity */

    } while (ABS(delta) >= 1);


    if (n >= 0) {
      
      REAL(alphaOut)[n] = alpha;
      REAL(deltaOut)[n] = delta;
      REAL(sigvOut)[n] = sigv;
      REAL(hOut)[n] = hPtr[T-2];

    }

  }
  
  PutRNGstate();

  SET_ELT(ans, 0, alphaOut);
  SET_ELT(ans, 1, deltaOut);
  SET_ELT(ans, 2, sigvOut);
  SET_ELT(ans, 3, carray_to_sexp(m));
  SET_ELT(ans, 4, carray_to_sexp(M));
  SET_ELT(ans, 5, ScalarReal(N));
  SET_ELT(ans, 6, ScalarReal(lambda));
  SET_ELT(ans, 7, h);
  SET_ELT(ans, 8, hOut);

  set_names(ans, ansNames);

  UNPROT2;

  return ans;

}
#undef STATE_LEN
#undef ANSLEN




#endif
