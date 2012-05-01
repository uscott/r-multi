





kalman1 = function(y, Htr, f, Q, R, x = NULL, Atr = NULL)
  ## Estimates loglikelihood and produces forecasts for
  ## state-space model with constant coefficients.
  ## Based on treatment of KF in Hamilton.
{

  .Call("kalman1", y, Htr, f, Q, R, x, Atr, PACKAGE = "kalman.uri")
  
}





ar.fit = function(y, p = 1, init = NULL)
{

  y = as.matrix(y)
  y = y - mean(y)
  n = length(y)
  
  Htr    = t(c(1, rep(0, p - 1)))
  Flower = cbind(diag(p - 1), 0) ## Full matrix F is p x p.
  Qmat   = matrix(0, p, p)
  R      = matrix(0, 1, 1)
  x      = matrix(1, n, 1)
  
  f = function(par) {

    Qmat[1] = par[p + 1]
    
    kalman1(y, Htr, rbind(par[1:p], Flower), Qmat, R, x = NULL, Atr = NULL)$ll
    
  }


  optim.uri(f,
            npars    = p + 1,
            init     = runif(3),
            method   = c("ga", "gr"),
            minimize = FALSE)
  
}






kal1.eval = function(y, a0, a1, b1, gamma, phi, c0, sigma, rho)
{

  par = c(a0, a1, b1, gamma, phi, c0, sigma, rho)

  .Call("kal1", y, par, PACKAGE = "kalman.uri")

}






kalman2 = function(y, Htr, Fmat, Q, R, x = NULL, Atr = NULL) {
  ## Estimates loglikelihood and produces forecasts for
  ## state-space model with constant coefficients.
  ## Based on treatment of KF in Hamilton.

  stopifnot(all(is.matrix(Htr), is.matrix(Fmat), is.matrix(Q),
                is.matrix(R), is.matrix(Atr) | is.null(Atr)))

  y = as.matrix(y)
  
  n = ncol(y)
  m = nrow(y)
  r = ncol(Htr)

  stopifnot(all(r == nrow(Fmat), r == ncol(Fmat), n == nrow(Htr),
                r == nrow(Q), r == ncol(Q), n == nrow(R),
                n == ncol(R)))

  if (!is.null(x)) {

    x = as.matrix(x)
    k = ncol(x)

    stopifnot(is.matrix(Atr))
    stopifnot(n == nrow(Atr) & k == ncol(Atr))
    
  } else {

    k = 1
    x = matrix(0, m, k)
    Atr = matrix(0, n, k)
    
  }

  xi   = matrix(0, m, r)         ## xi[s,] = xi(s | s - 1)
  res  = matrix(NA, m, n)
  gain = array(0, c(m, r, n))
  P    = array(0, c(m, r, r))
  H    = t(Htr)
  Ftr  = t(Fmat)
  
  
  P[1,,] = matrix(solve(diag(r^2) - (Fmat %x% Fmat)) %*% as.vector(Q), r, r)

  ll = n * m * log(2 * pi)
  
  for (s in seq(m - 1)) {

    B  = Htr %*% P[s,,] %*% H + R
    ll = ll + log(det(B))
    B  = solve(B)
    
    K           = Fmat %*% P[s,,] %*% H %*% B
    dim(K)      = c(r, n)
    res[s, ]    = y[s,] - Atr %*% x[s,] - Htr %*% xi[s, ]
    ll          = ll + t(res[s,]) %*% B %*% res[s,]
    xi[s + 1, ] = Fmat %*% xi[s,] + K %*% res[s, ]

    B = Fmat - K %*% Htr
    
    P[s + 1,,] = B %*% P[s, ,] %*% t(B) + K %*% R %*% t(K) + Q

    gain[s, ,] = K
  }

  B  = Htr %*% P[m,,] %*% H + R
  ll = ll + log(det(B))
  B  = solve(B)

  K           = Fmat %*% P[m,,] %*% H %*% B
  dim(K)      = c(r, n)
  res[m, ]    = y[m,] - Atr %*% x[m,] - Htr %*% xi[m, ]
  ll          = ll + t(res[m,]) %*% B %*% res[m,]  
  
  ll = -0.5 * ll

  gain[m, , ] = K
  
  list(xi = xi, res = res, P = P, K = gain, ll = ll)
}









