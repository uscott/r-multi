

gls.sumsq = function(y, x, beta, phi)
  ##
  ##  Description: Computes the residuals and sum of squared
  ##    residuals for the GLM y = x %*% beta + u where the u's
  ##    follow an AR process with coefficients phi.  The
  ##    residuals returned are the residuals from the AR
  ##    process.
  ##
  ##  Returns:
  ##    List of length 2 with elements
  ##        sum.of.squares = sum of squared residuals
  ##        residuals      = residuals from AR process.
  ##

{
  
  stopifnot(length(y) == nrow(x), length(beta) == ncol(x),
            length(y) > length(phi))
  
  .Call("gls_sumsq2", as.double(y), as.matrix(x),
        as.double(beta), as.double(phi))
}




gls.uri = function(y, x,
  intercept.add  = TRUE,
  p              = NULL,
  beta.init      = NULL,
  phi.init       = NULL,
  na.rm          = TRUE,
  weights        = NULL,
  method         = c("ga", "N", "ga"),
  bootstrap      = NULL, ...)
  ##
  ##  Description: Fits the GLM y = x %*% beta + e where the e's follow
  ##    an AR(p) process with coefficients phi.
  ##
  ##  Returns:
  ##    List of length 4 with elements
  ##        beta = fitted values of beta
  ##        phi  = fitted values of phi (coefficients of AR process)
  ##        res  = fitted residuals from AR process
  ##

{

  y = as.double(y)
  x = as.matrix(x)

  if (is.null(weights))
    weights = rep(1, length(y))

  stopifnot (length(weights) == length(y))
  
  if (na.rm) {

    ix = !is.na(y) & !apply(x, 1, function(z) any(is.na(z)))
    y  = y[ix]
    weights = weights[ix]
    x  = x[ix, , drop = FALSE]
    
  }
  
  if (is.null(colnames(x)))
    x.names = paste("x", seq(ncol(x)), sep = "")
  else
    x.names = colnames(x)
  
  if (intercept.add) {
    
    x       = cbind(1, x)
    x.names = c("intercept", x.names)
    
  }

  colnames(x) = x.names
  
  m = nrow(x)
  n = ncol(x)
  
  d        = data.frame(cbind(y, x))
  names(d) = c("y", x.names)

  
  if (!is.null(weights)) {
    
    d = d * sqrt(weights)
    y = y * sqrt(weights)
    x = x * sqrt(weights)

  }
  
  if (is.null(phi.init)) {
    
    g        = lm(y ~ 0 + ., data = d, ...)
    phi.init = ar(g$res, method = "mle", ord = p)$ar
    
  }

  if (is.null(beta.init))
    beta.init = lm(y ~ 0 + ., data = d, ...)$coef

  p = length(phi.init)

  stopifnot (n == length(beta.init))
  
  f = function(coefs) {
    
    .Call("gls_sumsq2",
          y, x,
          coefs[1:n],
          coefs[n + 1:p])$sum.of.squares
    
  }


  opt = optim.uri(f,
    init   = c(beta.init, phi.init),
    method = method,
    ...)

  j = which.max(opt$value)
  
  coefs = opt$par[, j]
  
  beta = coefs[1:n]
  phi  = coefs[n + 1:p]

  res = gls.sumsq(y, x, beta, phi)$res

  if (is.null(bootstrap))
    boot.coef = NULL
  else {
    
    boot.coef = matrix(NA, bootstrap, n)
    
    Y = y
    X = x

    for (i in 1:p) {
      
      Y[-seq(i)]   = Y[-seq(i)]   - phi[i] * y[seq(m - i)]
      X[-seq(i), ] = X[-seq(i), ] - phi[i] * x[seq(m - i), ]
      
    }

    h = lm(Y ~ 0 + X, ...)
    
    for (b in seq(bootstrap))
      boot.coef[b, ] = lm(h$fit + sample(h$res, rep = T) ~ 0 + X)$coef
    
  }
  
  list(beta = beta, phi = phi, res = res, y = y, x = x, boot.coef = boot.coef)
}










gls.uri.2 = function(y, x, intercept.add = TRUE, weights = NULL, method = "N", ...)
{

  y = as.double(y)
  x = as.matrix(x)
  m = nrow(x)
  
  if (intercept) x = cbind(rep(1, m), x)
  
  if (!is.null(weights)) {
    x = x * sqrt(weights)
    y = y * sqrt(weights)
  }

  g = lm(y ~ 0 + x, x = TRUE, y = TRUE)
  phi = ar(g$res, method = "mle")$ar
  p = length(phi)
  
  if (p) {

    y = y2 = g$y
    x = x2 = g$x
    
    for (i in seq(p)) {
      y2 = y2 - phi[i] * lag2(y, i)
      x2 = x2 - phi[i] * apply(x, 2, lag2, n = i)
    }

    g = lm(y2 ~ 0 + x2, x = TRUE, y = TRUE)
  }

  if (!is.null(weights)) g$res = g$res / sqrt(weights)
  
  list(beta = g$coef, phi = phi, res = g$res, regression = g)
}



