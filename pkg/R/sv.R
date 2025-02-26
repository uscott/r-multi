

gibbs.uri = function(pxy, pyx, n = 1e2, yInit = rep(0, n), stopPt = 1e3) {
  ## pxy = function(n, y) is a function that gives you n independent
  ## draws of x conditional on y.  Similarly pyx = function(n, x).

  y = yInit
  x = pxy(n, y)
  
  for (k in seq(stopPt)) {
    y = pyx(n, x)
    x = pxy(n, y)
  }
  
  x
}




######################################################################








getInitState = function(r) {
  ## r = returns time series
 
  TT = length(r)

  ## Step 1
  ## Draw alpha, delta, sigv from priors

  ## Prior for sigv is inverse gamma
  ## Notation for Bayesian linear regression
  ## to be consistent with Hamilton pp. 355, 356.
  N = 4 ## Corresponds to N/2 = 2 = nu_0 + 1
  lambda = 0.005^2 ## Corresponds to lambda = nu_0 * s_0^2
  sigv = sqrt(rigamma(1, N, lambda))
  
  delta = runif(1, min = -1, max = 1)
  alpha = rnorm(1, sd = 1)

  ## Step 2
  ## Simulate log(h) as AR(1) process.
  logh = arima.sim(model = list(ar = delta), n = TT, sd = sigv)
  logh = alpha/(1 - delta) + logh
  
  y = logh[-1]
  X = cbind(1, logh[-TT])

  m1 = rep(0, 2)
  M1i = matrix(0, 2, 2)
  b = as.double(solve(crossprod(X)) %*% crossprod(X,y))
  
  ## Bayesian update N and lambda
  N = N + TT - 1
  res = y - X %*% b
  lambda = lambda + crossprod(res)

  ## Bayesian update sigv
  sigv = sqrt(rigamma(1, N, lambda))
  
  ## Draw new alpha, delta conditional on new sigv
  M2i = crossprod(X) ## M2i is a 2x2 matrix
  M2  = solve(M2i)
  m2  = M2 %*% (M1i %*% m1 + crossprod(X,y))

  repeat {
    
    beta = rmvnorm(100, mean = m2, sigma = sigv^2 * M2)
    ix = which(abs(beta[,2]) < 1)
    ## Continue until |beta[2]| < 1 for stationarity
    if (length(ix)) {
      beta = beta[ix[1], ]
      break
    }
  }

  alpha = beta[1]
  delta = beta[2]

  list(alpha = alpha, delta = delta, sigv = sigv,
       m = m2, M = M2, N = N, lambda = lambda, h = exp(logh), r = r)
}





updateParams = function(currentState) {

  .Call("updateParams1", currentState)
  
}





updateParams.old = function(currentState) {

  h = currentState$h
  m = currentState$m
  M = currentState$M
  lambda = currentState$lambda
  N = currentState$N
  
  TT = length(h)

  y = log(h[-1])
  X = cbind(1, log(h[-TT]))
  XtX = crossprod(X)
  Xty = crossprod(X, y)
  
  b = as.double(solve(XtX) %*% Xty)
  
  ## Update N and lambda
  N = N + TT - 1
  
  res = y - X %*% b
  Mi = solve(M)
  Mnew = solve(Mi + XtX)
  
  lambda = lambda + crossprod(res)
  lambda = lambda + t(b - m) %*% Mi %*% Mnew %*% XtX %*% (b - m)

  ## Draw new sigv
  sigv = sqrt(rigamma(1, N, lambda))
  
  ## Draw new alpha, delta conditional on new sigv
  m = as.double(Mnew %*% (Mi %*% m + Xty))
  M = Mnew

  repeat {
    
    beta = rmvnorm(100, mean = m, sigma = sigv^2 * M)
    ix = which(abs(beta[,2]) < 1)
    
    ## Continue until |delta| < 1 for stationarity
    if (length(ix)) {
      beta = beta[ix[1], ]
      break
    }
    
  }
  
  alpha = beta[1]
  delta = beta[2]

  currentState$alpha = alpha
  currentState$delta = delta
  currentState$sigv = sigv
  currentState$m = m
  currentState$M = M
  currentState$N = N
  currentState$lambda = lambda
  
  currentState
}







hUpdate.old = function(currentState) {

  alpha = currentState$alpha
  delta = currentState$delta
  sigv = currentState$sigv
  h = currentState$h
  r = currentState$r
  
  TT = length(h)
  stopifnot(TT == length(r))
  
  mu = (alpha * (1 - delta) + delta * (log(lag2(h, -1)) + log(lag2(h, 1))))
  mu = mu / (1+delta^2)

  sigma = sqrt(sigv^2 / (1 + delta^2))

  phiLN = (1 - 2 * exp(sigma^2))/(1 - exp(sigma^2))
  phi = -0.5 + phiLN - 1
  thetaLN = (phiLN - 0.5) * exp(mu + 0.5 * sigma^2)
  theta1 = 0.5 * r^2
  theta = theta1 + thetaLN
  ## shape = -phi, scale = 1 / rate = 1 / theta
  ## q = pdf of gamma(shape, scale)
  ## mode of q is at theta/(phi + 1)
  xMode = theta / (phi + 1)

  cc = 1.1 * xMode^(-1.5)
  cc = cc * exp(-0.5*r^2/xMode) * exp(-0.5 * (log(xMode) - mu)^2 / sigma^2)
  cc = cc / dgamma(1/xMode, shape = phi, rate = theta)

  ix1 = seq(2, TT - 1)

  while (length(ix1)) {
    ## Draw from pdf q:
    h1 = h[ix1]
    h2 = 1/rgamma(length(ix1), shape = phi, rate = theta[ix1])
    u = runif(length(ix1))

    p2 = h2^(-1.5) * exp(-0.5 * r[ix1]^2 / h2)
    p2 = p2 * exp(-0.5 * (log(h2) - mu[ix1])^2 / sigma^2)
    
    q2 = dgamma(1/h2, shape = phi + 2, rate = theta[ix1])
    
    ix3 = which(u < p2 / (q2 * cc[ix1]))

    if (length(ix3)) {

      v = runif(length(ix3))
      ix2 = ix1[ix3]
      
      p1 = h1[ix3]^(-1.5) * exp(-0.5 * r[ix2]^2 / h1[ix3])
      p1 = p1 * exp(-0.5 * (log(h1[ix3]) - mu[ix2])^2 / sigma^2)
      
      q1 = dgamma(1/h1[ix3], shape = phi + 2, rate = theta[ix2])

      num = (p2 / pmin(p2, cc[ix1] * q2))[ix3]
      denom = p1 / pmin(p1, cc[ix2] * q1)
      
      ix4 = which(v < num / denom)

      if (length(ix4))
        h1[ix4] = h2[ix4]
      
      h[ix1[ix3]] = h1[ix3]

      ix1 = setdiff(ix1, ix2)
    }
    
  }

  h[c(1, TT)] = exp(alpha+delta*log(h[c(2, TT-1)]) + sigv*rnorm(2))

  currentState$h = h
  currentState
  
}





hUpdate = function(currentState, maxit) {

  .Call("hUpdate1", currentState, maxit);
  
}







algo1 = function(r, n = 2e2, initLoops = 10, output.h = FALSE, maxit=10) {

  require(mvtnorm)
  ## r = returns time series

  TT = length(r)
  alphaOut = deltaOut = sigvOut = rep(NA, n)
  mOut = matrix(NA, 2, n)
  MdistChange = rep(NA, n)
  if (output.h)
    hOut = matrix(NA, TT, n)
  else
    hOut = NULL
  
  r = c(NA, r, NA)

  ## Draw initial parameter values
  currentState = getInitState(r)
  
  ## Now that we've drawn alpha, delta and sigv apply
  ## MCMC method

  currentState = hUpdate(currentState,maxit=maxit)
  
  ################################################
  ## Repeat the above stuff a number of         ##
  ## times to get to stationary state hopefully ##
  ###############################################
  for (j in seq(initLoops)) {
    
    ## Now that we've drawn alpha, delta and sigv apply
    ## MCMC method

    currentState = updateParams(currentState)
    currentState = hUpdate(currentState,maxit=maxit)
    
  }

  ######################################################################
  ################## Continue sampling from stationary distribtuion ####
  ######################################################################

  Mold = currentState$M
  
  for (j in 1:n) {
    
    ## Now that we've drawn alpha, delta and sigv apply
    ## MCMC method

    currentState = updateParams(currentState)

    alphaOut[j] = currentState$alpha
    deltaOut[j] = currentState$delta
    sigvOut[j]  = currentState$sigv
    mOut[, j]   = currentState$m
    MdistChange[j] = sum((currentState$M - Mold)^2)

    Mold = currentState$M
    
    currentState  = hUpdate(currentState,maxit=maxit)

    if (output.h)
      hOut[, j] = currentState$h[1 + seq(TT)]
  }

  if (output.h)
    hOut = ts(hOut)
  
  list(alpha = alphaOut, delta = deltaOut,
       sigv = sigvOut, m = mOut, M = MdistChange, h = hOut)
}




algo2 = function(r, n = 2e2, initLoops = 10, output.h = FALSE, maxit=10) {

  info = list(n = n, initLoops = initLoops,
    output.h = output.h, maxit = maxit)
  
  r = c(NA, r, NA)

  ## Draw initial parameter values
  currentState = getInitState(r)
  
  ## Now that we've drawn alpha, delta and sigv apply
  ## MCMC method

  currentState = hUpdate(currentState,maxit=maxit)
  
  .Call("algo2", currentState, info, PACKAGE = "uri")
  
}



algo3 = function(r, loops = 2e2, preLoops = 0,
  diffusePrior = FALSE, initState = NULL, maxit=100) {

  info = list(numLoops=loops, initLoops=preLoops,
    maxit=maxit, diffusePrior = diffusePrior)
  
  .Call("algo3", r, info, initState, PACKAGE = "uri")
  
}




svSim1 = function(len, alpha, delta, sigv) {

  logh = arima.sim(list(ar=delta), n = len, n.start = len, sd = sigv)
  logh = logh + alpha/(1-delta)

  h = exp(logh)
  r = sqrt(h) * rnorm(len)

  data.frame(r = r, h = h)
  
}



rigamma = function(n, N, lambda)
  1/rgamma(n, shape = N/2, rate = lambda/2)



bayesLmTest = function(beta, sig, N = 100, trials = 100,
  option1 = FALSE) {

  require(mvtnorm)
  
  k = length(beta) - 1

  stopifnot(length(sig) == 1)

  betaOut = matrix(NA, trials, k+1)
  sigOut = rep(NA, trials)

  sigOut[1] = exp(runif(1))
  betaOut[1,] = runif(k+1)

  n = lambda = 0
  Mi = matrix(0, k+1, k+1)
  m = numeric(k+1)
  
  for (i in 1:trials) {

    X = cbind(1, matrix(rnorm(k*N), N, k))
    
    y = X %*% beta + sig * rnorm(N)

    g = lm(y ~ 0 + X)
    olsRes = g$res
    b = g$coef
    
    n = n + N

    Minew = crossprod(X)+Mi
    Mnew = solve(Minew)
    
    ss = crossprod(olsRes)
    lambda = lambda + ss
    
    if (option1)
      lambda = lambda + t(b-m)%*%Mi%*%Mnew%*%crossprod(X)%*%(b-m)
    
    m = Mnew %*% (Mi%*%m + crossprod(X,y))
    Mi = Minew

    
    sigOut[i] = sqrt(rigamma(1, n, lambda))
    betaOut[i,] = rmvnorm(1, m, sigOut[i]^2 * solve(Mi))
    
  }

  list(beta = betaOut, sig = sigOut)
  
}





