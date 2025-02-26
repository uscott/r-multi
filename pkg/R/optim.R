


optim.uri = function(fn,
  method         = "ga",
  init           = NULL,
  npar           = NULL,
  maxit          = rep(20, length(method)),
  lower          = NULL,
  upper          = NULL,
  minimize       = TRUE,
  tol            = 1e-8,
  enforce.bounds = FALSE,
  ...)
  ##
  ##  Description: General purpose optimization function.
  ##    Basically acts as a wrapper for various other
  ##    optimization functions.
  ##

{

  
  if (is.matrix(init)) {
    
    stopifnot(is.double(init))
    
  } else if(is.double(init) & is.null(dim(init))) {
    
    stopifnot(is.double(init))
    init = as.matrix(init)
    
  }
  
  method = substr(method, 1, 2)

  M           = length(method)
  maxit       = rep(maxit, M)[1:M]
  par         = init
  md          = Inf
  itnum       = 0
  valold      = Inf
  convergence = 1
  
  repeat {
    
    ## Loop through vector of optimization methods
    for (k in 1:M) {

      ## Choose optimization method
      if      ("ga" == method[k])
        optim.fn = optim.ga1
      else if ("gr" == method[k])
        optim.fn = optim.grad2
      else
        optim.fn = optim2

      ## Do optimization
      ans = try(
        optim.fn(fn             = fn,
                 init           = par,
                 npar           = npar,
                 lower          = lower,
                 upper          = upper,
                 minimize       = minimize,
                 maxit          = maxit[k],
                 enforce.bounds = enforce.bounds,
                 ...)
        )

      ## Update par & val, figure out max dist betweeen pairs of pars
      if (!is.null(ans$par) & !is.null(ans$val)) {

        ind = which.not.na(ans$val)
        par = as.matrix(ans$par[, ind, drop = FALSE])
        val = ans$val[ind]
        
        if (length(val)) { # Figure out max distance between val and valold

          md = max(abs(outer(valold, val, "-")), na.rm = TRUE)

          valold = val
          
        }

        convergence = as.integer(identical(md >= tol, TRUE))
        
      }

      if (!convergence) # Finish optimization early
        break

    }

    itnum = itnum + 1

    if (!convergence || itnum >= maxit[1])
      break
    
  }
  
  ans$convergence = convergence

  
  if (!convergence) {
    ## Only need 1 parameter set & corresponding value

    j = ifelse(minimize, which.min(ans$value), which.max(ans$value))
    
    ans$par   = ans$par[, j, drop = FALSE]
    ans$value = ans$value[j]

  }
    
  ans
}






optim2 = function(fn,
  init,
  method         = "Nelder-Mead",
  lower          = NULL,
  upper          = NULL,
  minimize       = TRUE,
  npar           = NULL,
  maxit          = NULL,
  tol            = NULL,
  enforce.bounds = FALSE,
  ...) 
  ##
  ##  Description: Generalization of optim function from base
  ##    package that simply allows for more than one initial
  ##    value.
  ##
{

  init  = as.matrix(init)

  if (is.null(maxit))
    maxit = 500

  if (is.null(tol))
    tol = sqrt(.Machine$double.eps)
  
  if (!is.null(npar))
    stopifnot (npar == nrow(init))
  
  npar  = nrow(init)
  Ninit = ncol(init)
  
  if((is.null(lower) & is.null(upper)) | !enforce.bounds) {
    
    tf = tfi = function(x) x
    f  = fn
    
  } else {
    
    tf  = transform.par(lower, upper, TRUE)
    tfi = transform.par(lower, upper)
    f   = function(x) fn(tf(x))
    
    ## tf  transforms unconstrained to constrained parameters
    ## tfi = inverse of tf, transforms constrained to unconstrained
  }

  
  control = list(fnscale = ifelse(minimize, 1, -1), maxit = maxit, reltol = tol)
    
  
  ans = list(par        = matrix(NA, npar, Ninit),
             value      = rep(NA, Ninit),
             convergence = rep(NA, Ninit),
             info        = list())

  init = tfi(init) # transform initial params to unconstrained

  
  for (k in seq(Ninit)) {
    
    tmp = try(optim(init[, k], f, method = method, control = control, ...))

    if (!is.null(tmp$val) & !is.null(tmp$par)) {
      
      ans$value[k]       = tmp$val
      ans$par[, k]       = tf(as.matrix(tmp$par))
      ans$convergence[k]  = tmp$convergence
      
    }

    ans$info[k] = list(tmp)
    
  }

  ans
}






optim.ga1 = function(fn,
  init = NULL,
  npar = NULL,
  pop.size = 5,
  maxit = 5,
  lower = NULL,
  upper = NULL,
  minimize = TRUE,
  tol = sqrt(.Machine$double.eps),
  enforce.bounds = FALSE,
  ...) 
  ##
  ##  Description: Optimizes function fn using genetic algorithm.
  ##
{
  init = init.population(pop.size = pop.size,
    npar  = npar,
    lower = lower,
    upper = upper,
    init  = init)


    if((is.null(upper) & is.null(lower)) | !enforce.bounds) {
      
      tf  = tfi = function(x) x
      f = fn
      
    } else {
      
      tf  = transform.par(lower, upper, TRUE)
      tfi = transform.par(lower, upper)
      f = function(x) fn(tf(x))
      
      ## tf  transforms unconstrained to constrained parameters
      ## tfi = inverse of tf, transforms constrained to unconstrained
      
    }

  init = tfi(as.matrix(init))
  
  ans = try(.Call("optim_ga1",
    f,
    init,
    as.integer(maxit),
    as.logical(minimize),
    as.double(tol),
    new.env(),
    PACKAGE = "uri"))

  if (is.double(ans$par))
    ans$par = tf(ans$par)

  ans
}








optim.grad1 = function(fn,
  init,
  lower = NULL,
  upper = NULL,
  maxit = 100,
  minimize = TRUE,
  npar     = length(init),
  tol      = sqrt(.Machine$double.eps),
  enforce.bounds = FALSE,
  ...)
{

  stopifnot(maxit >= 1 & is.double(init))
  
  if (!is.null(npar))
    stopifnot(npar == length(init))
  
  if((is.null(lower) & is.null(upper)) | !enforce.bounds) {
    
    tf = tfi = function(x) x
    f  = fn
  
  } else {
    
    tf  = transform.par(lower, upper, TRUE)
    tfi = transform.par(lower, upper, FALSE)
    f   = function(x) f(tf(x))
    
    ## tf transforms from unconstrained to constrained parameters
    ## tfi = inverse of tf, transforms back
    
  }

  init = tfi(init) # transform initial params to unconstrained

  ans = try(.Call("optim_gradient1",
    f,
    init,
    as.double(tol),
    as.integer(minimize),
    as.integer(maxit),
    new.env(),
    PACKAGE = "uri"))

  
  if (is.list(ans) & !is.null(ans$par) & !is.null(ans$value))
    ans$par = tfi(ans$par)

  ans
}






optim.grad2 = function(fn,
  init,
  lower    = NULL,
  upper    = NULL,
  maxit    = 100,
  minimize = TRUE,
  npar     = NULL,
  tol      = sqrt(.Machine$double.eps),
  enforce.bounds = FALSE,
  ...)

{
  init = as.matrix(init)
  
  stopifnot(maxit >= 1 & ncol(init) > 0)

  m = nrow(init)
  n = ncol(init)

  ans = list(par             = matrix(NA, m, n),
             value           = rep(NA, n),
             convergence      = rep(NA, n),
             gradient.lengths = rep(NA, n))

  for (k in 1:n) {
    
    ans1 = try(optim.grad1(fn = fn,
      init = init[, k],
      lower = lower,
      upper = upper,
      minimize = minimize,
      npar = npar,
      tol = tol,
      enforce.bounds = enforce.bounds,
      ...))

    if (is.list(ans1) & !is.null(ans1$val) & !is.null(ans1$par)) {
      
      ans$par[, k]      = ans1$par
      ans$value[k]      = ans1$val
      ans$convergence[k] = ans1$convergence
      ans$gradient.lengths[k] = ans1$gradient.length
      
    }
    
  }

  ans
}






transform.par = function(lower = NULL, upper = NULL,
                         to.constrained = FALSE) {
  
  ##*******************************************************************
  ##
  ##   Description: Returns function that transforms from constrained
  ##    to unconstrained parameters or vice versa.  to.constrained = FALSE
  ##    corresponds to transformation from constrained to unconstrained
  ##    parameters.
  ##
  ##*******************************************************************

  stopifnot(!is.null(lower) | !is.null(upper))

  if (is.null(lower)) lower = rep(-Inf, length(upper))
  if (is.null(upper)) upper = rep( Inf, length(lower))

  stopifnot(all(upper > lower, is.double(lower), is.double(upper),
                length(lower) == length(upper)))

  interval.length = upper - lower
  
  ind1 = which(-Inf < lower & upper < Inf)
  ind2 = which(-Inf < lower & upper == Inf)
  ind3 = which(lower == -Inf & upper < Inf)

  if (!to.constrained) {
    
    fn = function(x) {
      
      if (is.matrix(x) & length(x)) {
        ans = apply(x, 2, fn)
        dim(ans) = dim(x)
        dimnames(ans) = dimnames(x)
        ans
      } else if (is.list(x) & length(x)) {
        ans = lapply(x, fn)
        names(ans) = names(x)
        ans
      } else if (length(x)) {
        stopifnot(length(x) == length(lower))
        x[ind1] = qnorm(((x[ind1] - lower[ind1]) / interval.length[ind1]))
        x[ind2] =  log( x[ind2] - lower[ind2])
        x[ind3] = -log(-x[ind3] + upper[ind3])
        x
      } else x
    }

  } else {
    fn = function(x) {

      if (is.matrix(x) & length(x)) {
        ans = apply(x, 2, fn)
        dim(ans) = dim(x)
        dimnames(ans) = dimnames(x)
        ans
      } else if (is.list(x) & length(x)) {
        ans = lapply(x, fn)
        names(ans) = names(x)
        ans
      } else if (length(x)) {
        stopifnot(length(x) == length(lower))
        x[ind1] = pnorm(x[ind1]) * interval.length[ind1] + lower[ind1]
        x[ind2] = exp(x[ind2]) + lower[ind2]
        x[ind3] = upper[ind3] - exp(-x[ind3])
        x
      } else x
    }

  }

  fn
  
}



init.population = function(pop.size = 5,
  npar = NULL,
  lower = NULL,
  upper = NULL,
  init = NULL,
  mean = 1,
  sd = 1)
  ##
  ##  Description: Returns matrix whose columns consist of
  ##    all those from init the remaining from N(mean, sd), the latter
  ##    transformed to be row-wise in the intervals cbind(lower, upper).
  ##
{
  
## Last modified 10-Jan-2003.

  if (!is.null(init)) {
    
    init = as.matrix(init)
    npar = nrow(init)
    
  }

  if (is.null(npar))
    npar = ifelse(is.null(lower), length(upper), length(lower))
  
  if (is.null(lower))
    lower = rep(-Inf, npar)
  if (is.null(upper))
    upper = rep( Inf, npar)

  stopifnot(length(lower) == npar & length(upper) == npar)

  tf = transform.par(lower, upper, TRUE)
  
  if (is.matrix(init) & identical(ncol(init) > 0, TRUE)) {

    n        = ncol(init)
    pop.size = max(n, pop.size)
    par      = matrix(NA, npar, pop.size)
    
    par[,  seq(n)] = init
    par[, -seq(n)] = rnorm((pop.size - n) * npar)
    par[, -seq(n)] = tf(par[, -seq(n), drop = FALSE])
    
  } else if (is.null(init))
    par = tf(matrix(rnorm(pop.size * npar), npar, pop.size))
  

  par
}






