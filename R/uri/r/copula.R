



######################################################################
######################################################################
######################################################################
##
##  Following functions are general functions related to
##  empirical distributions & copulas.
##
######################################################################
######################################################################
######################################################################


emp.copula = function(x, y = NULL)
{
  ##
  ## Description: Returns the empirical copula function of x, y.
  ##
  
  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }
  
  stopifnot (is.matrix(x))

  S = nrow(x)
  M = ncol(x)

  x.os = apply(x, 2, sort)

  function(u) {

    stopifnot (all(0 <= u & u <= 1 & M == length(u)))
    u = trunc(u * S)
    mean(apply(x <= matrix(x.os[S * seq(0, M - 1) + u], S, M, TRUE), 1, prod))
    
  }
}







emp.cdf.fn = function(x) {
  ##
  ##  Description: Returns the empirical cdf function of x.
  ##

  x = x[!is.na(x)]
  
  function(y, na.rm = TRUE) {

    if (na.rm)
      y = y[!is.na(y)]
    
    .Call("emp_cdf_2", x, as.double(y))
  }
}



emp.cdf.eval = function(x, y, na.rm = TRUE) {
  ##
  ##  Description: Evaluates the empirical cdf of x
  ##    (with missing values removed always) on y
  ##    (where missing values of y are removed iff na.rm).
  ##
  if (na.rm) 
    y = y[!is.na(y)]

  .Call("emp_cdf_2", as.double(x[!is.na(x)]), as.double(y))
}







######################################################################
## The following functions extract the empirical uniform, normal
## and t-distributed random variables.

emp.unif = get.emp.unif = function(x, na.rm = FALSE)
{
  
  ans = .Call("empunif2", x)

  if (is.list(ans))
    ans = sapply(ans, as.double)

  if (na.rm) {

    if (is.list(ans))
      ans = sapply(ans, function(s) s[!is.na(s)])
    else if (is.matrix(ans))
      ans = apply(ans, 2, function(s) s[!is.na(s)])
    else if (is.double(ans))
      ans = ans[!is.na(ans)]
    
  }
  
  ans
}


emp.norm = get.emp.norm = function(x, na.rm = FALSE)
{

  ans = .Call("empnorm2", x)

  if (is.list(ans))
    ans = sapply(ans, as.double)

  if (na.rm) {

    if (is.list(ans))
      ans = sapply(ans, function(s) s[!is.na(s)])
    else if (is.matrix(ans))
      ans = apply(ans, 2, function(s) s[!is.na(s)])
    else if (is.double(ans))
      ans = ans[!is.na(ans)]
    
  }
  
  ans
  
}



emp.t = get.emp.t = function(x, df, na.rm = FALSE)
{

  ans = emp.unif(x, na.rm = na.rm)

  if (is.double(ans))
    ans = qt(ans, df = df)
  else if (is.list(ans))
    ans = sapply(ans, qt, df = df)
  else if (is.matrix(ans))
    ans = apply(ans, 2, qt, df = df)
  else
    ans = NA

  ans
           
}





######################################################################
######################################################################
######################################################################
##
##  Following functions calculate various dependence measures.
##
######################################################################
######################################################################
######################################################################




kendall.tau = function(x, y = NULL, use = "p") {

  if (!is.null(y))
    cor.test(x, y, method = "kendall", use = use)$estimate
  else {
    
    stopifnot (is.matrix(x) & ncol(x) >= 2)

    n = ncol(x)
    ans = diag(n)

    for (i in seq(1, n - 1)) {

      f = function(u) cor.test(x[, i], u, method = "kendall", use = use)$estimate
      
      ans[i, (i + 1):n] =
        ans[(i + 1):n, i] = apply(as.matrix(x[, (i + 1):n]), 2, f)

    }
    ans
  }
  
}



spearman.rho = function(x, y = NULL, use = "p") {

  if (!is.null(y))
    cor.test(x, y, method = "spearman", use = use)$estimate
  else {
    
    stopifnot (is.matrix(x) & ncol(x) >= 2)

    n   = ncol(x)
    ans = diag(n)

    
    for (i in seq(1, n - 1)) {

      f = function(u) cor.test(x[, i], u, method = "spearman")$estimate
      
      ans[i, (i + 1):n] =
        ans[(i + 1):n, i] = apply(as.matrix(x[, (i + 1):n]), 2, f)

    }
    ans
  }
}





rho.uri = function(x, y = NULL,
  method = c("pearson", "kendall", "spearman"),
  use    = "p") {
  ##
  ##  Description: Three methods of computing correlation.
  ##
  ##    pearson  - standard method.
  ##    kendall  - returns sin(0.5 * pi * kendall.tau(x, y)).
  ##    spearman - returns 2 * sin(pi / 6 * spearman.rho(x, y)).
  ##
  method = tolower(substr(method, 1, 1))

  if      ("p" == method)
    cor(x, y, use = use)
  else if ("k" == method)
    sin(0.5 * pi * kendall.tau(x, y, use = use))
  else if ("s" == method)
    2 * sin(pi / 6 * spearman.rho(x, y, use = use))
  else
    NA
}




######################################################################
######################### Plackett Copula ############################
######################################################################


C.plackett = function(u, v, phi) {

  stopifnot(length(u) == length(v))
  
  .C("C_plackett", as.double(u), as.double(v), as.double(phi),
     as.integer(length(u)), ans = double(length(u)))
}


get.C.plackett = function(phi) {
  ##
  ##   Description: Returns the Plackett copula function with
  ##    parameter phi.
  ##    Thus get.C.plackett(phi)(u, v) = C.plackett(u, v, phi).
  ##

  function(u, v) C.plackett(u, v, phi)
}



c.plackett = function(u, v, phi) {

  p = phi - 1
  
  ifelse(p, {
    phi * (1 + (u - 2 * u * v + v) * p) /
      ((1 + p * (u + v))^2 - 4 * u * v * phi * p)^1.5
  }, 1.0)
}



get.c.plackett = function(phi)
  function(u, v) c.plackett1(u, v, phi)


plackett.rho = function(phi)
  ifelse(1 == phi, 0, (phi + 1)/(phi - 1) - 2 * phi * log(phi) / (phi - 1)^2)




sim.plackett = function(n, phi, sim.type = 1) {

  if (1 == sim.type) {
    
    ans = .C("sim_plackett_1",
      u   = double(n),
      v   = double(n),
      phi = as.double(phi),
      n   = as.integer(n))

    cbind(ans$u, ans$v)
    
  }
  else {

    tmp = mvtunif.uri(n)
    u = tmp[,1]
    t = tmp[,2]

    d = (1 - 2 * t)^2
    e = ifelse(t <= 0.5, -1, 1)
    
    a0 = -(phi + 1)/(phi - 1)
    a1 = 2 * phi / (phi - 1) + 2 * phi * u
    b1 = -4 * u * phi
    b2 = 4 * u * phi + 4 * u^2 * phi * (phi - 1)
    
    c0 = a0^2 - d
    c1 = 2 * a0 * a1 - d * b1
    c2 = a1^2 - d * b2

    x = (-c1 + e * sqrt(c1^2 - 4 * c0 * c2)) / (2 * c0)
    cbind(u, (x - 1) / (phi - 1) - u)

    
  }
  
}




ll.plackett = function(u, v, phi) 
  ##
  ## Loglikelihood for bivariate Plackett copula
  ## where u, v are assumed to be uniform marginals
  ##
{
  stopifnot (length(u) == length(v))
  
  .C("ll_plackett",
     as.double(u),
     as.double(v),
     as.double(phi),
     as.integer(length(u)),
     ll = double(1))$ll
}





fit.plackett.cml = function(x, y,
  method = c("ga"),
  init   = NULL,
  maxit  = 50,
  lower  = 0,
  upper  = Inf,
  tol    = sqrt(.Machine$double.eps),
  enforce.bounds = FALSE,
  transf = TRUE)
{
  ##
  ## Fits Plackett copula using CML
  ##
  
  if (transf) {
    x = emp.unif(x)
    y = emp.unif(y)
  }
  
  ans = optim.uri(function(phi) ll.plackett(x, y, phi),
    method   = method,
    init     = init,
    npar     = 1,
    maxit    = maxit,
    lower    = lower,
    upper    = upper,
    minimize = FALSE,
    tol      = tol,
    enforce.bounds = enforce.bounds)

  ans$ll     = ans$values
  ans$values = NULL
  
  ans
}




######################################################################
################### Gaussian Copula ##################################  
######################################################################



ll.gaussian.mvt = function(x, y = NULL, rho, transf = FALSE) {

  if (!is.null(y)) {

    stopifnot(is.double(x))
    x = cbind(x, y)
    
  }

  if (!identical(all(abs(rho <= 1)), TRUE))
    -Inf
  else 
    .Call("ll_gaussian_mvt", x, rho, transf)

}


ll.gaussian.bvt = function(x, y, rho, transf = FALSE) {

  .Call("ll_gaussian_bvt", x, y, rho, transf)
  
}




fit.gaussian.cml = function(x, y = NULL, method = "fast",
  init = NULL, maxit = 50,  tol = 1e-4,
  enforce.bounds = TRUE, transf = TRUE, use = "p") {
  ##
  ## Fits Gaussian copula using CML
  ##
  if (is.null(y)) {

    x = as.matrix(x)
    stopifnot(ncol(x) > 1)

    n = ncol(x)

    if (2 == n) {
      ans = fit.gaussian.cml(x[,1], x[,2],
        method = method,
        init = init,
        maxit = maxit,
        tol  = tol,
        enforce.bounds = enforce.bounds,
        transf = transf)
    }
    else if ("fast" == method) {

      ans = list()

      if (transf)
        x = emp.norm(as.matrix(na.omit(as.data.frame(x))))

      ans$rho = cor(x)
      ans$ll  = ll.gaussian.mvt(x, rho = ans$rho, transf = FALSE)
      
    }
    else {
      
      if (transf)
        x = emp.norm(as.matrix(na.omit(as.data.frame(x))))
      
      if (is.null(init))
        init = cor(x)[upper.tri(diag(n))]

      npar  = 0.5 * n * (n - 1)
      
      ans = optim.uri(function(par) ll.gaussian.mvt(x,
        rho = {r = diag(n); r[upper.tri(r)] = r[lower.tri(r)] = par; r},
        transf = FALSE),
        method   = method, init  = init,           npar = npar,
        maxit    = maxit,  lower = rep(-1, npar), upper = rep(1, npar),
        minimize = FALSE,    tol = tol,  enforce.bounds = enforce.bounds)

      ans$ll     = ans$values
      ans$values = NULL
      ans$rho    = apply(ans$par, 2,
        function(par) {
          r = diag(n); r[upper.tri(r)] = r[lower.tri(r)] = par; r
        })
      
    }
    
  }
  else { # if y is not null

    if (transf) {
      x = emp.norm(x)
      y = emp.norm(y)
    }
    
    if ("fast" == method) {

      ans = list()
      ans$rho = cor(x, y)
      ans$ll  = ll.gaussian.bvt(x, y, ans$rho, transf = FALSE)
      
    }
    else {
      
      if (is.null(init))
        init = rho.uri(x, y, "k")

      ans = optim.uri(function(rho)
        ll.gaussian.bvt(x, y, rho, transf = FALSE),
        method = method, init = init, npar = 1,
        maxit = maxit, lower = -1, upper = 1,
        minimize = FALSE, tol = tol, enforce.bounds = enforce.bounds)

      ans$ll     = ans$values
      ans$values = NULL
      ans$rho    = ans$par
    }
  }
  
  ans
}





######################################################################
##################### Student Copula #################################
######################################################################



student.ll.1 = function(x, y = NULL, rho, df, transf = FALSE) {
  ## Loglikelihood for student copula
  
  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }

  stopifnot (is.matrix(x))
  
  S = nrow(x)
  M = ncol(x)
  
  a = pseudo.inverse(rho)
  g = function(v) as.double(t(v) %*% a %*% v)

  if (transf) x = apply(x, 2, emp.unif)

  x = apply(x, 2, qt, df = df)
  
  term1 = sum(-0.5 * (df + M) * log(1 + apply(x, 1, g) / df) +
    0.5 * (df + 1) * log(apply(1 + x^2 / df, 1, prod)))

  term2 = S * (-0.5 * log(det(rho)) + lgamma(0.5 * (df + M))
    + (M - 1) * lgamma(0.5 * df) - M * lgamma(0.5 * (df + 1)))

  term1 + term2

}




student.ll.2 = function(x, y, rho, df, transf = FALSE) {
  ## Canonical loglikelihood for bivariate student copula
  
  stopifnot(is.double(x) & is.double(y) & length(x) == length(y))

  S = length(x)

  if (transf) {
    x = emp.unif(x)
    y = emp.unif(y)
  }

  x = qt(x, df = df)
  y = qt(y, df = df)
  
  term1 = sum(
    - 0.5 * (df + 2) * log(1 + (x^2 - 2*rho*x*y + y^2) / ((1 - rho^2) * df))
    + 0.5 * (df + 1) * log((1 + x^2 / df) * (1 + y^2 / df)))

  term2 = S * (-0.5 * log(1 - rho^2) + lgamma(0.5 * (df + 2))
    + lgamma(0.5 * df) - 2 * lgamma(0.5 * (df + 1)))

  term1 + term2

}


sim.student = function(n, rho, df) {  
  z = rmvnorm(n, sigma = rho)
  s = rchisq(n, df)
  pt(sqrt(df / s) * z, df)
}



gaussian.cml.1 = function(x, y = NULL) {
  
  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }

  stopifnot(is.matrix(x))

  z = apply(x, 2, emp.norm)

  cor(z)
}



student.cml.1 = function(x, y = NULL, df.limit = 1e5) {
  
  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }

  M = ncol(x)

  to.rho = function(par) {
    rho = diag(M)
    rho[upper.tri(rho)] = par
    rho[lower.tri(rho)] = t(rho)[lower.tri(rho)]
    rho
  }
  
  obj.fn = function(par) {
    stopifnot (length(par) == M * (M - 1) / 2 + 1)
    student.ll.1(x, rho = to.rho(par[-1]), df = par[1])
  }

  lower = c(0,        rep(-1, M * (M - 1) / 2))
  upper = c(df.limit, rep( 1, M * (M - 1) / 2))

  ans = optim.uri(obj.fn, init = c(30, rep(0, M * (M - 1) / 2)),
            lower = lower, upper = upper,
            min = FALSE, method = "ga")
  
  ans = optim.uri(obj.fn, init = ans$par,
            lower = lower, upper = upper,
            min = FALSE, method = "N")
  
  ans = optim.uri(obj.fn, init = ans$par,
            lower = lower, upper = upper,
            min = FALSE, method = "ga")

  ans$rho = list()
  ans$df = 0
  for (k in seq(length(ans$par))) {
    ans$df[k] = ans$par[[k]][1]
    ans$rho[[k]] = to.rho(ans$par[[k]][-1])
  }
  ans$par = NULL
  ans$ll = ans$val
  ans$val = NULL
  ans$rho.best = ans$rho[[1]]
  ans$ll.best = ans$ll[1]

  ans
}



student.cml.2 = function(x, y, df.limit = 1e4) {
  
  stopifnot(is.double(x) & is.double(y) & length(x) == length(y))
  
  obj.fn = function(par) {
    stopifnot (2 == length(par))
    student.ll.2(x, y, rho = par[1], df = par[2])
  }

  lower = c(-1, 0)
  upper = c( 1, df.limit)

  ans = optim.uri(obj.fn, init = c(0, 30),
            lower = lower, upper = upper,
            min = FALSE, method = "ga")
  
  ans = optim.uri(obj.fn, init = ans$par,
            lower = lower, upper = upper,
            min = FALSE, method = "N")
  
  ans = optim.uri(obj.fn, init = ans$par,
            lower = lower, upper = upper,
            min = FALSE, method = "ga")

  ans$df = ans$rho = 0
  for (k in seq(length(ans$par))) {
    ans$df[k] = ans$par[[k]][2]
    ans$rho[k] = ans$par[[k]][1]
  }
  
  ans$par = NULL
  ans$ll = ans$val
  ans$val = NULL
  ans$rho.best = ans$rho[1]
  ans$ll.best = ans$ll[1]

  ans
}



student.cml.3 = function(x, y = NULL, df = 1e2, tol = 1e-3, maxit = 50) {
  ## Doesn't seem to work

  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }
  
  stopifnot (is.matrix(x))
  
  z = apply(x, 2, emp.norm)
  rho1 = cor(z)

  S = nrow(x)
  M = ncol(x)
  it.number = 0
  
  repeat {

    a = solve(rho1)

    b = apply(z, 1, function(r) as.double(t(r) %*% a %*% r)) / df
    
    rho2 = (df + M) / (S * df) * t(z) %*% (z / b)

    rho2 = rho2 / (abs(matrix(diag(rho2), M, M, byrow = TRUE)) *
      abs(matrix(diag(rho2), M, M)))
      
    if (sum((rho2 - rho1)^2) < tol | it.number > maxit)
      break
    else {
      it.number = it.number + 1
      rho1 = rho2
    }
  }

  rho2
  
}




######################################################################
########################## Gumbel Copula #############################
######################################################################



C.gumbel = function(u, v, a)
  exp(-((-log(u))^a + (-log(v))^a)^(1/a))

c.gumbel = function(u, v, a) {

  uu = -log(u)
  vv = -log(v)

  (uu * vv)^(a - 1) / (u * v) * (uu^a + vv^a)^(1/a - 2) *
    exp(-(uu^a + vv^a)^(1/a)) * ((uu^a + vv^a)^(1/a) + a - 1)
}




gumbel.ll.1 = function(x, y, a, transf = FALSE) {

  stopifnot (all(is.double(x), is.double(y), is.double(a),
                 length(x) == length(y), a >= 1))

  n = length(x);

  if (transf) {
    x = .C("get_emp_unif", as.double(x), as.integer(n), ans = double(n))$ans
    y = .C("get_emp_unif", as.double(y), as.integer(n), ans = double(n))$ans
  }
  
  .C("gumbel_ml", as.integer(n), as.double(x), as.double(y),
     as.double(a[1]), loglik = double(1))$loglik
}




gumbel.ll.2 = function(u, v, a) {

  stopifnot (all(is.double(u), is.double(v), is.double(a),
                 length(u) == length(v), a >= 1, 0 <= u,
                 0 <= v, u <= 1, v <= 1))

  n = length(u);

  .C("gumbel_ml", as.integer(n), as.double(u), as.double(v),
     as.double(a[1]), loglik = double(1))$loglik
}




######################################################################
######################### Old Functions ##############################
######################################################################


sim.plackett.1.old = function(n, phi) {
  ##
  ##   Description: Returns a matrix of dimension n x 2 simulated
  ##    from the Plackett copula.
  ##

  u = runif(n)
  t = runif(n)

  d = (1 - 2 * t)^2
  e = ifelse(t <= 0.5, -1, 1)
  
  a0 = -(phi + 1)/(phi - 1)
  a1 = 2 * phi / (phi - 1) + 2 * phi * u
  b1 = -4 * u * phi
  b2 =  4 * u * phi + 4 * u^2 * phi * (phi - 1)
  
  c0 = a0^2 - d
  c1 = 2 * a0 * a1 - d * b1
  c2 = a1^2 - d * b2

  x = (-c1 + e * sqrt(c1^2 - 4 * c0 * c2)) / (2 * c0)
  cbind(u, (x - 1) / (phi - 1) - u)
  
}





emp.unif.old = function(x) emp.cdf.old(x)(x)
emp.norm.old = function(x) qnorm(emp.cdf.old(x)(x))
emp.t.old    = function(x, df) qt(emp.cdf.old(x)(x), df = df)




emp.cdf.old = function(x, na.rm = TRUE) {
  ##
  ## Description: Returns the empirical cdf corresponding to the
  ##    vector x.
  ##

  if (na.rm) x = x[!is.na(x)]
  x = as.double(x)

  function(y)
    apply(outer(x, y, "<"), 2, mean) + 0.5 / length(x)
}




ll.plackett.old = function(u, v, phi) {
  ##
  ## Loglikelihood for bivariate Plackett copula
  ## where u, v are assumed to be uniform marginals
  ##
  ifelse(phi - 1, {
    a   = log(phi) + log(1 + (u - 2 * u * v + v) * (phi - 1))
    b   = 1.5 * log((1 + (phi - 1) * (u + v))^2 - 4 * u * v * phi * (phi - 1))
    ans = sum(a - b)
    ifelse(is.finite(ans), ans, -Inf)
  }, 0.0)
}





C.plackett.old = function(u, v, phi) {
  ##
  ##   Description: Evaluates the Plackett copula function with
  ##    parameter phi on u, v and returns the result.
  ##
  
  p   = phi - 1
  
  ifelse(p, {
    t1 = 1 + p * (u + v);
    t2 = sqrt((1 + p * (u + v))^2 - 4 * u * v * phi* p);
    0.5 / p * (t1 - t2)
  }, u * v)
  
}





ll.gaussian.mvt.old = function(x, y = NULL, rho, transf = FALSE) {
## Canonical loglikelihood for multivariate Gaussian copula
  if (!is.null(y)) {
    stopifnot (is.double(x))
    x = cbind(x, y)
  }

  if (!identical(all(abs(rho) <= 1), TRUE))
    -Inf
  else {
    
    S = nrow(x)
    M = ncol(x)
  
    a = try(solve(rho) - diag(M))

    if (is.matrix(a)) {
    
      g = function(v) as.double(t(v) %*% a %*% v)
      
      if (transf)
        x = apply(x, 2, emp.norm)
      
      -0.5 * S * log(det(rho)) - 0.5 * sum(apply(x, 1, g))
    
    }
    else -Inf
    
  }
}



ll.gaussian.bvt.old = function(x, y, rho, transf = FALSE) {
  ## Canonical loglikelihood for bivariate Gaussian copula

  if (!identical(abs(rho) <= 1, TRUE))
    -Inf
  else {
    
    stopifnot(is.double(x) & is.double(y) & length(x) == length(y))

    S = length(x)

    if (transf) {
      x = emp.norm(x)
      y = emp.norm(y)
    }
        
    -0.5*S*log(1 - rho^2) - 0.5*rho/(1 - rho^2)*sum(rho*x^2 - 2*x*y + rho*y^2)

  }
}

