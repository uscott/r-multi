


evalNgarch11 = function(x, pars, llonly = FALSE)
    ##
    ##  Description: Returns list including loglikelihood,
    ##    conditional variances and residuals for the given
    ##    time series x and parameters.
    ##
{
    .Call("ngarch11", x, pars, llonly, PACKAGE = "uri")
}


fitNgarch11 = function(x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, fitInit = FALSE, option = NULL)
{
    .Call( "fit_ngarch11", x, init, fitInit,
          popSize, tol,
          stopLags, minit, maxit, option,
          PACKAGE = "uri")
}

## modInfo should correspond to the return value of fitNgarch11

bootfitNgarch11 = function(
    modInfo, numBoots = 20, popSize = 5, stopLags = 5,
    minit = 5, maxit = 100,
    tol = 1e-5, fitInit = FALSE, option = NULL )
{
    z = modInfo$res
    len = length( modInfo$res )
    pars = modInfo$par[1:5,1]
      
    for ( i in 1:numBoots )
    {
        indices = sample( 1:len, len, rep = TRUE )
        x = simNgarch11( len, paths = 1, pars, x0 = NULL, h0 = NULL,
            z = z[ indices ] )$x
        modInfo2 = fitNgarch11( x, init = pars, popSize = popSize,
            stopLags = stopLags, minit = minit, maxit = maxit, tol = tol,
            fitInit = fitInit, option = option )
    }
}



wfitNgarch11 = function(x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, fitInit = FALSE, option = NULL)
{
    .Call("fit_ngarch11_w1", x, init, fitInit,
          popSize, tol,
          stopLags, minit, maxit, option,
          PACKAGE = "uri")
}



simNgarch11 = function(n, paths = 1, pars, x0 = NULL, h0 = NULL, z = NULL)
{
    .Call("sim_ngarch11", n, paths, pars, x0, h0, z, PACKAGE = "uri")
}



ngarch11prices = function(exDts, modelInfo,
    paths = 1e4, r = 0, q = 0,
    rsd = 0, ## Not implemented
    qsd = 0, ## Not implemented
    div = NULL, boot = FALSE, risk.neutral = TRUE)
    ## Function for extracting terminal distribution.
    ## modelInfo contains all the relevant info from a fitted
    ## NGARCH(1, 1) model.  div is a named vector of
    ## dividends where the names are the dividend dates
    ## (div argument currently not implemented).
    ## q is a cts dividend yield.  qsd is the standard
    ## deviation of q as a truncated normal RV.
{
    
    S0 = 1
  
    if (any(is.null(modelInfo$par),
            is.null(modelInfo$res),
            boot & is.null(modelInfo$bootPar),
            is.null(modelInfo$freq)))
        stop ("incomplete NGARCH(1, 1) model info given")

    ## Get time to expiration
    t0   = strptime(as.char(Sys.time()), format = "%Y-%m-%d %H:%M:%S")
    tau  = as.double(as.POSIXlt(exDts) - t0)
    names(tau) = names(exDts)
    tau  = tau[tau > 0]
    freq = modelInfo$freq # Frequency in days
    tau  = sort(tau) / freq

    tauInt  = trunc(tau)
    tauFrac = tau - tauInt

    m = length(tau)
    n = ceiling(max(tau))

    ## Simulate terminal prices
    N = ceiling(paths / nrow(modelInfo$bootPar))

    if (boot)
    {
        lambda = rep(modelInfo$bootPar$lambda, N)[1:paths]
        a0     = rep(modelInfo$bootPar$a0,     N)[1:paths] 
        a1     = rep(modelInfo$bootPar$a1,     N)[1:paths] 
        b1     = rep(modelInfo$bootPar$b1,     N)[1:paths] 
        gamma  = rep(modelInfo$bootPar$gamma,  N)[1:paths]
        h0     = rep(modelInfo$bootPar$h,      N)[1:paths]
    }
    else
    {
        lambda = modelInfo$par["lambda", 1]
        a0     = modelInfo$par["a0",     1]
        a1     = modelInfo$par["a1",     1]
        b1     = modelInfo$par["b1",     1]
        gamma  = modelInfo$par["gamma",  1]
        h0     = last(modelInfo$h)
    }

    x0 = last(modelInfo$x)
  
    if (risk.neutral)
    {
        gamma  = lambda + gamma
        lambda = 0
    }

    pars = list(lambda=lambda,a0=a0,a1=a1,b1=b1,gamma=gamma)
  
    rsim = simNgarch11(n = n, paths = paths,
        par = pars, z  = modelInfo$res, x0 = x0, h0 = h0)$x
  
    rsim = rsim[, !is.na(rsim[n, ]), drop = FALSE]

  
    r1  = apply(rsim, 2, cumsum)[tauInt, , drop = FALSE]
    r2  = rsim[tauInt + 1, , drop = FALSE]
    mr2 = apply(r2, 1, mean)
    r2  = (r2 - mr2) * sqrt(tauFrac) + mr2 * tauFrac
  
    S   = exp(r1 + r2)

    if (risk.neutral)
        S = S / apply(S, 1, mean) * exp(tau * (r - q) / 365 * freq)
    else
        S = S * exp(-tau * q / 365 * freq)
 
    gc()

    stopifnot (nrow(S) == m)
   
    rownames(S) = names(tau)

    list(S = t(S), tau = tau * freq, freq = freq, r = r, q = q, t0 = t0)
}



ngarch11vols = function(distInfo,  S0 =  1, strikes)
    ## Function for pricing options.
    ## distInfo is a list of length 2 containing
    ## the physical terminal distribution and the corresponding
    ## times to expiration
    ## NGARCH(1, 1) model.
    ##
{
    S    = t(distInfo$S) * S0
    tau  = distInfo$tau
    freq = distInfo$freq
    t0   = distInfo$t0
    r    = distInfo$r
    q    = distInfo$q
  
    stopifnot (is.matrix(S)&is.vector(tau)&is.double(freq)&is.double(S0))
    stopifnot (all(tau > 0, freq > 0, S0 > 0))
    stopifnot (all(length(strikes), strikes > 0))

    tau     = sort(tau)
    m       = length(tau)
    n       = ceiling(max(tau))
    tauInt  = trunc(tau)
    tauFrac = tau - tauInt

    if (1 == length(r))
        r = rep(r, m)
    if (1 == length(q))
        q = rep(q, m)
  
    rmq = r - q

    gc()
  
    if (is.null(dim(S))) S = t(S)
  
    stopifnot (nrow(S) == m)
  
    mu  = apply(S, 1, mean) / S0
  
    K = strikes; rm(strikes)
  
    if (is.vector(K))
        K = matrix(K, m, length(K), byrow = TRUE)
    else if (is.matrix(K) & identical(1 == nrow(K), TRUE))
        K = matrix(as.vector(K), m, ncol(K), byrow = TRUE)
  
    K = t(apply(K, 1, sort))
    n = ncol(K)

    stopifnot (is.matrix(K))
  
    calls  = matrix(NA, m, n)
    puts   = matrix(NA, m, n)
    prices = matrix(NA, m, n)
    deltas = matrix(NA, m, n)
    gammas = matrix(NA, m, n)
    vols   = matrix(NA, m, n)
  
    optTypes    = ifelse(K < S0, "p", "c")
    cix         = (K >= S0)
    pix         = (K <  S0)
    payoffSigns = ifelse(K >= S0, 1, -1)

    u = 1e-2

    tau = tau / 365.0
  
    for (i in seq(m))
    {
        discFac    = exp(-r[i] * tau[i])
        sgn        = matrix(payoffSigns[i,], ncol(S), n, byrow = TRUE)
        
        payoff     = pmax(sgn * outer(S[i,], K[i,], "-"), 0)
    
        prices[i,] = discFac * apply(payoff, 2, mean)
    
        payoff2    = pmax(sgn * outer((1 + u) * S[i,], K[i,], "-"), 0)
        payoff3    = pmax(sgn * outer((1 - u) * S[i,], K[i,], "-"), 0)
        
        prices2    = discFac * apply(payoff2, 2, mean)
        prices3    = discFac * apply(payoff3, 2, mean)
    
        deltas[i,] = 100 * (prices2 - prices3) / (2 * u * S0) 
        gammas[i,] = 100 * (prices2 - 2 * prices[i, ] + prices3) / (u * S0)^2

        deltas[i,] = round(deltas[i,], 2)
        gammas[i,] = round(gammas[i,], 2)

        vols[i,]   = impvol(365 * tau[i], S0, K[i, ],
                price = prices[i, ], r = r[i], q = q[i],
                opt = optTypes[i, ], tol = 1e-7, maxit = 1e2)

        calls[i, ] = prices[i, ]
        puts [i, ] = prices[i, ]

        parity = exp(-q[i]*tau[i])*S0 - discFac*K[i,]
    
        calls[i,pix[i,]] = calls[i,pix[i,]] + parity[pix[i,]]
        puts [i,cix[i,]] = puts [i,cix[i,]] - parity[cix[i,]]
        
    }
    
    if (!is.null(names(tau)))
        exp.names = names(tau)
    else if (!is.null(rownames(K)))
        exp.names = rownames(K)
    else
        exp.names = paste("expiration.date", seq(m), sep = "" )

    strike.names = paste(K, optTypes, sep = "")
    dim(strike.names) = dim(K)
    
    rownames(prices)       = rownames(calls)    = rownames(puts) = exp.names
    rownames(vols)         = rownames(optTypes) = rownames(pix)  = exp.names
    rownames(deltas)       = rownames(gammas)   = rownames(S)    = exp.names
    names(mu)              = rownames(K)        = names(q)       = exp.names
    rownames(strike.names) = rownames(cix)      = names(r)       = exp.names
  
    ans = list()
    for (i in exp.names)
    {
        ans[[i]] = matrix(NA, 8, n)
    
        rownames(ans[[i]]) =
            c("calls", "puts", "calls.ratio",
              "puts.ratio", "vols", "call.delta", "put.delta", "gamma")

        ans[[i]][1, ] = round(calls[i, ], 3)
        ans[[i]][2, ] = round(puts [i, ], 3)

        atmix = which.min(abs(K[i, ] - S0))
    
        patm = puts [i, atmix]
        catm = calls[i, atmix]

        ans[[i]][3, ] = round(catm / calls[i, ], 2)
        ans[[i]][4, ] = round(patm / puts [i, ], 2)
        ans[[i]][5, ] = round(vols[i, ], 2)
    
        ans[[i]][6, ] = deltas[i, ] +
            100 * ifelse(pix[i, ], rep(exp(-q[i] * tau[i]), n), rep(0, n))
    
        ans[[i]][7, ] = deltas[i, ] -
            100 * ifelse(cix[i, ], rep(exp(-q[i] * tau[i]), n), rep(0, n))
    
        ans[[i]][8, ] = gammas[i, ]
    
        ans[[i]] = as.data.frame(t(ans[[i]]))
        rownames(ans[[i]]) = K[i,]
      
  }

    ans$vert.sprds = list()
    ans$cs = ans$ps = list()
  
    for (i in exp.names)
    {
        
        ans$vert.sprds[[i]] = round(outer(prices[i,], prices[i,], "-"), 3)
        ans$cs[[i]] = round(outer(calls[i,],calls[i,],"-"), 3)
        ans$ps[[i]] = round(outer(puts[i,],puts[i,],"-"), 3)
    
        dimnames(ans$vert.sprds[[i]]) =
            list(paste("+", strike.names[i,], sep = ""),
                 paste("-", strike.names[i,], sep = ""))

        dimnames(ans$cs[[i]]) = dimnames(ans$ps[[i]]) =
            list(paste("+", K[i,], sep=""),
                 paste("-", K[i,], sep=""))
    
    }
  
    ans$prices       = drop(prices)
    ans$vols         = drop(vols)
    ans$deltas       = drop(deltas)
    ans$mu           = drop(mu)
    ans$disc.rates   = drop(r)
    ans$div.yields   = drop(q)
    ans$tau          = drop(tau)
    ans$strikes      = drop(K)
    ans$distribution = drop(S)
    ans$t0           = t0
    ans$S0           = S0
  
    ans
}



######################################################################
######################################################################
######################################################################
##
## This section has functions for fitting the AR(1)-NGARCH(1, 1)
## model
##
##  x[t] = mu + phi * x[t - 1] + sqrt(h[t]) * (lambda + z[t])
##  h[t] = a0 + h[t - 1] * (b1 + a1 * (z[t - 1] - gamma)^2)
##
######################################################################
######################################################################
######################################################################



ar1ngarch11 =
  function(x,
           mu = par[1],
           lambda = par[2],
           phi = par[3],
           a0 = par[4],
           a1 = par[5],
           b1 = par[6],
           gamma = par[7],
           xmin1 = par[8],
           h0 = par[9],
           par = NULL,
           ll.only = FALSE)
  ##
  ##  Description: Returns list including loglikelihood,
  ##    conditional variances and residuals for the given
  ##    time series x and parameters.
  ##

{
    n = length(x)

    if (ll.only)
        h = res = numeric(0)
    else
    {
        h   = double(n)
        res = double(n)
    }
    
    .C( "ar1ng11",
       n       = as.integer(n),
       x       = as.double(x),
       mu      = as.double(mu),
       lambda  = as.double(lambda),
       phi     = as.double(phi),
       a0      = as.double(a0),
       a1      = as.double(a1),
       b1      = as.double(b1),
       gamma   = as.double(gamma),
       h0      = as.double(h0),
       xmin1   = as.double(xmin1),
       h       = h,
       res     = res,      
       ll      = double(1),
       NAOK    = TRUE,
       PACKAGE = "uri")
}



fit.ar1ngarch11 = function( x,
    method = c("ga", "N", "ga"),
    init   = NULL,
    maxit  = 5,
    lower  = c(-Inf, -Inf, -1,   0,   0,   0, -Inf,   0, -Inf),
    upper  = c( Inf,  Inf,  1, Inf, Inf, Inf,  Inf, Inf,  Inf),
    tol    = sqrt(.Machine$double.eps),
    enforce.bounds = FALSE,
    fit.init.state = FALSE, ...)
    ##
    ##  Description: Fits AR(1)-NGARCH(1, 1) model to time
    ##    series x:
    ##
    ##    x[t] = mu + phi * x[t-1] + lambda * sqrt(h[t]) - 0.5 * h[t]
    ##              + sqrt(h[t])
    ##
{
    n = length(x)
  
    npar = ifelse(fit.init.state, 9, 7)

    ## Define the loglikelihood function
    llfn = function(par)
    {
        .C("ar1ng11",
           n       = as.integer(n),
           x       = as.double(x),
           mu      = par[1],
           lambda  = par[2],
           phi     = par[3],
           a0      = par[4],
           a1      = par[5],
           b1      = par[6],
           gamma   = par[7],
           h0      = par[8],
           xmin1   = par[9],
           h       = numeric(0),
           res     = numeric(0),
           ll      = double(1),
           NAOK    = TRUE,
           PACKAGE = "uri")$ll
    }

    ## Prepare starting parameters
    if (is.null(init))
    {
        require(tseries)
        init = runif(npar, 0, 1e-8)

        names(init) = c("mu", "lambda", "phi",
                 "a0" , "a1", "b1",
                 "gamma", "h0", "xmin1")[1:npar]

        init[c("a0", "a1", "b1")] =
            garch(x, trace = FALSE)$coef[c("a0", "a1", "b1")]
    }
    
    init = as.matrix(init)

    ## Fit the model
    fit = optim.uri(llfn,
        method   = method,
        init     = init,
        npar     = npar,
        maxit    = maxit,
        lower    = lower[1:npar],
        upper    = upper[1:npar],
        minimize = FALSE,
        tol      = tol,
        enforce.bounds = enforce.bounds,
        ...)

    ## Work with outputs a little
    if (!is.null(fit$par))
    {
    
        rownames(fit$par) =
            c("mu", "lambda", "phi", "a0", "a1", "b1", "gamma", "h0", "xmin1")[1:npar]

        j          = which.max(fit$value)
        fit$x      = x
        fit$ll     = fit$value
        fit$value  = NULL

        tmp     = ar1ngarch11(x, par = fit$par[, j])
        fit$h   = tmp$h
        fit$res = tmp$res
        fit$aic = aic(fit$ll, npar)
        fit$bic = bic(fit$ll, npar, n)

        names(fit$h) = names(fit$x)
    }

    fit
}



sim.ar1ngarch11 =
    function(n,
             paths = 1,
             mu       = 0,
             lambda   = 0,
             phi      = 0,
             a0       = 0,
             a1       = 0,
             b1       = 0,
             gamma    = 0,
             h0       = NULL,
             x.min1   = NULL,
             par      = NULL,
             V        = NULL,
             z        = numeric(0))
{
    if      (is.vector(par))
    {
        mu      = par[1]
        lambda  = par[2]
        phi     = par[3]
        a0      = par[4]
        a1      = par[5]
        b1      = par[6]
        gamma   = par[7]
        
        if (missing(h0))     h0      = par[8]
        if (missing(x.min1)) x.min1  = par[9]
    }
    else if (is.list(par))
    {
        mu      = par$mu
        lambda  = par$lambda
        phi     = par$phi
        a0      = par$a0
        a1      = par$a1
        b1      = par$b1
        gamma   = par$gamma
        if (missing(h0))     h0      = par$h0
        if (missing(x.min1)) x.min1  = par$x.min1
    }
    else if (is.matrix(par))
    {
        mu      = par[, 1]
        lambda  = par[, 2]
        phi     = par[, 3]
        a0      = par[, 4]
        a1      = par[, 5]
        b1      = par[, 6]
        gamma   = par[, 7]
        if (missing(h0))     h0      = par[, 8]
        if (missing(x.min1)) x.min1  = par[, 9]
    }
    else if (!is.null(par))
        stop ("don't know what to do with argument par")


    if (!is.null(V))
        a0 = V * (1 - b1 - a1 * (1 + gamma^2))
    
    .Call("sim_ar1ng11",
          as.integer(n),
          as.integer(paths),
          as.double(mu),
          as.double(lambda),
          as.double(phi),
          as.double(a0),
          as.double(a1),
          as.double(b1),
          as.double(gamma),
          as.double(h0),
          as.double(x.min1),
          as.double(z))
}



ar1ngarch11.prices = function(exDts, modelInfo,
    paths        = 1e4,
    r            = 0,
    q            = 0,
    rsd          = 0, ## Not implemented
    qsd          = 0, ## Not implemented
    div          = NULL,
    boot         = TRUE,
    risk.neutral = TRUE)
    ## Function for extracting terminal distribution.
    ## modelInfo contains all the relevant info from a fitted
    ## AR(1)-NGARCH(1, 1) model.  div is a named vector of
    ## dividends where the names are the dividend dates
    ## (div argument currently not implemented).
    ## q is a cts dividend yield.  qsd is the standard
    ## deviation of q as a truncated normal RV.
{
    S0 = 1
  
    if (any(is.null(modelInfo$par),
            is.null(modelInfo$res),
            boot & is.null(modelInfo$bootPar),
            is.null(modelInfo$freq)))
        stop ("incomplete AR(1)-NGARCH(1, 1) model info given")

    ## Get time to expiration
    t0   = strptime(as.char(Sys.time()), format = "%Y-%m-%d %H:%M:%S")
    tau  = as.double(as.POSIXlt(exDts) - t0)
    names(tau) = names(exDts)
    tau  = tau[tau > 0]
    freq = modelInfo$freq # Frequency in days
    tau  = sort(tau) / freq

    tauInt  = trunc(tau)
    tauFrac = tau - tauInt

    m = length(tau)
    n = ceiling(max(tau))

    ## Simulate terminal prices
    N = ceiling(paths / nrow(modelInfo$bootPar))

    if (boot)
    {
        mu     = rep(modelInfo$bootPar$mu,     N)[1:paths]
        lambda = rep(modelInfo$bootPar$lambda, N)[1:paths]
        phi    = rep(modelInfo$bootPar$phi,    N)[1:paths]
        a0     = rep(modelInfo$bootPar$a0,     N)[1:paths] 
        a1     = rep(modelInfo$bootPar$a1,     N)[1:paths] 
        b1     = rep(modelInfo$bootPar$b1,     N)[1:paths] 
        gamma  = rep(modelInfo$bootPar$gamma,  N)[1:paths] 
        h0     = rep(modelInfo$bootPar$h,      N)[1:paths]
        x.min1 = rep(last(modelInfo$x), paths)
    }
    else
    {
        mu     = modelInfo$par["mu",     1]
        lambda = modelInfo$par["lambda", 1]
        phi    = modelInfo$par["phi",    1]
        a0     = modelInfo$par["a0",     1]
        a1     = modelInfo$par["a1",     1]
        b1     = modelInfo$par["b1",     1]
        gamma  = modelInfo$par["gamma",  1]
        h0     = last(modelInfo$h)
        x.min1 = last(modelInfo$x)
    }

    if (risk.neutral)
    {
        gamma  = lambda + gamma
        lambda = lambda - lambda
        mu     = mu - mu
    }
  
    rsim = sim.ar1ngarch11(n = n,
        paths  = paths,
        mu     = mu,
        lambda = lambda,
        phi    = phi,
        a0     = a0,
        a1     = a1,
        b1     = b1,
        z      = modelInfo$res,
        h0     = h0,
        x.min1 = x.min1)$x

    rsim = rsim[, !is.na(rsim[n, ]), drop = FALSE]
  
    r1  = apply(rsim, 2, cumsum)[tauInt, , drop = FALSE]
    r2  = rsim[tauInt + 1, , drop = FALSE]
    mr2 = apply(r2, 1, mean)
    r2  = (r2 - mr2) * sqrt(tauFrac) + mr2 * tauFrac
  
    S   = exp(r1 + r2)

    if (risk.neutral)
        S = S / apply(S, 1, mean) * exp(tau * (r - q) / 365 * freq)
    else
        S = S * exp(-tau * q / 365 * freq)
  
    gc()
  
    stopifnot (nrow(S) == m)
    
    rownames(S) = names(tau)

    list(S = t(S), tau = tau * freq, freq = freq, r = r, q = q, t0 = t0)
}



ar1ngarch11.vols = function(distInfo,  S0 =  1, strikes)
    ## Function for pricing options.
    ## distInfo is a list of length 2 containing
    ## the physical terminal distribution and the corresponding
    ## times to expiration
    ## AR(1)-NGARCH(1, 1) model.
    ##
{
    S    = t(distInfo$S) * S0
    tau  = distInfo$tau
    freq = distInfo$freq
    t0   = distInfo$t0
    r    = distInfo$r
    q    = distInfo$q
    
    stopifnot (all(is.matrix(S),    is.vector(tau),
                   is.double(freq), is.double(S0)))
    stopifnot (all(tau > 0, freq > 0, S0 > 0))
    stopifnot (all(length(strikes), strikes > 0))

    tau     = sort(tau)
    m       = length(tau)
    n       = ceiling(max(tau))
    tauInt  = trunc(tau)
    tauFrac = tau - tauInt

    if (1 == length(r))
        r = rep(r, m)
    if (1 == length(q))
        q = rep(q, m)
  
    rmq = r - q

    gc()
  
    if (is.null(dim(S))) S = t(S)
  
    stopifnot (nrow(S) == m)
  
    mu  = apply(S, 1, mean) / S0
  
    K = strikes; rm(strikes)
  
    if (is.vector(K))
        K = matrix(K, m, length(K), byrow = TRUE)
    else if (is.matrix(K) & identical(1 == nrow(K), TRUE))
        K = matrix(as.vector(K), m, ncol(K), byrow = TRUE)
  
    K = t(apply(K, 1, sort))
    n = ncol(K)

    stopifnot (is.matrix(K))
  
    calls  = matrix(NA, m, n)
    puts   = matrix(NA, m, n)
    prices = matrix(NA, m, n)
    deltas = matrix(NA, m, n)
    gammas = matrix(NA, m, n)
    vols   = matrix(NA, m, n)
  
    optTypes    = ifelse(K < S0, "p", "c")
    cix          = (K >= S0)
    pix          = (K <  S0)
    payoffSigns = sign(K - S0)

    u = 1e-2

    tau = tau / 365.0
  
    for (i in seq(m))
    {
        discFac         = exp(-r[i] * tau[i])
        sgn        = matrix(payoffSigns[i,], ncol(S), n, byrow = TRUE)
        
        payoff     = pmax(sgn * outer(S[i,], K[i,], "-"), 0)
    
        prices[i,] = discFac * apply(payoff, 2, mean)
    
        payoff2    = pmax(sgn * outer((1 + u) * S[i,], K[i,], "-"), 0)
        payoff3    = pmax(sgn * outer((1 - u) * S[i,], K[i,], "-"), 0)

        prices2    = discFac * apply(payoff2, 2, mean)
        prices3    = discFac * apply(payoff3, 2, mean)
    
        deltas[i,] = 100 * (prices2 - prices3) / (2 * u * S0) 
        gammas[i,] = 100 * (prices2 - 2 * prices[i, ] + prices3) / (u * S0)^2

        deltas[i,] = round(deltas[i,], 2)
        gammas[i,] = round(gammas[i,], 2)


        vols[i,]   = impvol(365 * tau[i], S0, K[i, ],
                price = prices[i, ],
                r     = r[i],
                q     = q[i],
                opt   = optTypes[i, ],
                tol   = 1e-7,
                maxit = 1e2)

        calls[i, cix[i,]] = prices[i, cix[i,]]
        puts [i, pix[i,]] = prices[i, pix[i,]]

        calls[i, pix[i,]] =
            prices[i, pix[i,]] - discFac * K[i, pix[i,]] + exp(-q[i] * tau[i]) * S0

        puts [i, cix[i,]] =
            prices[i, cix[i,]] - exp(-q[i] * tau[i]) * S0 + discFac * K[i, cix[i,]]
    }

  
    if (!is.null(names(tau)))
        exp.names = names(tau)
    else if (!is.null(rownames(K)))
        exp.names = rownames(K)
    else
        exp.names = paste("expiration.date", seq(m))

    strike.names = paste(K, optTypes, sep = "")
    dim(strike.names) = dim(K)
    
    rownames(prices)       = rownames(calls)     = rownames(puts) = exp.names
    rownames(vols)         = rownames(optTypes)  = rownames(pix)  = exp.names
    rownames(deltas)       = rownames(gammas)    = rownames(S)    = exp.names
    names(mu)              = rownames(K)         = names(q)       = exp.names
    rownames(strike.names) = rownames(cix)       = names(r)       = exp.names
  
    ans = list()
    for (i in exp.names)
    {
        ans[[i]] = matrix(NA, 8, n)
    
        rownames(ans[[i]]) =
            c("calls", "puts", "calls.ratio",
              "puts.ratio", "vols", "call.delta", "put.delta", "gamma")

        ans[[i]][1, ] = round(calls[i, ], 4)
        ans[[i]][2, ] = round(puts [i, ], 4)

        atmix = which.min(abs(K[i, ] - S0))
    
        patm = puts [i, atmix]
        catm = calls[i, atmix]

        ans[[i]][3, ] = round(catm / calls[i, ], 2)
        ans[[i]][4, ] = round(patm / puts [i, ], 2)
        ans[[i]][5, ] = vols  [i, ]
    
        ans[[i]][6, ] = deltas[i, ] +
            100 * ifelse(pix[i, ], rep(exp(-q[i] * tau[i]), n), rep(0, n))
    
        ans[[i]][7, ] = deltas[i, ] -
            100 * ifelse(cix[i, ], rep(exp(-q[i] * tau[i]), n), rep(0, n))
    
        ans[[i]][8, ] = gammas[i, ]
    
        ans[[i]] = as.data.frame(t(ans[[i]]))
        rownames(ans[[i]]) = K[i,]
    }
    
    ans$vert.sprds = list()
  
    for (i in exp.names)
    {
        ans$vert.sprds[[i]] = round(outer(prices[i,], prices[i,], "-"), 3)

        dimnames(ans$vert.sprds[[i]]) =
            list(paste("+", strike.names[i,], sep = ""),
                 paste("-", strike.names[i,], sep = ""))
    }
  
    ans$prices       = drop(prices)
    ans$vols         = drop(vols)
    ans$deltas       = drop(deltas)
    ans$mu           = drop(mu)
    ans$disc.rates   = drop(r)
    ans$div.yields   = drop(q)
    ans$tau          = drop(tau)
    ans$strikes      = drop(K)
    ans$distribution = drop(S)
    ans$t0           = t0
    ans$S0           = S0
  
    ans
}



######################################################################
######################################################################
######################################################################
##
##  Following section contains functions for fitting the "decaying
##  influence" augmented
##  NGARCH(1, 1) model
##
##  x[t]  = (mu + lambda * sqrt(h[t]) - 0.5 * h[t]) * dt[t]
##              + sqrt(h[t] * dt[t]) * z[t]
##
##  h[t]  = a0[t] + a1 * h[t-1] * (z[t-1] - gamma)^2 + b1 * h[t-1]
##
##  a0[t] = a00 * (1 - beta^age[t]) + a01 * beta^age[t] * y[t]
##
######################################################################
######################################################################
######################################################################



decay.ngarch11 = function(x, dt, y, age, par, ll.only = FALSE)
    ##
    ##  Description: Returns list including loglikelihood,
    ##    conditional variances and residuals for the given
    ##    time series x and parameters.
    ##
    ##    mu      = par[1],
    ##    lambda  = par[2],
    ##    a00     = par[3],
    ##    a01     = par[4],
    ##    beta    = par[5],
    ##    a1      = par[6],
    ##    b1      = par[7],
    ##    gamma   = par[8],
    ##    h0      = par[9],
    ##
{
    .Call("ngarch11_decay", x, dt, y, age, par, ll.only, PACKAGE = "uri")
}



fit.decay.ngarch11 = function(x, dt, y, age,
    method = c("ga", "N", "ga"),
    init   = NULL,
    maxit  = 100,
    tol    = sqrt(.Machine$double.eps),
    zero.mean = FALSE,
    enforce.bounds = FALSE,
    fit.init.state = FALSE, ...)
{
    npar = 8 + as.logical(fit.init.state) - as.logical(zero.mean)
    ## Define the loglikelihood function
    llfn = function(par)
    {
        if (zero.mean) par = c(0, par)
    
        .Call("ngarch11_decay",
              x, dt, y, age, par,
              ll.only = TRUE,
              PACKAGE = "uri")$ll
    }
  
    ansNames =
        c("mu", "lambda", "a00", "a01", "beta",
          "a1", "b1", "gamma", "h0")[1:npar]

    if (zero.mean)
        ansNames = ansNames[-1]
  
    ## Prepare initial parameters
    if (is.null(init))
    {
        require(tseries)
        init = runif(npar, 0, 1e-8)

        names(init) = ansNames

        init[c("a00", "a1", "b1")] =
            garch(x / sqrt(dt), trace = FALSE)$coef[c("a0", "a1", "b1")]
    }
  
    init = as.matrix(init)

    ## Fit the model
    fit = optim.uri(llfn,
        method   = method,
        init     = init,
        npar     = npar,
        maxit    = maxit,
        lower    = NULL,
        upper    = NULL,
        minimize = FALSE,
        tol      = tol,
        enforce.bounds = enforce.bounds,
        ...)

    ## Work with the output a little
    if (!is.null(fit$par))
    {
        rownames(fit$par) = ansNames

        j          = which.max(fit$val)
        fit$x      = x
        fit$ll     = fit$value
        fit$value  = NULL

        tmp     = decay.ngarch11(x, dt, y, age, par = fit$par[, j])
        fit$h   = tmp$h
        fit$res = tmp$res
        fit$aic = aic(fit$ll, npar)
        fit$bic = bic(fit$ll, npar, n)

        names(fit$h) = names(fit$x)
    }

    fit
}



sim.decay.ngarch11 = 
    function(n, paths = 1, dt, y, age,
             par      = NULL,
             V        = NULL,
             z        = numeric(0))
{
    .Call("sim_decay_ng11", n, paths, dt, y, age, par, as.double(z),
          PACKAGE = "uri")
}



######################################################################
######################################################################
######################################################################
##
##  Following section contains functions for fitting the augmented
##  NGARCH(1, 1) model
##
##  x[t]  = mu + lambda * sqrt(h[t]) - 0.5 * h[t] + sqrt(h[t]) * z[t]
##  h[t]  = a0[t] + a1 * h[t-1] * (z[t-1] - gamma)^2 + b1 * h[t-1]
##  a0[t] = a00 + a01 * y[t]
##
######################################################################
######################################################################
######################################################################



aug.ngarch11 = function(x, y, par, ll.only = FALSE, zeroIntercept = TRUE)
    ##
    ##  Description: Returns list including loglikelihood,
    ##    conditional variances and residuals for the given
    ##    time series x and parameters.
    ##
    ##    mu      = par[1],
    ##    lambda  = par[2],
    ##    a00     = par[3],
    ##    a01     = par[4],
    ##    a1      = par[5],
    ##    b1      = par[6],
    ##    gamma   = par[7],
    ##    h0      = par[8],
    ##
{
    if (zeroIntercept) par = c(0, par)
  
    .Call("ngarch11_aug", x, y, par, ll.only, PACKAGE = "uri")
}



fit.aug.ngarch11 = function(x, y,
    method = c("ga", "gr"),
    init   = NULL,
    maxit  = 20,
    lower  = c(-Inf, -Inf,   0,   0,   0,   0, -Inf,   0),
    upper  = c( Inf,  Inf, Inf, Inf, Inf, Inf,  Inf, Inf),
    tol    = sqrt(.Machine$double.eps),
    zeroIntercept  = TRUE,
    enforce.bounds = FALSE,
    fit.init       = FALSE, ...)
{
    n = length(x)

    x = as.double(x)
    y = as.double(y)
  
    npar = 7 + as.logical(fit.init) - as.logical(zeroIntercept)

    ## Define the loglikelihood function
    llfn = function(par)
    {
        if (zeroIntercept) par = c(0, par)
    
        .Call("ngarch11_aug", x, y, par, TRUE, PACKAGE = "uri")$ll
    
    }

    ansNames = c("mu", "lambda", "a00", "a01", "a1", "b1",  "gamma", "h0")

    if (zeroIntercept) ansNames = ansNames[-1]

    ansNames = ansNames[1:npar]
  
    ## Prepare initial parameters
    if (is.null(init))
    {
        require(tseries)
        init = runif(npar, 0, 1e-8)

        names(init) = ansNames

        init[c("a00", "a1", "b1")] =
            garch(x, trace = FALSE)$coef[c("a0", "a1", "b1")]
    }

    init = as.matrix(init)

    ## Fit the model
    fit = optim.uri(llfn,
        method   = method,
        init     = init,
        npar     = npar,
        maxit    = maxit,
        lower    = lower[1:npar],
        upper    = upper[1:npar],
        minimize = FALSE,
        tol      = tol,
        enforce.bounds = enforce.bounds,
        ...)

    ## Work with the output a little
    if (!is.null(fit$par))
    {
        rownames(fit$par) = ansNames

        j          = which.max(fit$value)
        fit$x      = x
        fit$ll     = fit$value[j]
        fit$par    = fit$par[, j]
        fit$value  = NULL

        tmp     = aug.ngarch11(x, y, par = fit$par, zeroIntercept = zeroIntercept)
        fit$h   = tmp$h
        fit$res = tmp$res
        fit$aic = aic(fit$ll, npar)
        fit$bic = bic(fit$ll, npar, n)

        names(fit$h) = names(fit$x)
    }

    fit
}



sim.aug.ngarch11 = function(pathLen, numPaths, y, par, x0, h0, z = numeric(0))
{
    .Call("sim_aug_ng11",
          pathLen, numPaths, y, par, x0, h0, as.double(z),
          PACKAGE = "uri")
}



######################################################################
######################################################################
######################################################################
##
##  Following section contains functions for fitting the seasonal
##  NGARCH(1, 1) model
##
##  x[t]  = mu + lambda * sqrt(h[t]) - 0.5 * h[t] + sqrt(h[t]) * z[t]
##  h[t]  = a0[t] + a1 * h[t-1] * (z[t-1] - gamma)^2 + b1 * h[t-1]
##  a0[t] = a00 * exp(a01 * cos(2*pi/365*yday[t]))
##
######################################################################
######################################################################
######################################################################



seasonal.ngarch11 = function(x, yday, par, ll.only = FALSE)
  ##
  ##  Description: Returns list including loglikelihood,
  ##    conditional variances and residuals for the given
  ##    time series x and parameters.
  ##
  ##    mu      = par[1],
  ##    lambda  = par[2],
  ##    a00     = par[3],
  ##    a01     = par[4],
  ##    omega   = par[5],
  ##    a1      = par[6],
  ##    b1      = par[7],
  ##    gamma   = par[8],
  ##    h0      = par[9],
  ##
{
    .Call("ngarch11_seasonal", x, yday, par, ll.only, PACKAGE = "uri")
}



fit.seasonal.ngarch11 =
  function(x,
           yday,
           method = c("ga", "N", "ga"),
           init   = NULL,
           maxit  = 20,
           lower  = c(-Inf, -Inf,   0,   0,   0,   0,   0, -Inf,   0),
           upper  = c( Inf,  Inf, Inf, Inf, Inf, Inf, Inf,  Inf, Inf),
           tol    = sqrt(.Machine$double.eps),
           enforce.bounds = FALSE,
           fit.init.state = FALSE, ...)
{
    n = length(x)

    x    = as.double(x)
    yday = as.double(yday)
  
    npar = ifelse(fit.init.state, 9, 8)

    ## Define the loglikelihood function
    llfn = function(par)
      {
    
          .Call("ngarch11_seasonal", x, yday, par, TRUE, PACKAGE = "uri")$ll
    
      }

  
    ansNames =
      c("mu", "lambda", "a00", "a01", "omega",
        "a1", "b1",  "gamma", "h0")[1:npar]

  
    ## Prepare initial parameters
    if (is.null(init))
      {

          require(tseries)
          init = runif(npar, 0, 1e-8)

          names(init) = ansNames
      
          init[c("a00", "a1", "b1")] =
            garch(x, trace = FALSE)$coef[c("a0", "a1", "b1")]
          
      }

  
    init = as.matrix(init)

    ## Fit the model
    fit = optim.uri(llfn,
      method   = method,
      init     = init,
      npar     = npar,
      maxit    = maxit,
      lower    = lower[1:npar],
      upper    = upper[1:npar],
      minimize = FALSE,
      tol      = tol,
      enforce.bounds = enforce.bounds,
      ...)

    ## Work with the output a little
    if (!is.null(fit$par))
      {
    
          rownames(fit$par) = ansNames

          j          = which.max(fit$val)
          fit$x      = x
          fit$ll     = fit$value
          fit$value  = NULL
          
          tmp     = seasonal.ngarch11(x, yday, par = fit$par[, j])
          fit$h   = tmp$h
          fit$res = tmp$res
          fit$aic = aic(fit$ll, npar)
          fit$bic = bic(fit$ll, npar, n)

          names(fit$h) = names(fit$x)
          
      }

    fit
}

