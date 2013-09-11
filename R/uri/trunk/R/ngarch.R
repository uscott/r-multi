


evalNgarch11 = function(x, pars, llonly = FALSE)
    ##
    ##  Description: Returns list including loglikelihood,
    ##    conditional variances and residuals for the given
    ##    time series x and parameters.
    ##
{
    .Call( "ngarch11", x, pars, llonly, PACKAGE = "uri")
}



fitNgarch11 = function( x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, fitInit = FALSE, option = NULL)
{
    .Call( "fit_ngarch11", x, init, fitInit,
          popSize, tol,
          stopLags, minit, maxit, option,
          PACKAGE = "uri")
}



bootfitNgarch11 = function(x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, fitInit = FALSE, option = NULL, numBoots = 20 )
{
    .Call( "BootfitNgarch11", x, init, fitInit, popSize, tol,
          stopLags, minit, maxit, option, numBoots,
          PACKAGE = "uri" )
}



simNgarch11 = function(n, paths = 1, pars, x0 = NULL, h0 = NULL, z = NULL)
{
    .Call( "sim_ngarch11", n, paths, pars, x0, h0, z, PACKAGE = "uri")
}



bootsimNgarch11 = function( n, paths = 1, pars, x0 = NULL, h0 = NULL, z = NULL )
{
    .Call( "sim_ngarch11_2", n, paths, pars, x0, h0, z, PACKAGE = "uri" )
}



get.dist.ng11 = function(
    exp.dts,
    modelInfo,
    paths = 1e4,
    r = 0, ## discount rates
    q = 0, ## dividend yields
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
  
    if ( any( is.null( modelInfo$par ),
             is.null( modelInfo$res ),
             boot & is.null( modelInfo$bootPar ),
             is.null( modelInfo$freq )))
        stop( "incomplete NGARCH(1, 1) model info given" )

    t0   = as.POSIXct( Sys.time())
    tau  = as.double( difftime( as.POSIXct( exp.dts ), t0, "days" ))
    
    names(tau) = names( exp.dts )
    
    tau  = tau[ tau > 0 ]
    freq = modelInfo$freq # Frequency in days
    tau  = sort( tau ) / freq

    tauInt  = trunc( tau )
    tauFrac = tau - tauInt

    m = length( tau )
    n = ceiling( max( tau ))

    ## Simulate terminal prices
    
    if ( boot )
    {
        N        = ceiling( paths / ncol( modelInfo$bootPar ))

        bootPars = matrix(
            rep( modelInfo$bootPar, N * ncol( modelInfo$bootPar )),
            nrow( modelInfo$bootPar ),
            N * ncol( modelInfo$bootPar ))
        
        bootPars             = bootPars[, 1:paths ]
        rownames( bootPars ) = rownames( modelInfo$bootPar )

        h0 = modelInfo$bootVar[  nrow( modelInfo$bootVar ),  ]
        x0 = modelInfo$bootPath[ nrow( modelInfo$bootPath ), ]
        h0 = rep( h0, N )[ 1:paths ]
        x0 = rep( x0, N )[ 1:paths ]
        
        if ( risk.neutral )
        {
            bootPars[ "gamma", ] = bootPars[ "gamma", ] + bootPars[ "lambda", ]
            bootPars[ "lambda",] = 0
        }

        rsim = bootsimNgarch11( n = n, paths = paths, bootPars, z = modelInfo$res,
            x0 = x0, h0 = h0 )$x        
    }
    else
    {
        lambda = modelInfo$par[ "lambda", 1 ]
        a0     = modelInfo$par[ "a0",     1 ]
        a1     = modelInfo$par[ "a1",     1 ]
        b1     = modelInfo$par[ "b1",     1 ]
        gamma  = modelInfo$par[ "gamma",  1 ]
        h0     = last( modelInfo$h )
        x0     = last( modelInfo$x )
        
        if ( risk.neutral )
        {
            gamma  = lambda + gamma
            lambda = 0
        }
        
        pars = list( lambda = lambda, a0 = a0, a1 = a1, b1 = b1, gamma = gamma)
        
        rsim = simNgarch11( n = n, paths = paths,
            par = pars, z = modelInfo$res, x0 = x0, h0 = h0 )$x
    }
    
    rsim = rsim[, !is.na( rsim[ n, ] ), drop = FALSE ]        
    r1   = apply( rsim, 2, cumsum )[ tauInt, , drop = FALSE ]
    r2   = rsim[ tauInt + 1, , drop = FALSE ]
    mr2  = apply( r2, 1, mean)
    r2   = ( r2 - mr2 ) * sqrt( tauFrac ) + mr2 * tauFrac
    S    = exp( r1 + r2 )

    if ( risk.neutral )
        S = S / apply( S, 1, mean ) * exp( tau * (r - q) / 365 * freq )
    else
        S = S * exp( -tau * q / 365 * freq )
 
    stopifnot( nrow( S ) == m )
   
    rownames( S ) = names( tau )

    retval = list(
        S       = t(S),
        tau     = tau * freq,
        freq    = freq,
        r       = r,
        q       = q,
        t0      = t0,
        exp.dts = exp.dts )

    return( retval )
}



get.vols.ng11 = function(
    distInfo,
    underlying =  1,
    strikes )
    ## Function for pricing options.
    ## distInfo is a list of length 2 containing
    ## the physical terminal distribution and the corresponding
    ## times to expiration
    ## NGARCH(1, 1) model.
    ##
{
    S0   = underlying
    S    = t( distInfo$S ) * S0
    tau  = distInfo$tau
    freq = distInfo$freq
    t0   = distInfo$t0
    r    = distInfo$r
    q    = distInfo$q
  
    stopifnot(
        is.matrix( S )   &
        is.vector( tau ) &
        is.double( freq )&
        is.double( S0 ))

    stopifnot( all( tau > 0, freq > 0, S0 > 0 ))
    stopifnot( all( length(strikes), strikes > 0 ))

    tau     = sort( tau )
    m       = length( tau )
    n       = ceiling( max( tau ))
    tauInt  = trunc( tau )
    tauFrac = tau - tauInt

    if ( 1 == length( r ))
        r = rep( r, m )
    if ( 1 == length( q ))
        q = rep( q, m )
  
    rmq = r - q

    if ( is.null( dim( S )))
        S = t(S)
  
    stopifnot( nrow( S ) == m )
  
    mu = apply( S, 1, mean ) / S0
    K  = strikes
    rm(strikes)
  
    if( is.vector( K ))
        K = matrix( K, m, length( K ), byrow = TRUE )
    else if ( is.matrix( K ) & identical( 1 == nrow( K ), TRUE ))
        K = matrix( as.vector( K ), m, ncol( K ), byrow = TRUE )
  
    K = t( apply( K, 1, sort ))
    n = ncol( K )

    stopifnot( is.matrix( K ))
  
    calls  = matrix( NA, m, n )
    puts   = matrix( NA, m, n )
    prices = matrix( NA, m, n )
    deltas = matrix( NA, m, n )
    gammas = matrix( NA, m, n )
    vols   = matrix( NA, m, n )
  
    optTypes    = ifelse( K < S0, "p", "c")
    cix         = ( K >= S0 )
    pix         = ( K <  S0 )
    payoffSigns = ifelse( K >= S0, 1, -1 )

    u   = 1e-2
    tau = tau / 365.0
  
    for( i in 1 : m )
    {
        discFac    = exp( -r[ i ] * tau[ i ] )
        sgn        = matrix( payoffSigns[ i, ], ncol(S), n, byrow = TRUE )        
        payoff     = pmax( sgn * outer( S[ i, ], K[ i, ], "-" ), 0 )
        prices[i,] = discFac * apply( payoff, 2, mean )
        payoff2    = pmax( sgn * outer( (1 + u) * S[i,], K[i,], "-"), 0 )
        payoff3    = pmax( sgn * outer( (1 - u) * S[i,], K[i,], "-"), 0 )
        prices2    = discFac * apply( payoff2, 2, mean )
        prices3    = discFac * apply( payoff3, 2, mean )
    
        deltas[i,] = 100 * ( prices2 - prices3 ) / ( 2 * u * S0 ) 
        gammas[i,] = 100 * ( prices2 - 2 * prices[ i, ] + prices3 ) / (u * S0)^2

        deltas[i,] = round( deltas[ i, ], 2 )
        gammas[i,] = round( gammas[ i, ], 2 )

        vols[i,]   = impvol(
                365 * tau[i], S0, K[ i, ],
                price = prices[ i, ],     r = r[i], q     = q[ i ],
                opt   = optTypes[ i, ], tol = 1e-7, maxit = 1e3 )

        calls[ i, ] = prices[ i, ]
        puts [ i, ] = prices[ i, ]

        parity = exp( -q[ i ] * tau[ i ] ) * S0 - discFac * K[i,]
    
        calls[ i, pix[ i, ]] = calls[ i, pix[ i,]] + parity[ pix[ i, ]]
        puts [ i, cix[ i, ]] = puts [ i, cix[ i,]] - parity[ cix[ i, ]]        
    }
    
    if( !is.null( names( tau )))
        exp.names = names( tau )
    else if( !is.null( rownames( K )))
        exp.names = rownames( K )
    else
        exp.names = paste( "expi", 1 : m, sep = "" )

    strike.names = paste( K, optTypes, sep = "" )
    dim(strike.names) = dim(K)
    
    rownames( prices )       = rownames( calls )    = rownames( puts ) = exp.names
    rownames( vols )         = rownames( optTypes ) = rownames( pix )  = exp.names
    rownames( deltas )       = rownames( gammas )   = rownames( S )    = exp.names
    names( mu )              = rownames( K )        = names( q )       = exp.names
    rownames( strike.names ) = rownames( cix )      = names( r )       = exp.names
    
    colnames( prices ) = colnames( vols ) = colnames( deltas ) = K[ 1, ]
    
    ans = list()
    
    for( i in exp.names )
    {
        ans[[ i ]] = matrix( NA, 8, n )
    
        rownames(ans[[ i ]]) =
            c("calls", "puts", "calls.ratio",
              "puts.ratio", "vols", "call.delta", "put.delta", "gamma")

        ans[[ i ]][ 1, ] = round( calls[ i, ], 3 )
        ans[[ i ]][ 2, ] = round( puts [ i, ], 3 )

        atmix = which.min( abs( K[ i, ] - S0 ))
    
        patm = puts [ i, atmix ]
        catm = calls[ i, atmix ]

        ans[[ i ]][ 3, ] = round( catm / calls[i, ], 3 )
        ans[[ i ]][ 4, ] = round( patm / puts [i, ], 3 )
        ans[[ i ]][ 5, ] = round( vols[ i, ], 3 )
    
        ans[[ i ]][ 6, ] = deltas[ i, ] +
            100 * ifelse(
                      pix[ i, ],
                      rep( exp( -q[ i ] * tau[ i ]), n ),
                      rep( 0, n ))
    
        ans[[ i ]][ 7, ] = deltas[ i, ] -
            100 * ifelse(
                      cix[ i, ],
                      rep( exp( -q[ i ] * tau[ i ]), n ),
                      rep( 0, n ))
    
        ans[[ i ]][ 8, ] = gammas[ i, ]
    
        ans[[ i ]] = as.data.frame( t( ans[[ i ]] ))
        
        rownames( ans[[ i ]] ) = K[ i, ]
    }

    ans$vert.sprds  = list()
    ans$cs = ans$ps = list()
  
    for ( i in exp.names )
    {        
        ans$vert.sprds[[ i ]] =
            round( outer( prices[ i, ], prices[ i, ], "-" ), 3 )
        
        ans$cs[[ i ]] = round( outer( calls[ i, ], calls[ i,], "-"), 3 )
        ans$ps[[ i ]] = round( outer( puts[  i, ], puts[  i,], "-"), 3 )
    
        dimnames( ans$vert.sprds[[ i ]] ) =
            list( paste( "+", strike.names[ i, ], sep = "" ),
                  paste( "-", strike.names[ i, ], sep = "" ))

        dimnames( ans$cs[[ i ]] ) = dimnames( ans$ps[[ i ]] ) =
            list(paste( "+", K[ i, ], sep = "" ),
                 paste( "-", K[ i, ], sep = "" ))    
    }
  
    ans$prices       = drop( prices )
    ans$vols         = drop( vols )
    ans$deltas       = drop( deltas )
    ans$mu           = drop( mu )
    ans$disc.rates   = drop( r )
    ans$div.yields   = drop( q )
    ans$tau          = drop( tau ) * 365
    ans$strikes      = drop( K )
    ans$distribution = drop( S )
    ans$t0           = t0
    ans$underlying   = S0
    ans$exp.dts      = distInfo$exp.dts
    
    return( ans )
}


fitNgarch11.b = function( x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, fitInit = FALSE, option = NULL)
{
    .Call(
        "fit_ngarch11_B",
        x,
        init,
        fitInit,
        popSize,
        tol,
        stopLags,
        minit,
        maxit,
        option,
        PACKAGE = "uri" )
}


simNgarch11.b = function(n, paths = 1, pars, x0 = NULL, h0 = NULL, z = NULL)
{
    .Call( "sim_ngarch11_B", n, paths, pars, x0, h0, z, PACKAGE = "uri")
}
