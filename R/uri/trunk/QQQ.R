###############################################################################
# Load Systematic Investor Toolbox (SIT)
# http://systematicinvestor.wordpress.com/systematic-investor-toolbox/
###############################################################################
setInternet2(TRUE)
con = gzcon(url('http://www.systematicportfolio.com/sit.gz', 'rb'))
source(con)
close(con)

library( uri )

##*****************************************************************
##Load historical data
##****************************************************************** 
load.packages( 'quantmod' )	


get.price.data = function( data, tickers, from.date )
{
    getSymbols( tickers,
               src         = 'yahoo',
               from        = from.date,
               env         = data,
               auto.assign = TRUE )

    for( i in tickers )
        data[[ i ]] = adjustOHLC( data[[ i ]], use.Adjusted = TRUE )		

    bt.prep( data, align = 'keep.all', dates = '1991::')
}



run.ngarch = function(
    data,
    tickers,
    window.size = 252*2,
    last.prices = NULL,
    last.date   = Sys.time() )
{
    prices  = matrix( NA, nrow( data[[ tickers[1] ]] ), length( tickers ))

    for( j in 1:length( tickers ))
        prices[, j] = data[[ tickers[j] ]][,6]

    colnames( prices ) = tickers

    dates = as.POSIXct( paste( data[[ "dates" ]], "20:00 GMT" ))
    
    if( !is.null( last.prices ))
    {
        prices = rbind( prices, last.prices )
        dates  = c( dates, as.POSIXct( last.date ))
    }
    
    
    ## Get returns
    returns  = matrix( NA, nrow( prices ) - 1, ncol( prices ))
    daysdiff = as.numeric(
        difftime( dates[ -1 ], dates[ -length( dates )], unit = "days" ))
    
    for( j in 1 : ncol( prices ))
        returns[ , j ] = diff( log( prices[ , j ] )) / sqrt( daysdiff )

    window.size  = min( window.size, nrow( returns ))
    indx.include = nrow( returns ) - seq( window.size - 1, 0 )

    returns = returns[ indx.include, ]

    if( is.null( data[[ "fits" ]] ))
        fits = list()
    else
        fits = data[[ "fits" ]]

    for( j in 1:ncol( returns ))
    {
        x = returns[ , j ]
        x = x[ !is.na( x ) ]
        
        fits[[ tickers[ j ]]] = fitNgarch11.b(
                x,
                init = c( var( x ), runif( 3, 0, 1e-4 )),
                stopLag = 50, minit = 100, maxit = 5e4, fitInit = TRUE,
                opt = list( grad = TRUE ))

        fits[[ tickers[ j ]]]$par = fits[[ j ]]$par[ 1:4 ]
    
        fits[[ tickers[ j ]]] = fitNgarch11.b(
                x,
                init = fits[[ j ]]$par,
                pop = 10,
                stopLag = 1, minit = 100,
                maxit = 5e3, fitInit = TRUE,
                opt = list( grad = TRUE ))

        fits[[ j ]]$par = fits[[ j ]]$par[ 1:4 ]
    
        fits[[ j ]]$par = as.matrix( c( 0, fits[[ j ]]$par ))

        row.names( fits[[ j ]]$par ) = c( "lambda", "a0", "a1", "b1", "gamma" )
    }

    data[[ "fits" ]] = fits

    return( fits )
}



get.ng11.dist = function( fits, exp.dts, freq = 1 )
{
    dists = list()

    for( j in 1:length( fits ))
    {
        fits[[ j ]]$freq = freq
    
        dists[[ j ]] = ngarch11prices( exp.dts, fits[[ j ]], paths = 100e3 )
        dists[[ j ]]$exp.dts = exp.dts
    }
    names( dists ) = names( fits )

    return( dists )
}


get.ng11.vols = function( dists, und.prices, strikes )
{
    vols = list()

    for( j in 1:length( dists ))
    {
        vols[[ j ]] = ngarch11vols(
                dists[[ j ]],
                S0 = und.prices[ j ],
                strikes = strikes[,j] )    
    }
    names( vols ) = tickers

    return( vols )
}


data      = new.env()
from.date = '2003-01-01'
tickers   = c( "SPY", "QQQ", "GLD" )

get.price.data( data, tickers, from.date )

window.size = 252*3

last.prices = c( 168.57, 77.73, 131.44 )
##last.prices = NULL
fits  = run.ngarch( data, tickers,
    window.size,
    last.price = last.prices )

exp.dts = c(
    "2013-09-13 20:00:00 GMT",
    "2013-09-20 20:00:00 GMT",
    "2013-09-27 20:00:00 GMT",
    "2013-10-18 20:00:00 GMT",
    "2013-11-15 20:00:00 GMT",
    "2013-12-31 20:00:00 GMT" )


dists = get.ng11.dist( fits, exp.dts )

und.prices = c( 168.73, 77.85, 131.50 )
strikes    = t(round( und.prices * 2, 0 ) / 2 + t(cbind( -10:10, -10:10, -10:10 ) / 2 ))
vols       = get.ng11.vols( dists, und.prices, strikes )

##vols = get.ng11.vols( dists[1], und.prices[1], strikes )

tmp = vols[[3]]$vol[4,]

