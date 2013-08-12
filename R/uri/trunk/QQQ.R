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

run.ngarch = function( data, tickers, from.date, window.size = 252*2 )
{
    getSymbols( tickers,
               src         = 'yahoo',
               from        = from.date,
               env         = data,
               auto.assign = TRUE )

    for( i in ls( data ))
        data[[ i ]] = adjustOHLC( data[[ i ]], use.Adjusted = TRUE )		

    bt.prep( data, align = 'keep.all', dates = '1991::')

    prices  = matrix( NA, nrow( data[[ tickers[1] ]] ), length( tickers ))

    for( j in 1:length( tickers ))
        prices[, j] = data[[ tickers[j] ]][,6]

    colnames( prices ) = tickers

    dates = as.POSIXct( data[[ "dates" ]] )

    ## Remove NAs
    tmp       = apply( prices, 1, sum )
    indx.omit = which( is.na( tmp ))

    if( length( indx.omit ))
    {
        prices = prices[ -indx.omit, ]
        dates  = dates[ -indx.omit ]
    }

    ## Get returns
    returns  = matrix( NA, nrow( prices ) - 1, ncol( prices ))
    daysdiff = as.numeric( diff( dates ))
    for( j in 1 : ncol( prices ))
    {
        tmp            = prices[ , j ]
        returns[ , j ] = log( tmp[ -1 ] / tmp[ -length( tmp ) ] ) / sqrt( daysdiff )
    }

    indx.include = nrow( returns ) - seq( window.size - 1, 0 )

    returns = returns[ indx.include, ]

    fits = list()

    for( j in 1:ncol( returns ))
    {
        fits[[ j ]] = fitNgarch11.b( returns[ , j ], #init = fits2[[j]]$par[-1],
                stopLag = 50, minit = 100, maxit = 5e4, fitInit = TRUE,
                opt = list( grad = TRUE ))

        fits[[ j ]]$par = fits[[ j ]]$par[ 1:4 ]
    
        fits[[ j ]] = fitNgarch11.b( returns[ , j ], init = fits[[ j ]]$par,
                pop = 10,
                stopLag = 1, minit = 100, maxit = 5e3, fitInit = TRUE,
                opt = list( grad = TRUE ))

        fits[[ j ]]$par = fits[[ j ]]$par[ 1:4 ]
    
        fits[[ j ]]$par = as.matrix( c( 0, fits[[ j ]]$par ))

        row.names( fits[[ j ]]$par ) = c( "lambda", "a0", "a1", "b1", "gamma" )
    }

    names( fits ) = tickers

    return( fits )
}

get.ng11.dist = function( fits, exp.dts, freq = 1 )
{
    dists = list()

    for( j in 1:length( fits ))
    {
        fits[[ j ]]$freq = freq
    
        dists[[ j ]] = ngarch11prices( exDts, fits[[ j ]], paths = 100e3 )
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
from.date = '2006-01-01'
tickers   = c( "SPY", "QQQ", "GLD" )

window.size = 252*2

exp.dts = c(
    "2013-08-16 20:00:00 GMT",
    "2013-08-23 20:00:00 GMT",
    "2013-09-20 20:00:00 GMT",
    "2013-11-15 20:00:00 GMT",
    "2013-12-31 20:00:00 GMT" )

und.prices = c( 169.10, 76.73, 129.23 )
##current.price = c( last( prices[,1] ), last( prices[,2] ), last( prices[,3] ))
strikes = t(round( und.prices * 2, 0 ) / 2 + t(cbind( -10:10, -10:10, -10:10 ) / 2 ))

fits  = run.ngarch( data, tickers, from.date, window.size )
dists = get.ng11.dist( fits, exp.dts )
vols  = get.ng11.vols( dists, und.prices, strikes )
