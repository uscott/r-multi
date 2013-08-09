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


from.date = '2006-01-01'

tickers   = c( "SPY", "QQQ", "GLD" )
data      = new.env()

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

window.size = 252*2
indx.include = nrow( returns ) - seq( window.size - 1, 0 )

returns = returns[ indx.include, ]

fits = list()

for( j in 1:ncol( returns ))
{
    fits[[ j ]] = fitNgarch11.b( returns[ , j ], 
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

exDts = c(
    "2013-08-09 20:00:00 GMT",
    "2013-08-16 20:00:00 GMT",
    "2013-08-23 20:00:00 GMT",
    "2013-12-31 20:00:00 GMT" )
dists = list()

for( j in 1:length(fits))
{
    fits[[ tickers[ j ]]]$freq = 1
    
    dists[[ j ]] = ngarch11prices( exDts, fits[[ j ]], paths = 100e3 )
}


current.price = c( 169.97, 76.77, 126.84 )
strikes = t(round( current.price * 2, 0 ) / 2 + t(cbind( -6:6, -6:6, -6:6 ) / 2 ))
vols = list()

for( j in 1:length(fits))
{
    vols[[ j ]] = ngarch11vols(
            dists[[ j ]],
            S0 = current.price[j],
            strikes = strikes[,j] )    
}



