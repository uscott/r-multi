## Function for minimum variance portfolio calculation

mvp = function( returns, use = "everything" )
    ## returns is a matrix of asset returns with each column corresponding to an asset
{
    N = ncol( returns )

    if ( N <= 0 )
        return( numeric( 0 ))

    A = cov( returns, use = use ) # covariance matrix

    Id = matrix( 0, N, N ) # identity matrix

    for( i in 1:N )
        Id[ i, i ] = 1

    B    = solve( A, Id ) # inverse of covariance matrix A
    ones = matrix( 1, N, 1 ) # column vector of 1s
    w    = ( B %*% ones ) / as.double( t(ones) %*% B %*% ones ) # weights

    rownames( w ) = colnames( returns )
    
    return( w )
}
