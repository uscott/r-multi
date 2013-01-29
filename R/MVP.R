## Function for minimum variance portfolio calculation

mvp = function( returns )
    ## returns is a matrix of asset returns with each column corresponding to an asset
{
    N = ncol( returns )

    if ( N <= 0 )
        return( numeric(0))

    A = matrix( NA, N, N ) # covariance matrix

    for( i in 1:N )
        for( j in 1:i )
        {
            A[ i, j ] = cov( returns[ , i ], returns[ , j ] )

            if ( i != j )
                A[ j, i ] = A[ i, j ]
        }

    Iden = matrix( 0, N, N ) # identity matrix

    for( i in 1:N )
        Iden[ i, i ] = 1

    B    = solve( A, Iden ) # inverse of covariance matrix A
    
    ones = matrix( 1, N, 1 ) # column vector of 1s

    w    = ( B %*% ones ) / as.double( t(ones) %*% B %*% ones )

    return( w )
}
