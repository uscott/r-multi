
interp.crv = function( x.known, y.known, x.arg )
{
    if( length( x.known ) != length( y.known ))
        stop( 'Length mismatch' )
    
    return( .Call(
                'R_interp_crv',
                as.double( x.known ),
                as.double( y.known ),
                as.double( x.arg ),
                PACKAGE = 'uri' ))
}


interp.sfc = function( t.in, m.in, v.in, t.arg, m.arg )
{
    if( length( t.in ) != length( m.in ) | length( m.in ) != length( v.in ))
        stop( 'Length mismatch' )

    if(      length( t.arg ) == 1 & length( m.arg ) > 1 )
        t.arg = rep( t.arg, length( m.arg ))
    else if( length( m.arg ) == 1 & length( t.arg ) > 1 )
        m.arg = rep( m.arg, length( t.arg ))
    
    return( .Call(
                'R_interp_sfc',
                as.double( t.in ),
                as.double( m.in ),
                as.double( v.in ),
                as.double( t.arg ),
                as.double( m.arg ),
                PACKAGE = 'uri' ))
}

