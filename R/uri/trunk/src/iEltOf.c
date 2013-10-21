
#include "uri.h"


int iEltOf( long val, SEXP t, long *index )
{
    int retval = 0;
    long tlen, i;

    if ( !index ) 
        error ("null pointer");

    PROTECT( t = AS_INTEGER( t ));
    tlen = length( t );

    for( i = 0; i < tlen && !retval; ++i )
        if ( retval = ( val == INTEGER( t )[ i ] )) 
            *index = i;

    UNPROTECT_PTR( t );

    return retval;
}
