#include "uri.h"
#include<stdio.h>

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This file contains functions for calculating the BS formula, greeks,
**  implied vol and other stuff.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#define IS_VALID_OPTION_TYPE( o ) ( (o) == 'c' || (o) == 'p' || (o) == 's' )
#define INV_SQRT_2PI (0.398942280401)


static double d2_from_d1( double d1, double T, double vol)
{
	return d1 - vol * sqrt( T );
}


static double get_d1( double T, double S, double K, double vol, double r, double q)
{
    return ( log( S / K ) + ( r - q + 0.5 * vol * vol ) * T ) / ( vol * sqrt( T ));
}


static double get_d2( double T, double S, double K, double vol, double r, double q)
{
    return ( log( S / K ) + ( r - q - 0.5 * vol * vol ) * T ) / ( vol * sqrt( T ));
}


static void get_d1_and_d2(
	double T, double S, double K, double vol, double r, double q, double *d1, double *d2
)
{
	if ( d1 != NULL )
		*d1 = get_d1( T, S, K, vol, r, q );

	if ( d2 != NULL )
	{
		if ( d1 != NULL )
			*d2 = d2_from_d1( *d1, T, vol );
		else
			*d2 = get_d2( T, S, K, vol, r, q );
	}
}


static double intrinsic( double T, double S, double K, double r, double q, char o )
{
	double p = exp( -q * T ) * S - exp( -r * T ) * K;

	switch ( o )
	{
	case 'c':
		return MAX( 0, +p);
	case 'p':
		return MAX( 0, -p);
	}
	return ABS( p );
}


static double zero_spot_delta( double T, double q, char opt )
{
	switch ( opt )
	{
	case 'c':
		return 0;
	case 'p':
	case 's':
		return -exp( -q * T );
	}

	return NA_REAL;
}


static double zero_spot_price( double T, double K, double r, char o )
{
	switch ( o )
	{
	case 'c':
		return 0;
	case 'p':
	case 's':
		return exp( -r * T ) * K;
	}

	return NA_REAL;
}


static double zero_spot_theta( double T, double K, double r, char o )
{
	switch ( o )
	{
	case 'c':
		return 0;
	case 'p':
	case 's':
		return r * K * exp( -r * T );
	}

	return NA_REAL;
}


static double zero_strike_delta( double T, double q, char opt )
{
	switch ( opt)
	{
	case 'p':
		return 0;
	case 'c':
	case 's':
		return exp( -q * T );
	}

	return NA_REAL;
}


static double zero_strike_price( double T, double S, double q, char opt)
{
	switch ( opt )
	{
	case 'p':
		return 0;
	case 'c':
	case 's':
		return exp( -q * T ) * S;
	}

	return NA_REAL;
}


double zero_strike_theta( double T, double S, double q, char o )
{
	switch ( o )
	{
	case 'p':
		return 0;
	case 'c':
	case 's':
		return q * S * exp( -q * T );
	}

	return NA_REAL;
}


static double zero_vol_delta( double T, double S, double K, double r, double q, char opt )
{
	double dfq;
	dfq = exp( -q * T );
	S *= dfq;
	K *= exp( -r * T );

	switch ( opt )
	{
	case 'c':
		return S < K ? 0 : dfq;
	case 'p':
		return S < K ? -dfq : 0;
	case 's':
		return S < K ? -dfq : dfq;
	}

	return NA_REAL;
}


static double zero_vol_gamma( double T, double S, double K, double r, double q )
{
	return exp( -q * T ) * S - exp( -r * T ) * K != 0 ? 0 : NA_REAL;
}



double zero_vol_theta( double T, double S, double K, double r, double q, char o )
{
	S *= exp( -q * T );
	K *= exp( -r * T );

	switch ( o )
	{
	case 'c':
		if ( S > K )
			return q * S - r * K;
		else if ( S < K )
			return 0;
		return q * S;
	case 'p':
		if ( S > K )
			return 0;
		else if ( S < K )
			return r * K - q * S;
		return r * K;
	case 's':
		if ( S > K )
			return q * S - r * K;
		else if ( S < K )
			return r * K - q * S;
		return q * S + r * K;
	}

	return NA_REAL;
}



void bs_value(
	const double T,
	const double S,
	const double K,
	const double vol,
	const double r,
	const double q,
    char opt,
	char ret,
	double *ans
)
/*******************************************************************
 *
 *  Description: Computes the BS price, delta or vega and sets *ans to it.
 *
 *  Assumes T is given in *years* (or whatever time units vol, r and q are given in).
 *
 *******************************************************************/
{
    int ok;
    double d1, d2;

    opt = tolower( opt );
    ret = tolower( ret );

    if ( ans == NULL )
        error ( "null pointer passed" );

    *ans = NA_REAL;

	if ( bs_price_chkargs( T, S, K, vol, r, q ))
		return;

	get_d1_and_d2( T, S, K, vol, r, q, &d1, &d2 );

	if ( ret == 'p' )
	{
		if ( T > 0 )
		{

		}
	}
    switch ( ret )
    {
        case 'p':

            if( 0.0 < T )
                *ans = exp(-q * T)* S * NORM_DIST( d1 ) - exp( -r * T)* K * NORM_DIST( d2 );
            else
                *ans = DMAX( S - K, 0 );

            if( opt == 'p' )
                *ans -= exp( -q * T ) * S - exp( -r * T ) * K;
            else if( opt == 's' )
            {
                *ans *= 2.0;
                *ans -= exp( -q * T ) * S - exp( -r * T ) * K;
            }
            break;

        case 'd':
            if( 0.0 < T )
                *ans = exp( -q * T ) * NORM_DIST( d1 );
            else
                *ans = ( S > K ? 1.0 : 0.0 );

            if( opt == 'p' )
                *ans -= exp( -q * T );
            else if( opt == 's' )
            {
                *ans *= 2.0;
                *ans -= exp( -q * T );
            }
            break;

        case 'v':

            if( 0.0 < T )
                *ans = exp( -q * T )* S * sqrt( T ) * NORM_DENS( d1 );
            else
                *ans = 0.0;

            if( opt == 's' )
                *ans *= 2.0;

            break;

        case 'g':

            if( 0.0 < T )
                *ans = exp( -q * T ) * NORM_DENS( d1 ) / (S*vol*sqrt(T));
            else
                *ans = S == K ? R_PosInf : 0.0;

            if( opt == 's' )
                *ans *= 2.0;
			break;

        default:
            *ans = NA_REAL;
            break;
    }
}


double bs_price(
	double T,
	double S,
	double K,
	double vol,
	double q,
	double r,
	char o
)
{
	if ( T == 0 || vol == 0 )
		return intrinsic( T, S, K, r, q, o );

	if ( K == 0 )
		return zero_strike_price( T, S, q, o );

	if ( S == 0 )
		return zero_spot_price( T, K, r, o );

	double d1, d2;
	get_d1_and_d2( T, S, K, vol, r, q, &d1, &d2 );

	S *= exp( -q * T );
	K *= exp( -r * T );

	double val =  S * NORM_DIST( d1 ) - K * NORM_DIST( d2 );

	return val;
}


int bs_price_chkargs(
	const double T,
	const double S,
	const double K,
	const double vol,
	const double r,
	const double q
)
{
	return ( S < 0 || K < 0 || !R_FINITE( q ) || vol < 0 || !R_FINITE( T ) || !R_FINITE( K )
			|| !R_FINITE( vol ) || !R_FINITE( r ) || T < 0.0 );
}


double bs_delta( double T, double S, double K, double vol, double r, double q, char opt )
{
	if ( vol < 0.0 )
	{
		return 2 * zero_vol_delta( T, S, K, r, q, opt) - bs_delta( T, S, K, -vol, r, q, opt);
	}

	if ( S <= 0.0 )
		return zero_spot_delta( T, q, opt );
	if ( K <= 0.0 )
		return zero_strike_delta( T, q, opt );
	if ( vol == 0.0 )
		return zero_vol_delta( T, S, K, r, q, opt);

	// double Nd1 = NORM_DIST
	return NA_REAL;
}


double bs_gamma( double T, double S, double K, double vol, double r, double q, char o)
{
	double d1 	= get_d1( T, S, K, vol, r, q );
    double g 	= exp( -q * T ) * NORM_DENS( d1 ) / ( S * vol * sqrt( T ));
	switch ( o )
	{
		case 'c':
		case 'p':
			return g;
		case 's':
			return 2 * g;
	}
	return NA_REAL;
}


double bs_theta( double T, double S, double K, double vol, double r, double q, char o )
{
	double d1, d2;
	get_d1_and_d2( T, S, K, vol, r, q, &d1, &d2 );
	S *= exp( -q * T );
	K *= exp( -r * T );

	double theta;

	switch ( o )
	{
	case 'c':
		theta = -S * NORM_DENS( d1 ) * vol / 2 / sqrt( T );
		theta -= r * K * NORM_DIST( d2 );
		theta += q * S * NORM_DIST( d1 );
		return theta;
	case 'p':
		return NA_REAL;
	case 's':
		return NA_REAL;
	}


	return NA_REAL;
}


double bs_vega( double T, double S, double K, double vol, double r, double q, char o )
{
	if ( vol < 0 )
	{
		return -bs_vega( T, S, K, -vol, r, q, o );
	}

	if ( vol == 0 || T == 0 || S == 0 || K == 0 )
	{
		return 0;
	}

	double d1 	= get_d1( T, S, K, vol, r, q );
	double v 	= exp( -q * T ) * S * NORM_DENS( d1 ) * sqrt( T );

	switch ( o )
	{
	case 'c':
	case 'p':
		return v;
	case 's':
		return 2 * v;
	}

	return NA_REAL;
}

#define N_NUMARGS 6
static void bs_chkargs(SEXP *tau, SEXP *S, SEXP *K,
                       SEXP *vol, SEXP *r, SEXP *q,
                       SEXP *opt, SEXP *ret, int  *numprot)
{
    int i, N;

    SEXP *args[ N_NUMARGS + 1 ]  = { tau, S, K, vol, r, q, opt };
    char *names[ N_NUMARGS + 1 ] = { "tau", "S", "K", "vol", "r", "q", "opt" };

    /* Coerce types  */
    for ( i = 0; i < N_NUMARGS; i++ )
        ENSURE_NUMERIC( *args[i], *numprot );

    ENSURE_CHAR( *opt, *numprot );

    if ( length( *ret ) && !isNull( *ret ))
        ENSURE_CHAR( *ret, *numprot );

    /* Find maximal length */
    for (i = N = 0; i < N_NUMARGS + 1; i++)
        N = MAX(N, length(*args[i]));

    if ( N > 1 )
    { /* Repeat arguments of length == 1. */

        for ( i = 0; i < N_NUMARGS; i++ )
            if ( 1 == length( *args[ i ] ))
                PROT2( *args[ i ] = numrep( *args[ i ], N ), *numprot );

        if ( 1 == length( *opt ))
            PROT2( *opt = charrep( *opt, N ), *numprot );
    }

    /* Check lengths */
    for ( i = 0; i < N_NUMARGS + 1; i++ )
        if ( N != length( *args[ i ] ))
            error( "Argument %s of wrong length", names[ i ] );
}
#undef N_NUMARGS







int is_valid_option_type(char o)
{
	return IS_VALID_OPTION_TYPE( o );
}



SEXP bs(SEXP tau, SEXP S, SEXP K, SEXP vol, SEXP r, SEXP q, SEXP opt, SEXP ret)
/*******************************************************************
 *
 *  Description: Returns SEXP of type matrix with columns
 *    corresponding to the BS price, delta & vega.
 *
 *  Assumes tau is in CALENDAR DAYS
 *
 *******************************************************************/
{
    int numprot = 0;
    long i, j, m, retlen = 0;
    SEXP ans, dimnames;
    carray ans_matrix;

    bs_chkargs( &tau, &S, &K, &vol, &r, &q, &opt, &ret, &numprot );

/* By now all arguments should have same length */
    m      = length( tau );
    retlen = length( ret );

    ans_matrix = make_zero_matrix(m, retlen);

    for (j = 0; j < retlen; j++)
        for (i = 0; i < m; i++)
            bs_value(
				REAL( tau )[ i ] / 365.0, REAL( S )[ i ], REAL( K )[ i ],
				REAL( vol )[ i ], REAL( r )[ i ], REAL( q )[ i ],
                *CHAR( CHAR_PTR( opt )[ i ] ),
                *CHAR( CHAR_PTR( ret )[ j ] ),
                ARRAY2( ans_matrix )[ i ] + j
			);

    PROT2( ans = carray_to_sexp( ans_matrix ), numprot );
    PROT2( dimnames = NEW_LIST( 2 ), numprot );

    SET_ELT( dimnames, 0, GET_NAMES( tau ));
    SET_ELT( dimnames, 1, ret );

    setAttrib( ans, R_DimNamesSymbol, dimnames );

    UNPROTECT( numprot );

    return ans;
}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Functions for getting the strike corresponding to a certain delta.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


static void dtoK(
	const double T,
	const double S,
	const double vol,
	const double del,
	const char opt,
	const double r,
	const double q,
	double *ans
)
/*******************************************************************
 *
 *  Assumes T is in years.
 *
 *******************************************************************/
{
    register double t1, t2;

    if (T < 0 || S < 0 || vol < 0)
    {
        *ans = NA_REAL;
        return;
    }

    switch(opt)
    {
        case 'c':
            t1   = -vol * sqrt(T) * qnorm(exp(q * T) * del,0, 1, 1, 0);
            t2   = (r - q + 0.5 * SQR(vol)) * T;
            *ans = S * exp(t1 + t2);
            break;

        case 'p':
            t1   = -vol * sqrt(T) * qnorm(1.0 + exp(q * T) * del, 0, 1, 1, 0);
            t2   = (r - q + 0.5 * SQR(vol)) * T;
            *ans = S * exp(t1 + t2);
            break;

        case 's':
            t1   = -vol * sqrt(T) * qnorm(0.5 + 0.5 * exp(q * T) * del, 0, 1, 1, 0);
            t2   = (r - q + 0.5 * SQR(vol)) * T;
            *ans = S * exp(t1 + t2);
            break;

        default:
            *ans = NA_REAL;
            break;
    }
}



#define NARGS 7
static void deltaToStrike_chk(
	SEXP *T, SEXP *S, SEXP *vol, SEXP *del, SEXP *optionType, SEXP *r, SEXP *q, int *numprot
)
{
    SEXP *args[NARGS] = { T, S, vol, del, optionType, r, q };
    const SEXPTYPE type[ NARGS ] = { REALSXP, REALSXP, REALSXP, REALSXP, STRSXP, REALSXP, REALSXP };
    long lens[ NARGS ] = { 0 };
    long maxLen;
    int i;

    for ( i = 0; i < NARGS; i++ )
    {
        PROT2( *args[ i ] = coerceVector( *args[ i ], type[ i ]), *numprot );
        lens[ i ] = length( *args[ i ] );
    }

    maxLen = lmax( lens, NARGS );

    if ( maxLen > 1 )
        for ( i = 0; i < NARGS; i++ )
        {
            if ( 1 == lens[ i ] )
                PROT2( *args[ i ] = rep( *args[ i ], type[ i ], maxLen), *numprot);
            else if ( lens[ i ] != maxLen )
                error ("mismatched arg lengths");
        }
}


SEXP deltaToStrike(
	SEXP T, SEXP S, SEXP vol, SEXP del, SEXP optionType, SEXP r, SEXP q
)
{
    int numprot = 0;
    long len, i;
    SEXP K;

    deltaToStrike_chk( &T, &S, &vol, &del, &optionType, &r, &q, &numprot );

    len = length( T);
    PROT2( K = NEW_NUMERIC( len), numprot );

    for ( i = 0; i < len; i++ )
        dtoK(
			REAL( T )[ i ] / 365.0, REAL( S )[ i ], REAL( vol )[ i ],
            REAL( del )[ i ], *CHAR( GET_ELT( optionType, i )),
            REAL( r )[ i ], REAL( q )[ i ],
            REAL(K) + i );

    UNPROTECT(numprot);

    return K;
}

#undef NARGS
#undef IS_VALID_OPTION_TYPE
