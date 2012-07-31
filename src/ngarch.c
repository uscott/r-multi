#include "uri.h"
#include <time.h>

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
** This section treats the basic NGARCH(1, 1) model.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



#define NUM_NG11_PAR 6
#define LAMBDA_INDEX 0
#define A0_INDEX 1
#define A1_INDEX 2
#define B1_INDEX 3
#define GAMMA_INDEX 4
#define H1_INDEX 5


static void sim_ng11_path(const long T, const double *par, double *x0,
                          double *h0, double *z, double *x, double *h)
{
    const double *Lh = NULL;
    double const *hE = h + T;

    double z0;
    double lam, a0, a1, b1, gam;
    double a[ 1 ] = { 0.0 };
    double b[ 1 ] = { 0.0 };

    int ok;

    ok = par && x && h;

    if ( !ok )
        error( "sim_ng11_path: Null pointer passed where not allowed" );

    lam = par[ LAMBDA_INDEX ];
    a0  = par[ A0_INDEX ];
    a1  = par[ A1_INDEX ];
    b1  = par[ B1_INDEX ];
    gam = par[ GAMMA_INDEX ];

    ok = x0 && h0 && R_FINITE( *x0 ) && R_FINITE( *h0 );

    if ( !ok ) 
    {
        x0 = a;
        h0 = b;

        *h0 = a0 / (1 - b1 - a1 * (1 + SQR(gam)));
        *x0 = -0.5 * (*h0) + lam * sqrt(*h0);

        warning( "no starting points given" );
    }

    GetRNGstate();

    z0  = *x0 - ( -0.5 * ( *h0 ) + lam * sqrt( *h0 ));
    z0 /= sqrt( *h0 );
    
    *h = a0 + ( *h0 ) * ( b1 + a1 * SQR( z0 - gam ));
    *x = -0.5 * ( *h ) + sqrt( *h ) * ( lam + *z );
    
    h  += 1;
    x  += 1;
    Lh  = h - 1;
  
    for ( ; h < hE && R_FINITE(*h); x++, h++, Lh++ ) 
    {
        *h  = a0 + ( *Lh ) * ( b1 + a1 * SQR( *z - gam ));
         z += 1;
        *x  = -0.5 * ( *h ) + sqrt( *h ) * ( lam + *z );
    }

    PutRNGstate();
}


static void ng11(const long T, const double *x, const double *par,
                 int llonly, double *h, double *z, double *ll)

/********************************************************************
 *
 *   Description: Evaluation of loglikelihood, conditional variance
 *    and residuals for the given time series x[0], ..., x[T - 1]
 *    and the given coefficients par[0], ..., par[NUM_NG11_PAR-1].
 *
 ********************************************************************/
{
    int inc = 0;

    const double *Lh = NULL;
    double const *xT = x + T;
    double a[1]      = {0.0};
    double b[1]      = {0.0};
    double        
        lam, a0, a1, b1, gam, h1;
    int           
        ok;

    inc = !( llonly = llonly || !h || !z );

    llonly ? h = (double *) a, z = (double *) b : 0;

    ok = x && par && ll;

    if ( !ok )
        error ( "null pointer passed where not allowed" );

    lam = par[ LAMBDA_INDEX ];
    a0  = par[ A0_INDEX ];
    a1  = par[ A1_INDEX ];
    b1  = par[ B1_INDEX ];
    gam = par[ GAMMA_INDEX ];
    h1  = par[ H1_INDEX ];

    ok = 
        R_FINITE( T )  && R_FINITE( lam ) &&
        R_FINITE( a0 ) && R_FINITE( b1 )  &&
        R_FINITE( a1 ) && R_FINITE( gam );

    if ( !ok ) 
    {
        *ll = HUGE_NEGATIVE;
        warning("improper arguments, assuming loglikelihood = %.0f", HUGE_NEGATIVE);

        return;
    }

    if ( !R_FINITE( h1 )) 
    {
        warning( "Assuming h[1] = E(h)" );
        h1 = a0 / (1.0 - b1 - a1 * ( 1.0 + SQR( gam )));
    }

    *h  = h1;
    *z  = *x / sqrt( *h ) - lam + 0.5 * sqrt( *h );
    *ll = T * log( M_2PI ) + log( *h ) + SQR( *z );

    Lh  = h;
    x  += 1;
    h  += inc;
    
    for ( ; x < xT && R_FINITE( *ll ); x++, h += inc, Lh += inc) 
    {
        *h   = a0 + Lh[ 0 ] * ( b1 + a1 * SQR( *z - gam ));
         z  += inc;
        *z   = *x / sqrt( *h ) - lam + 0.5 * sqrt( *h );
        *ll += log( *h ) + SQR( *z );
    }

    *ll = R_FINITE( *ll ) ? -0.5 * ( *ll ) : HUGE_NEGATIVE;
}



static void ngarch11_chk(SEXP *x, SEXP *par, double *parArray, int *numprot)
{
    ENSURE_NUMERIC( *x, *numprot );

    if ( isNewList( *par ) ) 
    {
        parArray[ LAMBDA_INDEX ] = asReal( getListElt( *par, "lambda" ));
        parArray[ A0_INDEX ]     = asReal( getListElt( *par, "a0" ));
        parArray[ A1_INDEX ]     = asReal( getListElt( *par, "a1" ));
        parArray[ B1_INDEX ]     = asReal( getListElt( *par, "b1" ));
        parArray[ GAMMA_INDEX ]  = asReal( getListElt( *par, "gamma" ));
        parArray[ H1_INDEX ]     = asReal( getListElt( *par, "h1" ));
    } 
    else 
    {
        ENSURE_NUMERIC( *par, *numprot );
        ENSURE_LENGTH(  *par, NUM_NG11_PAR, *numprot );

        parArray[ LAMBDA_INDEX ] = REAL( *par )[ LAMBDA_INDEX ];
        parArray[ A0_INDEX ]     = REAL( *par )[ A0_INDEX ];
        parArray[ A1_INDEX ]     = REAL( *par )[ A1_INDEX ];
        parArray[ B1_INDEX ]     = REAL( *par )[ B1_INDEX ];
        parArray[ GAMMA_INDEX ]  = REAL( *par )[ GAMMA_INDEX ];
        parArray[ H1_INDEX ]     = REAL( *par )[ H1_INDEX ];
    }

    if ( SET_NA_TO_ZERO( parArray[ LAMBDA_INDEX ] )) 
        warning( "Assuming lambda = 0" );
    if ( SET_NA_TO_ZERO( parArray[ A1_INDEX ] )) 
        warning( "Assuming a1 = 0" );
    if ( SET_NA_TO_ZERO( parArray[ B1_INDEX ] )) 
        warning( "Assuming b1 = 0" );
    if ( SET_NA_TO_ZERO( parArray[ GAMMA_INDEX ] )) 
        warning( "Assuming gamma = 0" );
}



#define ANS_LEN 3
SEXP ngarch11(SEXP x, SEXP par, SEXP ll_only)
{
    int 
        numprot = 0;
    long 
        xLen, len;
    double 
        ll;
    char 
        *names[ ANS_LEN ]        = { "ll", "h", "res" };
    const int 
        llonly                   = asInteger( ll_only );
    double 
        parArray[ NUM_NG11_PAR ] = { 0 };

    SEXP ans, h, res;

    ngarch11_chk( &x, &par, parArray, &numprot );

    xLen = length(x);
    len  = !llonly * xLen;

    PROT2( h   = NEW_NUMERIC( len ),  numprot );
    PROT2( res = NEW_NUMERIC( len ),  numprot );
    PROT2( ans = NEW_LIST( ANS_LEN ), numprot );

    if ( len == xLen ) 
    {
        SET_NAMES( h,   GET_NAMES( x ));
        SET_NAMES( res, GET_NAMES( x ));
    }

    ng11( xLen, REAL( x ), parArray, llonly, REAL( h ), REAL( res ), &ll );

    SET_ELT( ans, 0, ScalarReal( ll ));
    SET_ELT( ans, 1, h );
    SET_ELT( ans, 2, res );

    set_names( ans, names );

    UNPROTECT(numprot);

    return ans;
}
#undef ANS_LEN



static double ll_ngarch11(SEXP par, SEXP x)
/*******************************************************************
 *
 *  Description: Used in optimization to fit NGARCH(1, 1) model.
 *    Modification of function ngarch11.
 *
 *******************************************************************/
{
    int numprot = 0;
    int ok;
    double ll = 0.0;

    ok = IS_NUMERIC( x ) && IS_NUMERIC( par );

    if ( !ok ) error( "bad args to ll_ngarch11" );

    ENSURE_LENGTH( par, NUM_NG11_PAR, numprot );
  
    ng11( length( x ), REAL( x ), REAL( par ), 
          /* llonly = */ 1, 
          /* h = */ (double *) NULL, 
          /* z = */ (double *) NULL, 
          &ll);

    UNPROTECT(numprot);

    return ll;
}



static void adjustInitPar(SEXP *initPar, int parLen, int *numprot)
{
    SEXP par = R_NilValue;
  
    if ( isNewList( *initPar )) 
    {
        PROT2( par = NEW_NUMERIC( parLen ), *numprot );

        REAL( par )[ LAMBDA_INDEX ] = asReal( getListElt( *initPar, "lambda" ));
        REAL( par )[ A0_INDEX ]     = asReal( getListElt( *initPar, "a0" ));
        REAL( par )[ A1_INDEX ]     = asReal( getListElt( *initPar, "a1" ));
        REAL( par )[ B1_INDEX ]     = asReal( getListElt( *initPar, "b1" ));
        REAL( par )[ GAMMA_INDEX ]  = asReal( getListElt( *initPar, "gamma" ));

        if ( H1_INDEX < parLen )
            REAL( par )[ H1_INDEX ] = asReal( getListElt( *initPar, "h1" ));

        SET_NA_TO_ZERO( REAL( par )[ LAMBDA_INDEX ]);
        SET_NA_TO_ZERO( REAL( par )[ A1_INDEX ]);    
        SET_NA_TO_ZERO( REAL( par )[ B1_INDEX ]);    
        SET_NA_TO_ZERO( REAL( par )[ GAMMA_INDEX ]);

        *initPar = par;
    }
}



#define ANS_LEN 7
SEXP fit_ngarch11( SEXP x, SEXP initPar, SEXP fitInit, 
                   SEXP basePopSize, SEXP tol, 
                   SEXP stopLags, SEXP minit, SEXP maxit, SEXP options)
{
    const int    parLen       = NUM_NG11_PAR-1 + !!asInteger( fitInit );
    const int    ibasePopSize = asInteger( basePopSize );
    const long   istopLags    = asInteger( stopLags );
    const long   imaxit       = asInteger( maxit );
    const long   iminit       = asInteger( minit );
    const double dtol         = asReal( tol );

    int  numprot = 0;
    int  gradToo;
    long gradMaxIt;
    char *names[ ANS_LEN ] = 
        { "par", "ll", "convergence", "x", "h", "res", "gradconv" };
    SEXP ans, fit1, fit2, tmp, dimNames, parNames;

    ENSURE_NUMERIC( x, numprot );

    adjustInitPar( &initPar, parLen, &numprot );

    /* Do the fitting! */
    fit1 = optimGa2( ll_ngarch11, initPar, x, parLen, 
                     ibasePopSize, istopLags, iminit, imaxit, 
                     0, dtol, 0);

    PROT2( fit1, numprot );
    PROT2( ans      = NEW_LIST( ANS_LEN ), numprot );
    PROT2( parNames = NEW_STRING( parLen ), numprot );
    PROT2( dimNames = NEW_LIST( 2 ), numprot );

    gradToo = asInteger( getListElt( options, "grad" ));
    gradToo = SET_NA_TO_ZERO( gradToo );

    if ( gradToo ) 
    {
        gradMaxIt = asInteger( getListElt( options, "maxit" ));
    
        if ( !R_FINITE( (double)gradMaxIt ) || gradMaxIt <= 0)
            gradMaxIt = imaxit;

        fit2 = optimGradient2( ll_ngarch11, 
                               getListElt( fit1, "par" ), 
                               x, dtol, 0, 0, gradMaxIt );
        PROT2( fit2, numprot );
    } 
    else 
    {
        fit2 = fit1;
    }

    PROT2( tmp = ngarch11( x, getListElt( fit2, "par" ), ScalarInteger( 0 )), 
           numprot);

    SET_ELT( ans, 0, getListElt( fit2, names[ 0 ] ));
    SET_ELT( ans, 1, getListElt( fit2, "value" ));
    SET_ELT( ans, 2, getListElt( fit1, names[2] ));
    SET_ELT( ans, 3, x );
    SET_ELT( ans, 4, getListElt( tmp, names[4] ));
    SET_ELT( ans, 5, getListElt( tmp, names[5] ));
    SET_ELT( ans, 6, getListElt( fit2, "convergence" ));

    set_names( ans, names );

    CHAR_PTR( parNames )[ LAMBDA_INDEX ] = mkChar( "lambda" );
    CHAR_PTR( parNames )[ A0_INDEX ]     = mkChar( "a0" );
    CHAR_PTR( parNames )[ A1_INDEX ]     = mkChar( "a1" );
    CHAR_PTR( parNames )[ B1_INDEX ]     = mkChar( "b1" );
    CHAR_PTR( parNames )[ GAMMA_INDEX ]  = mkChar( "gamma" );

    if ( NUM_NG11_PAR == parLen )
        CHAR_PTR( parNames )[ H1_INDEX ] = mkChar( "h1" );

    SET_ELT( dimNames, 0, parNames );
    SET_ELT( dimNames, 1, R_NilValue );

    setAttrib( GET_ELT(ans, 0), R_DimNamesSymbol, dimNames );

    UNPROTECT(numprot);

    return ans;
}
#undef ANS_LEN


#define ANS_LEN 11
SEXP BootfitNgarch11( SEXP x, SEXP initPar, SEXP fitInit, 
                      SEXP basePopSize, SEXP tol, 
                      SEXP stopLags, SEXP minit, SEXP maxit, SEXP options,
                      SEXP numBoots )
{
    SEXP 
        baseFit = R_NilValue, thisFit = R_NilValue, res, thisRes,
        bootPars, thisPars, x_new, ans, bootRes, bootPaths, bootVar,
        parNames, dimNames;
    int nboots = asInteger( numBoots ), numprot = 0, i, xlen;
    double *z, *pars, *thisPath, *thisVar;
    const int parLen = NUM_NG11_PAR-1 + !!asInteger( fitInit );
    char *names[ ANS_LEN ] = 
        { "par", "ll", "convergence", "x", "h", "res", "gradconv", "bootPar", 
          "bootPaths", "bootVar", "bootRes" };

    ENSURE_NUMERIC( x, numprot );
    
    xlen = length( x );
    adjustInitPar( &initPar, parLen, &numprot );
    
    PROT2( baseFit = fit_ngarch11( x, initPar, fitInit, basePopSize, tol,
                                   stopLags, minit, maxit, options ),
           numprot );

    PROT2( ans = NEW_LIST( ANS_LEN ), numprot );    
    PROT2( res = getListElt( baseFit, "res" ), numprot );
    PROT2( thisRes = NEW_NUMERIC( xlen ), numprot );    
    PROT2( bootPars = allocMatrix( REALSXP, parLen, nboots + 1 ), numprot );
    PROT2( bootPaths = allocMatrix( REALSXP, xlen, nboots + 1 ), numprot );
    PROT2( bootVar = allocMatrix( REALSXP, xlen, nboots + 1 ), numprot );    
    PROT2( bootRes = allocMatrix( REALSXP, xlen, nboots + 1 ), numprot );    
    PROT2( x_new = NEW_NUMERIC( xlen ), numprot );
    PROT2( parNames = NEW_STRING( parLen ), numprot );
    PROT2( dimNames = NEW_LIST( 2 ), numprot );
    
    z = REAL( thisRes );    
    pars = matcol1( bootPars, 0 );    
    
    memcpy( pars, REAL( getListElt( baseFit, "par" ) ), parLen * sizeof( double ));    
    memcpy( matcol1( bootPaths, 0 ), REAL( x ), xlen * sizeof( double ));
    memcpy( matcol1( bootVar, 0 ), REAL( getListElt( baseFit, "h" )), xlen * sizeof( double ));    
    memcpy( matcol1( bootRes, 0 ), REAL( getListElt( baseFit, "res" )), xlen * sizeof( double ));    

    thisPath = (double *) malloc( xlen * sizeof( double ));    
    thisVar  = (double *) malloc( xlen * sizeof( double ));
    
    for ( i = 1; i <= nboots; ++i )
    {        
        pars = matcol1( bootPars, i );
        
        get_random_sample( REAL( res ), z, xlen, xlen );
        
        sim_ng11_path( xlen, matcol1( bootPars, 0 ), NULL, NULL, z, thisPath, thisVar );
     
        memcpy( REAL( x_new ), thisPath, xlen * sizeof( double ));
        
        PROTECT( thisFit = fit_ngarch11( x_new, initPar, fitInit, basePopSize, tol,
                                         stopLags, minit, maxit, options ));
        
        
        memcpy( pars, REAL( getListElt( thisFit, "par" )), parLen * sizeof( double ));
        memcpy( matcol1( bootPaths, i ), thisPath, xlen * sizeof( double ));
        memcpy( matcol1( bootVar, i ), thisVar, xlen * sizeof( double ));
        memcpy( matcol1( bootRes, i ), z, xlen * sizeof( double ));
        

        UNPROTECT( 1 );        
    }

    CHAR_PTR( parNames )[ LAMBDA_INDEX ] = mkChar( "lambda" );
    CHAR_PTR( parNames )[ A0_INDEX ]     = mkChar( "a0" );
    CHAR_PTR( parNames )[ A1_INDEX ]     = mkChar( "a1" );
    CHAR_PTR( parNames )[ B1_INDEX ]     = mkChar( "b1" );
    CHAR_PTR( parNames )[ GAMMA_INDEX ]  = mkChar( "gamma" );

    if ( NUM_NG11_PAR == parLen )
        CHAR_PTR( parNames )[ H1_INDEX ] = mkChar( "h1" );

    SET_ELT( dimNames, 0, parNames );
    SET_ELT( dimNames, 1, R_NilValue );

    setAttrib( bootPars, R_DimNamesSymbol, dimNames );
    
    for ( i = 0; i < ANS_LEN - 4; ++i )
        SET_ELT( ans, i, getListElt( baseFit, names[ i ] ));

    SET_ELT( ans, ANS_LEN - 4, bootPars );
    SET_ELT( ans, ANS_LEN - 3, bootPaths );
    SET_ELT( ans, ANS_LEN - 2, bootVar );
    SET_ELT( ans, ANS_LEN - 1, bootRes );    
    
    set_names( ans, names );
    
    free( thisPath );
    free( thisVar );
    
    UNPROTECT( numprot );    

    return ans;    
}
#undef ANS_LEN



#define NPAR  (NUM_NG11_PAR-1)
#define NARGS 6
static void sim_ngarch11_chk(SEXP *T, SEXP *npaths, SEXP *par,
                             SEXP *x0, SEXP *h0, SEXP *z,
                             int  *numprot)
{
    int ok, i;

    SEXP *args[  NARGS ] = { T, npaths, par, x0, h0, z };
    char *names[ NARGS ] = { "T", "npaths", "par", "x0", "h0", "z" };

    for ( i = 0; i < NARGS - 3; i++ )
        if ( isNull( *args[ i ] ))
            error( "Argument %s not allowed to be NULL", names[ i ] );

    ok = asInteger( *T ) > 0 && asInteger( *npaths ) > 0;

    if ( !ok )
        error( "Non-positive number of paths or path length." );

    if ( NPAR != length( *par ))
        error( "Wrong number of parameters passed" );

    for ( i = 2; i < NARGS; i++ )
        if ( !IS_NUMERIC( *args[ i ] ))
            PROT2( *args[ i ] = AS_NUMERIC( *args[ i ] ), *numprot );
}
#undef NPAR
#undef NUM_NG11_PAR



SEXP sim_ngarch11(SEXP T, SEXP npaths, SEXP par,
                  SEXP x0, SEXP h0, SEXP z)
{
    int zlen; 
    int numprot = 0;
    long pathLen, numPaths;
    long j;
    double 
        *uPtr                = NULL, 
        *xPtr                = NULL, 
        *hPtr                = NULL;
    char *names[ 2 + NARGS ] = { "x", "h", "T", "npaths", "par", "x0", "h0", "z" };
    SEXP ans, x, h;

    sim_ngarch11_chk(&T, &npaths, &par, &x0, &h0, &z, &numprot);

    pathLen  = (long) asInteger( T );
    numPaths = (long) asInteger( npaths );

    zlen = length(z);
 
    PROT2( ans = NEW_LIST( 2 + NARGS ), numprot );

    PROT2( x = allocMatrix( REALSXP, pathLen, numPaths ), numprot );
    PROT2( h = allocMatrix( REALSXP, pathLen, numPaths ), numprot );

    GetRNGstate();
 
    uPtr = (double *) R_alloc(pathLen, sizeof(double));

    for ( j = 0; j < numPaths; j++ ) 
    {
        xPtr = matcol1( x, j );
        hPtr = matcol1( h, j );

        get_random_sample2( REAL( z ), uPtr, zlen, pathLen );

        sim_ng11_path( pathLen, REAL( par ), REAL( x0 ),
                       REAL(h0), uPtr, xPtr,
                       hPtr );
    }
 
    PutRNGstate();

    SET_ELT(ans, 0, x);
    SET_ELT(ans, 1, h);
    SET_ELT(ans, 2, T);
    SET_ELT(ans, 3, npaths);
    SET_ELT(ans, 4, par);
    SET_ELT(ans, 5, x0);
    SET_ELT(ans, 6, h0);
    SET_ELT(ans, 7, z);

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}
#undef NARGS

#undef LAMBDA_INDEX 
#undef A0_INDEX 
#undef A1_INDEX 
#undef B1_INDEX 
#undef GAMMA_INDEX 
#undef H1_INDEX
