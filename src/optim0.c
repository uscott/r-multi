
#include "uri.h"

#include <time.h>


#define THREE 3
#define BIG_POP_SIZE( n ) (2 * (n) + (2 * (n) - 1) * (n))
#define C2( n ) (((n) * ((n) - 1)) / 2)



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  This section contains the code for the gradient-descent
**  optimizer to called inside other C functions.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



#define NUM_TRIES 2
static int NumDeriv( OptimFunc f, 
                     double *par, 
                     long parlen, 
                     void *context, 
                     double *gradient,
                     double *eps,
                     double *par_new)

/*******************************************************************
 *
 *  Description: Numerically estimates gradient vector at point par
 *    of function par -> f( par, parlen, context )
 *
 *******************************************************************/
{
    int /* ok = 0, */ try_num = 0;
    long i;
    const long len = parlen;    
    double val = 0.0, val_new = 0.0;

    val = f( par, len, context );    

    if ( !R_FINITE( val )) return 0;
  
    dsetna( gradient, len ); /* Initialize elements of gradient to NA */
    memcpy( par_new, par, len * sizeof(double));

    GetRNGstate();

    for ( i = 0; i < len; ++i )
    {
        if (0.0 == ABS( eps[ i ] ) || !R_FINITE( eps[ i ] ))
            eps[ i ] = sqrt( DOUBLE_EPS );
        
        if ( unif_rand() < 0.5 )
            eps[ i ] *= -1.0;

        do 
        {
            par_new[ i ] = par[ i ] + eps[ i ];
            eps[ i ]    *= runif( 0.9, 1 );
            try_num     += 1;

            val_new = f( par_new, len, context );

        } while ( !R_FINITE( val_new ) && try_num < NUM_TRIES );

        if ( !R_FINITE( val_new ))
            return 0;

        gradient[ i ] = ( val_new - val ) / eps[ i ];
        par_new[ i ] = par[ i ];
        try_num = 0;

        if ( 0.0 == gradient[ i ] ) /* perturbation might be too small */
            eps[ i ] *= runif(1, 2);

    }

    PutRNGstate();

    return 1;

}
#undef NUM_TRIES



#define OPTIM_ANS_LEN 3
#define OPTIM_ANS_NAMES { "par", "value", "convergence" }

SEXP OptimGradient00( OptimFunc f, 
                      double *initPar, 
                      long parlen,
                      void *context,
                      const double tol, 
                      const int relTol, 
                      const int minimize, 
                      const long maxit)

  /*******************************************************************
   *
   *  Description: Basic gradient-descent optimizer intended for
   *    calling from inside other C functions and not from R.
   *    The function being optimized is x |-> f(x, controlPar).
   *
   *  Status: Not finished.
   *
   *******************************************************************/
{	
    const double sgn  = minimize ? 1.0 : -1.0;

    int numprot = 0;

    long i, it;
    int convergence = 1;
    double stepsize = 1e-1;
    double val1, val2, tolMult;
    double delta;

    char *ansNames[ OPTIM_ANS_LEN ] = OPTIM_ANS_NAMES;

    SEXP ans = R_NilValue /* par1 = R_NilValue, par2 = R_NilValue, 
        gradient = R_NilValue, eps = R_NilValue, 
        *x = &par1, *xNew = &par2, *tmp = NULL */ ;

    SEXP tmp;

    double *par1, *par2, *eps, *gradient, **p1, **p2, **p3;

    PROT2( ans      = NEW_LIST( OPTIM_ANS_LEN ), numprot );
    set_names( ans, ansNames );

    gradient = (double *) R_alloc( parlen, sizeof( double ));    
    par1     = (double *) R_alloc( parlen, sizeof( double ));
    par2     = (double *) R_alloc( parlen, sizeof( double ));
    eps      = (double *) R_alloc( parlen, sizeof( double ));

    p1       = &par1;
    p2       = &par2;
    
    PROT2( tmp = NEW_NUMERIC( parlen ), numprot );    

    val1 = f( initPar, parlen, context );

    if ( !R_FINITE( val1 ))
    {
        memcpy( REAL( tmp ), initPar, parlen * sizeof( double ));
        
        SET_ELT( ans, 0, tmp );
        SET_ELT( ans, 1, ScalarReal( val1 ));
        SET_ELT( ans, 2, ScalarInteger( convergence ));

        UNPROTECT( numprot );

        return ans;
    }

    memcpy( par1, initPar, parlen * sizeof(double));
    memcpy( par2, initPar, parlen * sizeof(double));

    for ( i = 0; i < parlen; ++i )
        eps[ i ] = sqrt( DOUBLE_EPS );

    GetRNGstate();

    for ( it = 0, convergence = 1; it < maxit && convergence; ++it ) 
    {
        if ( !NumDeriv( f, 
                        *p1,
                        parlen, 
                        context, 
                        gradient, 
                        eps, 
                        *p2 ))
            continue;

        stepsize *= runif( 1, 1.1 );
        
        for( i = 0; i < parlen; ++i )
        {
            delta      = -sgn * gradient[ i ] * stepsize;
            (*p2)[ i ] = (*p1)[ i ] + delta;            
        }
        
        val2 = f( *p2, parlen, context );

        if ( R_FINITE( val2 ) && sgn * val2 < sgn * val1 ) 
        {
            tolMult = relTol ? tol + ABS( val1 ) : 1.0;

            convergence = ( ABS( val2 - val1 ) >= tol * tolMult );

            val1 = val2;

            p3 = p1;
            p1 = p2;
            p2 = p3;
        }
        else
            stepsize /= runif( 1, 3 );
    }

    PutRNGstate();

    memcpy( REAL( tmp ), *p1, parlen * sizeof( double ));    

    SET_ELT(ans, 0, tmp );
    SET_ELT(ans, 1, ScalarReal( val1 ));
    SET_ELT(ans, 2, ScalarInteger( convergence ));

    UNPROTECT(numprot);

    return ans;

}



SEXP OptimGradient01( OptimFunc f, 
                      SEXP  initpar, 
                      long  parlen,
                      void *context,
                      const double tol, 
                      const int relTol, 
                      const int minimize, 
                      const long maxit)
{  
    int 
        numprot = 0;
    long 
        j, m, n;
    char 
        *ans_names[ OPTIM_ANS_LEN ] = OPTIM_ANS_NAMES;
    SEXP 
        par  = R_NilValue, val  = R_NilValue, conv = R_NilValue,
        ans  = R_NilValue, tmp2 = R_NilValue, tmp3 = R_NilValue;
    double
        *tmp1;    
  
    PROT2( ans = NEW_LIST( OPTIM_ANS_LEN ), numprot );
    set_names( ans, ans_names );

    if ( IS_MATRIX( initpar )) 
    {
        m = nrows( initpar );
        n = ncols( initpar );
    
        PROT2( par  = allocMatrix( REALSXP, m, n ), numprot );
        PROT2( val  = NEW_NUMERIC( n ), numprot );
        PROT2( conv = NEW_INTEGER( n ), numprot );

        tmp1 = (double *) malloc( m * sizeof( double ));        

        for ( j = 0; j < n; ++j ) 
        {
            memcpy( tmp1, matcol1( initpar, j ), m * sizeof( double ));
      
            tmp2 = OptimGradient00( f, 
                                    tmp1,
                                    m,
                                    context,
                                    tol, relTol, minimize, maxit );

            PROTECT( tmp2 );

            PROTECT( tmp3 = getListElt( tmp2, "par" ));

            if ( IS_NUMERIC( tmp3 ) && length( tmp3 ) == m )
                memcpy( matcol1( par, j ), REAL( tmp3 ), m * sizeof(double));
            else
                dsetna( matcol1( par, j ), m );
      
            REAL( val )[ j ]     = asReal( getListElt( tmp2, "value" ));
            INTEGER( conv )[ j ] = asInteger(getListElt(tmp2, "convergence" ));

            UNPROTECT( 2 );
        }

        free( tmp1 );

        SET_ELT(ans, 0, par);
        SET_ELT(ans, 1, val);
        SET_ELT(ans, 2, conv);
    } 
    else if ( isNewList( initpar )) 
    {
        n = length( initpar );
        PROT2( par  = NEW_LIST( n ), numprot );
        PROT2( val  = NEW_NUMERIC( n ), numprot );
        PROT2( conv = NEW_INTEGER( n ), numprot );
        
        for ( j = 0; j < n; ++j ) 
        {
            m = length( GET_ELT( initpar, j ));
            tmp1 = (double *) malloc( m * sizeof( double ));
            memcpy( tmp1, 
                    REAL( GET_ELT( initpar, j )), 
                    m * sizeof( double ));            

            tmp2 = OptimGradient00( f, 
                                    tmp1, 
                                    m, 
                                    context, 
                                    tol, 
                                    relTol, 
                                    minimize, 
                                    maxit);
            PROTECT( tmp2 );

            SET_ELT( par, j, getListElt( tmp2, "par" ));
            REAL( val )[ j ] = asReal( getListElt( tmp2, "value" ));
            INTEGER( conv )[ j ] = asInteger( getListElt( tmp2,"convergence" ));

            UNPROTECT( 2 );
            
            free( tmp1 );            
        }

        SET_ELT(ans, 0, par);
        SET_ELT(ans, 1, val);
        SET_ELT(ans, 2, conv);
        
    } else 
    {
        ans = OptimGradient00( f, 
                               REAL( initpar ), 
                               parlen,
                               context, 
                               tol, 
                               relTol, 
                               minimize, 
                               maxit);
        PROT2(ans, numprot);
    }

    UNPROT2;

    return ans;

}



/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
**  Following functions are for genetic algo optimizer intended for
**  being called from inside C functions and not from R.
*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/



static void OptimGa0_check( SEXP *initpar, 
                            const int parlen, 
                            const int basePopSize, 
                            int *numprot )
{
    long i = 0;
    long initLen = 0;
    SEXP tmp;

    if ( basePopSize < 1 || parlen < 1 )
        error( "Must have population size > 0 and parameter length > 0");

    GetRNGstate();

    if( isNull( *initpar ) || !length( *initpar )) 
    {
        PROT2( *initpar = allocMatrix( REALSXP, parlen, basePopSize ), 
               *numprot);        

        for ( i = 0; i < parlen; ++i )
            REAL( *initpar )[ i ] = 1e-4;        

        for ( i = parlen; i < parlen * basePopSize; ++i )
            REAL( *initpar )[ i ] = unif_rand();
    }
    else if( isMatrix( *initpar )) 
    {
        if( parlen != nrows( *initpar ))
            error( "Inconsistent parameter length info" );

        PROT2( tmp = allocMatrix( REALSXP, parlen, basePopSize ), *numprot );

        memcpy( REAL( tmp ), 
                REAL( *initpar ), 
                parlen * basePopSize * sizeof( double ));

        for ( i = length( *initpar ); i < parlen * basePopSize; ++i )
            REAL( tmp )[ i ] = unif_rand();

        *initpar = tmp;
    }
    else 
    {
        PROT2( *initpar = AS_NUMERIC( *initpar ), *numprot );

        initLen = length( *initpar );

        if ( initLen != parlen )
            error ("Inconsistent parameter length info");

        PROT2( tmp = allocMatrix( REALSXP, parlen, basePopSize ), *numprot);
    
        memcpy( REAL( tmp ), 
                REAL( *initpar ), 
                parlen * basePopSize * sizeof( double ));

        for ( i = initLen; i < parlen * basePopSize; ++i )
            REAL( tmp )[ i ] = unif_rand();

        *initpar = tmp;    
    }

    PutRNGstate();
}



static void Mutate( SEXP pop, const int basePopSize )
  /*******************************************************************
   *
   *  Description: Mutates population of column vectors 0 throught
   *    basePopSize - 1 in par putting the results in column vectors
   *    basePopSize throught 2 * basePopSize - 1.
   *
   *******************************************************************/

{
    int ok, j;
    const int parlen = nrows( pop );
    double *ptr1 = NULL, *ptr2 = NULL;

    ok = isMatrix( pop ) && basePopSize > 0 && ncols( pop ) >= 2 * basePopSize;

    if (!ok)
        error ("bad argument passed to mutate");

    GetRNGstate();

    ptr1 = REAL( pop );
    ptr2 = REAL( pop ) + parlen * basePopSize;

    for (j = 0; j < parlen * basePopSize; j++, ptr1++, ptr2++)
        *ptr2 = ISNAN(*ptr1) ? norm_rand() : rnorm(*ptr1, 1+ABS(*ptr1));

    PutRNGstate();
}




static void Breed( SEXP pop, const int basePopSize )
  /*******************************************************************
   *
   *  Description: Breeds population of column vectors 0 through
   *
   *******************************************************************/
{
    const int parlen = nrows(pop);
    const int bigPopSize = ncols(pop);

    int counter;

    int i, j, k, ok;
    double w;
    double *ptr1, *ptr2, *kidPtr;

    ok = 
        isMatrix(pop) && parlen &&
        bigPopSize == BIG_POP_SIZE(basePopSize);

    if ( !ok )
        error ( "Bad arguments in breed" );

    GetRNGstate();

    kidPtr = REAL( pop ) + 2 * parlen * basePopSize;

    counter = 0;

    for ( i = 0; i < 2*basePopSize - 1; ++i )
        for ( j = i + 1; j < 2*basePopSize; ++j ) 
        {
            ptr1 = matcol1( pop, i );
            ptr2 = matcol1( pop, j );

            for(k = 0; k < parlen; ++k ) 
            {
                w = norm_rand() + 0.5;

                *kidPtr++ = w * (*ptr1++) + (1.0 - w) * (*ptr2++);
            }

            if ( ++counter >= bigPopSize )
                error ( "Breed ... something wrong here" );
        }

    PutRNGstate();
}



static void GetNextGen( OptimFunc f, 
                        long parlen,
                        void *context, 
                        SEXP pop,
                        SEXP fargs, 
                        SEXP fvals, 
                        SEXP ix,
                        SEXP tmp, 
                        const int basePopSize, 
                        const int minimize)
{
    const int bigPopSize = ncols( pop );
    /* const int parlen = nrows( pop ); */
    const double sgn = minimize ? 1.0 : -1.0;

    int ok, *ixPtr;
    long i, j;
    double val;

    ok = 
        isMatrix( pop ) && isMatrix( tmp ) && IS_INTEGER( ix ) &&
        IS_NUMERIC( fargs ) && IS_NUMERIC( fvals ) && 
        bigPopSize == length( ix ) &&
        bigPopSize == length( fvals ) && 
        bigPopSize == BIG_POP_SIZE( basePopSize ) && 
        parlen == length( fargs ) &&
        parlen == nrows( tmp ) &&
        basePopSize == ncols( tmp );

    if ( !ok )
        error ( "Bad arguments to GetNextGen" );
  
    for( i = 0; i < bigPopSize; ++i ) 
        INTEGER( ix )[ i ] = i;

    Mutate( pop, basePopSize );
    Breed ( pop, basePopSize );

    for ( i = 0; i < bigPopSize; i++ ) 
    {
        memcpy( REAL( fargs ), matcol1( pop, i ), parlen * sizeof( double ));

        val = f( REAL( fargs ), parlen, context );

        REAL( fvals )[ i ] = R_FINITE( val ) ? sgn * val : HUGE;
    }

    /* Move top basePopSize population members to front */
    rsort_with_index( REAL( fvals ), INTEGER( ix ), bigPopSize );
    
    ixPtr = INTEGER( ix );
		
    for ( j = 0; j < basePopSize; j++ )
        memcpy( matcol1( tmp, j ), 
                matcol1( pop, *ixPtr++ ), 
                parlen * sizeof( double ));

    memcpy( REAL( pop ), REAL( tmp ), parlen * basePopSize * sizeof(double));
    /* Done sorting population */

    /* Scale back 1st basePopSize function values */
    for ( i = 0; i < basePopSize; i++ )
        REAL( fvals )[ i ] *= sgn;
}



SEXP OptimGa0( OptimFunc f, 
               SEXP initPar,
               const long parlen,
               void *context,
               const int basePopSize, 
               const long stopLags, 
               const long minit, 
               const long maxit, 
               const int minimize, 
               const double tol, 
               const int relTol)

  /*******************************************************************
   *
   *  Description: Uses genetic algorithms with gradient descent
   *    at the end to optimize function x |-> f(x, controlPar).
   *
   *******************************************************************/
{
    int numprot = 0;  
    const int bigPopSize = BIG_POP_SIZE(basePopSize);
    int gen = 1, notConvergedYet = 1, numAnsCols = 1;
    double fmin, fmax;
    long lagNum = 0;
    char *names[] = OPTIM_ANS_NAMES;
    int *dimPtr = NULL;

    SEXP fargs, ans, pop, fvals, tmp, ix, popDim;

    /* Check arguments and perform necessary adjustments */
    OptimGa0_check( &initPar, parlen, basePopSize, &numprot );

    PROT2(fvals = NEW_NUMERIC( bigPopSize ), numprot);
    PROT2(fargs = NEW_NUMERIC( parlen ),     numprot);
    PROT2(ix    = NEW_INTEGER( bigPopSize ), numprot);
    PROT2(ans   = NEW_LIST( OPTIM_ANS_LEN ), numprot);
    PROT2(pop   = allocMatrix( REALSXP, parlen, bigPopSize ),  numprot);
    PROT2(tmp   = allocMatrix( REALSXP, parlen, basePopSize ), numprot);

    GetRNGstate();

    memcpy( REAL( pop ), 
            REAL( initPar ), 
            parlen * basePopSize * sizeof(double));

    for ( gen = 0; gen < maxit && notConvergedYet; ++gen ) 
    {
        GetNextGen( f, parlen, context, pop, fargs, 
                    fvals, ix, tmp, basePopSize, minimize );

        fmin = REAL( fvals )[ 0 ];
        fmax = REAL( fvals )[ basePopSize - 1 ];

        notConvergedYet = gen < minit || ( ABS( fmax - fmin ) >= tol );

        lagNum = notConvergedYet ? 0 : lagNum + 1;

        if ( !notConvergedYet && lagNum < stopLags && gen < maxit - 1 )
            notConvergedYet = 1;     
    }

    PutRNGstate();

    numAnsCols = notConvergedYet ? basePopSize : 1;

    PROT2( popDim = NEW_INTEGER( 2 ), numprot );
    dimPtr = INTEGER(popDim);

    dimPtr[ 0 ] = parlen;
    dimPtr[ 1 ] = numAnsCols;
  
    PROT2( SET_LENGTH( pop, dimPtr[0] * dimPtr[1] ), numprot );
    setAttrib( pop, R_DimSymbol, popDim );
        
    SET_ELT( ans, 0, pop );
    SET_ELT( ans, 1, SET_LENGTH( fvals, numAnsCols ));  
    SET_ELT( ans, 2, ScalarInteger( notConvergedYet ));

    set_names(ans, names);

    UNPROTECT(numprot);

    return ans;
}


#undef OPTIM_ANS_LEN
#undef OPTIM_ANS_NAMES
#undef BIG_POP_SIZE
#undef C2



#ifdef aoalseijf

static double TestOptimFunc0(SEXP x, SEXP controlPar)
{
    int numprot = 0;
    long i;
    double ans = 0;

    ENSURE_NUMERIC(x, numprot);

    for (i = 0; i < length(x); i++)
        ans += SQR(REAL(x)[i]);

    UNPROT2;

    return ans;

}




SEXP TestOptim00(SEXP initPar, SEXP controlPar, SEXP parlen, SEXP basePopSize, 
                 SEXP stopLags, SEXP minit, SEXP maxit,
                 SEXP minimize, SEXP tol, SEXP relTol, SEXP optimType)
{
  int type = asInteger(optimType);
  SEXP ans;

  switch (type) 
  {
      case 0:

          ans = OptimGa0( testOptimFun, initPar, controlPar,
                          asInteger(parlen), asInteger(basePopSize),
                          asInteger(stopLags), asInteger(minit), asInteger(maxit), 
                          asInteger(minimize), asReal(tol), asInteger(relTol));
          break;
    
      case 1:

          ans = OptimGradient0( testOptimFun, initPar, controlPar,
                                asReal(tol), asInteger(relTol),
                                asInteger(minimize), asInteger(maxit));

          break;


      default:

          ans = R_NilValue;

          break;

  }

  return ans;

}


#endif
