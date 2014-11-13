#include "uri.h"

void interp_vol_sfc( long    len,
                     double *t,
                     double *m, 
                     double *v, 
                     double  m0, 
                     double  t0, 
                     double *v_out )
{
    const double 
        TOL = 1e-10;    
    double 
        t_sprd_max = -1, m_sprd_max = -1,
        tmp, sprd_lower = -1, sprd_upper = -1,
        t1, t2, m11, m12, m21, m22, v11, v12, v21, v22;    
    int
        found_lower_t = 0, found_upper_t = 0,
        found_lower_m = 0, found_upper_m = 0;
    long 
        i, k1, k2;
    
    /* Get max sprd from t0 */
    for( i = 0; i < len; ++i )
    {        
        if( ( tmp = ABS( t[ i ] - t0 )) > t_sprd_max )
            t_sprd_max = tmp;    

        if( ( tmp = ABS( m[ i ] - m0 )) > m_sprd_max )
            m_sprd_max = tmp;        
    }
    
    for( i = 0, sprd_lower = sprd_upper = t_sprd_max + 1; i < len; ++i )
    {
        tmp = ABS( t[ i ] - t0 );
        
        if( t[ i ] <  t0 && tmp < sprd_lower )
        {
            t1 = t[ i ];
            sprd_lower = tmp;
            found_lower_t = 1;            
        }
    
        if( t[ i ] >= t0 && tmp < sprd_upper )
        {
            t2 = t[ i ];
            sprd_upper = tmp;
            found_upper_t = 1;            
        }        
    }
    
    if( found_lower_t )
    {
        for( i = 0, sprd_lower = sprd_upper = m_sprd_max + 1; i < len; ++i )
        {
            if( ABS( t[ i ] - t1 ) > TOL )
                continue;
            
            tmp = ABS( m[ i ] - m0 );            

            if( m[ i ] <  m0 && tmp < sprd_lower )
            {
                m11 = m[ i ];
                sprd_lower = tmp;
                found_lower_m = 1;                
            }
            
            if( m[ i ] >= m0 && tmp < sprd_upper )
            {
                m12 = m[ i ];
                sprd_upper = tmp;
                found_upper_m = 1
            }            
        }
        
        if( found_lower_m && found_upper_m )
        {
        }
        
    }
    
    if( found_upper_t )
    {
    }
    
    if( !found_lower_t )
    {
    }
    
    if( !found_upper_t )
    {
    }
    
}
