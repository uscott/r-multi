

eval.ng11 = function(x, pars, llonly = FALSE)
    ##
    ##  Description: Returns list including loglikelihood,
    ##    conditional variances and residuals for the given
    ##    time series x and parameters.
    ##
{
    .Call( "ngarch11", x, pars, llonly, PACKAGE = "uri")
}



fit.ng11 = function( x, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, option = NULL)
{
    .Call( "fit_ngarch11", x, init,
          popSize, tol,
          stopLags, minit, maxit, option,
          PACKAGE = "uri")
}


fit.ng11.opcl = function( op, cl, init = NULL, popSize = 5,
    stopLags = 5, minit = 5, maxit = 10,
    tol = 1e-5, option = NULL)
{
    .Call( "fit_ngarch11_opcl", op, cl, init, 
          popSize, tol,
          stopLags, minit, maxit, option,
          PACKAGE = "uri")
}


sim.ng11 = function(n, paths = 1, pars, x0 = NULL, h0 = NULL, z = NULL)
{
    .Call( "sim_ngarch11", n, paths, pars, x0, h0, z, PACKAGE = "uri")
}



