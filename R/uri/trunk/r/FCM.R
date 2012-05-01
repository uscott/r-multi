## Implementation of Fast Convolution Method for pricing derivatives.

## Implementation for standard European call option in Black-Scholes world.
## Note that tau is in days.
Call.FCM.old = function(tau, S, K, vol, r = 0, q = r, nTimeSteps = 1, nGrdPts = 50)
{
    nGrdPts = as.integer(nGrdPts[1])
    nTimeSteps = as.integer(nTimeSteps[1])

    if (nGrdPts < 1 | nTimeSteps < 1)
        stop ("Must have nGrdPts >= 1 and nTimeSteps >= 1")

    N = nTimeSteps # Copying notation in paper.
    M = nGrdPts # Copying notation in paper.
    
    ## Compute time step size:
    dt = tau / 365.0 / N

    mu = r - q - 0.5 * vol^2
    discFac = exp(-r * dt)
    
    ## Choose spacial step size:
    h = 5 * vol * sqrt(tau / 365.0) / M

    ## Define quantity omega (denoted w here):
    w = vol * sqrt(dt) / h
    
    ## Define alpha and beta coefficients for linear interpolations,
    ## denoted here with a and b:
    a = b = rep(NA, 2 * M)
    ## Define vectors z and V for use in holding log-values and call option values:
    z = V = rep(NA, 2 * M + 1)

    ## Initiate z and V:
    for (i in seq(0, 2 * M)) {

        z[i+1] = mu * tau / 365.0 + (-M + i) * h + log(S)
        V[i+1] = max(exp(z[i+1]) - K, 0)
        
    }
    
    for (n in 1:N) {

        ## Get linear interpolation coefficients a and b.
        ## a[i] and b[i] are the slope and intercept respectively
        ## of the piecewise linear interpolation between V[i] and
        ## V[i+1].
        
        for (i in seq(1, 2 * M)) {

            y1 = V[i]; y2 = V[i+1]
            x1 = exp(z[i]); x2 = exp(z[i+1])

            a[i] = (y2 - y1) / (x2 - x1)
            b[i] = (x2 * y1 - x1 * y2) / (x2 - x1)
            
        }

        for (j in seq(0, 2 * M)) {

            V[j+1] = 0
            x = exp(z[j+1] + vol^2 * dt / 2)

            ## Update V
            for (i in seq(0, 2 * M - 1)) {

                V[j+1] = V[j+1] + discFac * a[i+1] * x * (pnorm((i-j+1)/w - vol*sqrt(dt)) - pnorm((i-j)/w - vol*sqrt(dt))) +
                    b[i+1] * (pnorm((i-j+1)/w) - pnorm((i-j)/w))
                
            }

        }

        ## Update z
        z = z - mu * dt
        
    }

    V[M + 1]
}



Call.FCM = function(tau, S, K, vol, nTimeSteps = 1, nGrdPts = 10)
{
    .Call("fcm_EuroCall0", tau, S, K, vol, nGrdPts, nTimeSteps)
}




