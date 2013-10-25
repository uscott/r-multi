library(uri)

######################################################################

exDts = c( "2013-11-01", "2013-12-01" )

lambda = 0*0.05
a0     = 2.5e-5 + 0*0.30^2 / 365
a1     = 1*0.06
b1     = 1*0.7
gamma  = 1.0

pars = c( a0, a1, b1, gamma )


n   = 2e3
sim = sim.ng11( n, paths = 2, par = pars, z = NULL )
x   = sim$x[,1]
y   = sim$x[,2]

matplot( apply( sim$x, 2, diffinv ), type = "l", lty = 1 )

fit1 = fit.ng11( x, init = c( var( x ), rep( .01, 3 )),
    stopLag = 50, maxit = 5e3, tol = 1e-3,
    option = list( grad = TRUE ))

fit1 = fit.ng11( x,
    init = 0*c( var( x ), rep( .01, 3 ))+1*fit1$par,
    pop = 10,
    stopLag = 1, minit = 100, maxit = 5e3, 
    opt = list( grad = TRUE ))

fit2 = fit.ng11( y, init = c( var( y ), rep( .01, 3 )),
    stopLag = 50, maxit = 5e3, tol = 1e-3,
    option = list( grad = TRUE ))

fit2 = fit.ng11( y,
    init = 0*c( var( y ), rep( .01, 3 ))+1*fit2$par,
    pop = 10,
    stopLag = 1, minit = 100, maxit = 5e3, 
    opt = list( grad = TRUE ))



fit0 = fit.ng11.opcl( x, y,
    ##init = 0*c( var( c( x, y )), rep( .01, 3 )) + 1*fit$par,
    init = as.vector( 0.5 * fit1$par + 0.5 * fit2$par ),
    pop = 10, stopLag = 1, minit = 100, maxit = 5e3, 
    opt = list( grad = TRUE ))


basefit      = tmp1
modInfo      = tmp1
modInfo$freq = 1

if ( FALSE )
{
    modInfo$par[ "lambda", 1] = 0
    modInfo$par[ "a0"    , 1] = 0.30^2 / 365
    modInfo$par[ "a1"    , 1] = 0
    modInfo$par[ "b1"    , 1] = 0
    modInfo$par[ "gamma" , 1] = 0
    modInfo$par[ "h1"    , 1] = 0.30^2 / 365
}


prcs  = ngarch11prices( exDts, modInfo )
vols1 = ngarch11vols( prcs, S0 = 100, strikes = c( 90, 95, 100, 105, 110 ))

modInfo      = tmp1
modInfo$freq = 1

prcs  = ngarch11prices( exDts, modInfo )
vols2 = ngarch11vols( prcs, S0 = 100, strikes = c( 90, 95, 100, 105, 110 ))

bootfit = bootfitNgarch11( x, init = tmp1$par, pop = 10,
    stopLag = 1, minit = 100, maxit = 5e3, fitInit = TRUE, opt = list( grad=FALSE ),
    numBoots = 10 )

bootfit$freq = 1

##debug( ngarch11prices )

bootPrcs = ngarch11prices( exDts, bootfit, paths = 1e5, boot = TRUE )
bootVols = ngarch11vols( bootPrcs, S0 = 100, strikes = c( 90, 95, 100, 105, 110 ))

##tmp2 = wfitNgarch11(x, minit=100, pop=10,maxit = 5e3, tol = 1e-3,
##  fitInit=T, option=list(grad=T))

##tmp2 = wfitNgarch11(x, init=tmp2$par, minit=100, pop=10,maxit = 5e2, tol = 1e-3,
##  fitInit=T, option=list(grad=T))



######################################################################
## Test ar1ngarch11a:
n = 1e3

m0 = 0.0
lambda = 0.0
phi = 0
a0 = 1e-5
b1 = .94
a1 = .04
gamma = .5


par.ng = c(m0 = m0, lambda = lambda, phi = phi, a0 = a0,
  a1 = a1, b1 = b1, gamma = gamma)

tmp1 = sim.ar1ngarch11(n, paths = 1, par = par.ng, z = numeric(0))
##tmp2 = sim.test(n, par = c(par.ng, h0 = 8.75e-5, x.min1 = -8e-3))

x = tmp1$x[,1]
h = tmp1$h

tmp3 = ar1ngarch11(x, par = par.ng)
tmp4 = ar1ngarch11(x, par = par.ng, ll.only = T)

fit = fit.ar1ngarch11(x)


fit = fit.ar1ngarch11(x, method = c("ga", "N", "C", "ga", "gr"),
  maxit = c(50, 5, 20, 20, 1e2), tol = 1e-3)



######################################################################
## Test augmented NGARCH:

n = 1e3

m0 = 0.0
lambda = 0.0

y = cos(2*pi/365 * (1:n))

a00 = 1e-5
a01 = 1e-5
b1 = .94
a1 = .04
gamma = .5

h0 = 1e-4
# simulate path

par = c(m0, lambda, a00, a01, a1, b1, gamma, h0)

tmp = sim.aug.ngarch11(n = n, paths = 1, y = y, par = par)

tmp2 = fit.aug.ngarch11(tmp$x, y, init = par[-1], fit.init = T, method = c("ga", "gr"), maxit = 100)


######################################################################
# Test seasonal NGARCH:


n = 2e3

m0 = 0.0
lambda = 0.0


a00 = 1e-5
a01 = .2
omega = 30
b1 = .94
a1 = .04
gamma = .5

# simulate path
a0 = h = x = rep(NA, n)
yday = as.double(seq(0, n - 1)) %% 365
z = rnorm(n)

a0[1] = a00 * exp(a01 * cos(2*pi/365 * (yday[1] - omega)))
h[1]  = a0[1] / (1-b1-a1*(1+gamma^2))
x[1]  = sqrt(h[1]) * z[1]


for (i in 2:n) {

  a0[i] = a00 * exp(a01 * cos(2*pi/365 * (yday[i] - omega)))
  h[i]  = a0[i] + a1 * h[i-1]*(z[i-1]-gamma)^2 + b1 * h[i-1]
  x[i]  = m0 + lambda * sqrt(h[i]) - 0.5 * h[i] + sqrt(h[i]) * z[i]
  
}


par = c(m0, lambda, a00, a01, omega, a1, b1, gamma)

tmp1 = seasonal.ngarch11(x, yday, par = par, ll.only = T)
tmp2 = seasonal.ngarch11(x, yday, par = par, ll.only = F)


tmp = fit.seasonal.ngarch11(x, yday, maxit = 1e2)

tmp1 = seasonal.ngarch11(x, yday, par = par)
tmp2 = seasonal.ngarch11(x, yday, par = tmp$par)



######################################################################
## Test augmented NGARCH:

n = 750

lambda = 0.03

y = 1 + cos(2*pi/365 * (1:n))

a00 = 1e-4
a01 = 1e-4
b1 = .94
a1 = .04
gamma = .5

h0 = 1e-4
x0 = 0
# simulate path

pars = c(mu = 0.0, lambda = lambda, a00 = a00, a01 = a01,
  a1 = a1, b1 = b1, gamma = gamma)

z = rt(1e5, df = 10)

tmp = sim.aug.ngarch11(n, 1, y, par = pars, x0 = x0, h0 = .001, z = z)
tmp1 = fit.aug.ngarch11(tmp$x, y)


######################################################################
## Test augmented decaying NGARCH:


n = 2e3

m0 = 0.1/365
lambda = 0.1

y = 1 + cos(2*pi/365 * (1:n))
y[-(1:1e3)] = y[1e3]
age = c(rep(1, 1e3), 1:n)[1:n]
dt = rep(1, n)
dt[1e3] = 100

beta = 0.99

a00 = 1e-5
a01 = 1e-4
b1 = .94
a1 = .04
gamma = .5

h0 = 1e-4
# simulate path


par = c(m0, lambda, a00, a01, beta, a1, b1, gamma, h0)

tmp = sim.decay.ngarch11(n = n, paths = 1, dt = 1, y = y,
  age = age, par = par)

tmp2 = fit.decay.ngarch11(tmp$x, dt, y, age, init = par[-1], fit.init = T)
