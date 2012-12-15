library( uri )

qqq = read.csv( "QQQ.csv" )
len = nrow(qqq)

qqq = qqq[ ( len - 1*252 ):len, ]
x   = qqq$Return

plot( exp( diffinv( x )))

tmp1 = fitNgarch11( x, stopLag = 50, maxit = 5e3, tol = 1e-3,
    fitInit = TRUE, option = list( grad = TRUE ))

tmp1 = fitNgarch11( x, init = tmp1$par, pop = 10,
    stopLag = 1, minit = 100, maxit = 5e3, fitInit = TRUE,
    opt = list( grad  = TRUE ))

tmp2 = bootfitNgarch11( x, init = tmp1$par, pop = 10,
    stopLag = 1, minit = 500, maxit = 5e3, fitInit = TRUE,
    opt = list( grad = FALSE ),
    numBoots = 19 )

indices = which( tmp2$bootPar["a1",] < 0 | tmp2$bootPar["b1",] < 0 )

tmp2$bootPar   = tmp2$bootPar[,-indices]
tmp2$bootVar   = tmp2$bootVar[,-indices]
tmp2$bootRes   = tmp2$bootRes[,-indices]
tmp2$bootPaths = tmp2$bootPaths[,-indices]

tmp1$freq = 1
tmp2$freq = 1

modInfo      = tmp1
modInfo$freq = 1

exDts = c( "2013-01-18", "2013-02-15", "2013-03-15")

prcs  = ngarch11prices( exDts, modInfo, paths = 20e3 )
prcs1 = prcs

bootPrcs =  ngarch11prices( exDts, tmp2, paths = 50e3 )
bootPrcs1 = bootPrcs

divd    = 0*.24
strikes = 64:66

S0        = 64.805 - divd + 0*qqq[ nrow(qqq), "Adj.Close" ]
vols1     = ngarch11vols( prcs1, S0 = S0, strikes = strikes)
bootVols1 = ngarch11vols( bootPrcs1, S0 = S0, strikes = strikes )

atmVol = .145
a0     = modInfo$par["a0",1]
h1     = modInfo$par["h1",1]
a0     = a0 * atmVol^2 / vols1$expiration.date1["69","vols"]^2
h1     = h1 * atmVol^2 / vols1$expiration.date1["69","vols"]^2
modInfo$par["a0",1] = a0
modInfo$par["h1",1] = h1

a0     = tmp2$bootPar["a0",]
h1     = tmp2$bootPar["h1",]
a0     = a0 * atmVol^2 / bootVols1$expiration.date1["69","vols"]^2
h1     = h1 * atmVol^2 / bootVols1$expiration.date1["69","vols"]^2
tmp2$bootPar["a0",] = a0
tmp2$bootPar["h1",] = h1


prcs2 = ngarch11prices( exDts, modInfo, paths = 20e3 )
vols2 = ngarch11vols( prcs2, S0 = S0,strikes = strikes)

bootPrcs2 =  ngarch11prices( exDts, tmp2, paths = 50e3 )
bootVols2 =  ngarch11vols( bootPrcs2, S0 = S0, strikes = strikes )
