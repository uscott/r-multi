library( uri )

qqq = read.csv( "QQQ.csv" )
len = nrow(qqq)

qqq = qqq[ ( len - 1*252 ):len, ]
x   = qqq$Return

tmp1 = fitNgarch11(x, stopLag=50,maxit = 5e3, tol = 1e-3, fitInit=TRUE, option=list(grad=T))

tmp1 = fitNgarch11(x, init = tmp1$par, pop=10,
  stopLag=1,minit = 100, maxit = 5e3, fitInit=T,opt=list(grad=TRUE))

modInfo = tmp1
modInfo$freq = 1

exDts = c( "2012-08-17", "2012-09-21" )

prcs = ngarch11prices( exDts, modInfo, paths = 20e3 )

strikes = 60:66

S0   = 62.4 + 0*qqq[ nrow(qqq), "Adj.Close" ]
vols = ngarch11vols( prcs, S0=S0,strikes=strikes)
