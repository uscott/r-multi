



library(uri)

n = 1e2
vols = rep(.3, n)
tau = n:1
S = sim.gbm(vols[1], n)
names(S) = paste("t", 1:length(S) - 1, sep = "")
names(tau) = dropLast(names(S))

tmp1 = bsDeltaHedge1(tau, S, S[1], vols, posn = 1, opt = "s", r = c(.05, .1))
tmp2 = bsDeltaHedge2(tau, S, S[1], vols, posn = 1, retDly = F, opt = "s")






test1 = function(numDays = 100, trials = 1e2,
  r = 0.05, period = 1,
  tcosts = .01, reltc = FALSE)
  
{

  
  m    = numDays
  vols = rep(.3, m)
  tau  = m:1
  df   = exp(-r * 1:m / 365)
  
  err1 = err2   = tc = del = rep(NA, trials)

  S   = sim.gbm(vols[1], m, paths = trials, S = 10)
  
  for (tr in 1:trials) {
    
    n    = sample(seq(8), 1)
    opt  = sample(c("c", "p", "s"), n, rep = TRUE)
    posn = sample(c(-(5:1), 1:5),   n, rep = TRUE)
    
    K = round(sort(runif(n, 0.5, 1.5)) * S[1, tr], 1)
    
    f    = rep(NA, n) ## Position premia
    po   = rep(NA, n) ## Position payoffs
    dels = matrix(0, m + 1, n) ## Deltas calculated directly
    
    pfn = function(x, k, opt) { ## Payoff function
      
      if ("s" == opt)
        abs(x - k)
      else if ("c" == opt)
        pmax( x - k, 0)
      else if ("p" == opt)
        pmax(-x + k, 0)
      else
        NA
      
    }
    
    for (j in seq(n)) {
      
      f[j]  = bs(m, S[1, tr], K[j], vols[1], opt = opt[j], r = r)
      po[j] = exp(-r * m / 365) * pfn(S[m + 1, tr], K[j], opt[j])
      
      dels[1:m , j] =
        posn[j] *
          bs(tau, S[1:m, tr], K[j], vols[1], r = r, ret = "d", opt = opt[j])
      
      dels[m + 1, j] =
        posn[j] *
          sign(S[m + 1, tr] - K[j]) * sign(pfn(S[m + 1, tr], K[j], opt[j]))
    }

    f    = sum(posn * f) ## Net premium
    po   = sum(posn * po) ## Net payoff
    dels = apply(dels, 1, sum) ## Net position delta
    
    hpl = bsDeltaHedge1(tau, S[, tr], K, vols,
      posn = posn, r = r, hedgePeriod = period,
      opt = opt, tcosts = tcosts, reltc  = reltc)

    err1[tr]   = (po - f + hpl$PL) / S[1, tr]
    err2[tr]   = (po - f - sum(df * dropLast(dels) * diff(S[, tr]))) / S[1, tr]
    
    del[tr] = sum((hpl$del - dels)^2)

    if (reltc)
      tmp = (sum(df * abs(diff(dels)) * S[-1, tr]) + abs(dels[1]) * S[1, tr]) * tcosts
    else
      tmp = (sum(df * abs(diff(dels))) + abs(dels[1])) * tcosts
    
    tc[tr]  = hpl$totalTransCosts - tmp
    
  }

  hedgeErr1 <<- err1 ## Should have mean approx. = 0
  hedgeErr2 <<- err2 

  del <<- del ## Should be close to 0
  tc  <<- tc ## Should be close to 0

  dels <<- dels
  tmp  <<- tmp
  hpl  <<- hpl
}



 






test2 = function(tauLen = 1e2, r = 0.05, period = 1, tcosts = .01, reltc = FALSE)
{
  ## Function bsdhpl2 returns a list consisting of pl, tc, dlypl, deltas all of length
  ## equal to the length of argument tau.  'pl' and 'tc' are vectors and dlypl,
  ## deltas are lists.
  
  vols = rep(.3, tauLen)
  tau  = tauLen:1
  discFac   = exp(-r * 1:tauLen / 365) # Discount factors
  PL        = rep(NA, tauLen)
  hedgeErr1 = rep(NA, tauLen)
  hedgeErr2 = rep(NA, tauLen)
  totalTransCosts = rep(NA, tauLen)
  
  S   = sim.gbm(vols[1], tauLen, paths = 1, S = 10)

  names(S)   = paste("t", 1:length(S) - 1, sep = "")
  names(tau) = dropLast(names(S))
  
  numStrikes = sample(seq(8), 1)
  opt        = sample(c("c", "p", "s"), numStrikes, rep = TRUE)
  posn       = sample(c(-(5:1), 1:5),   numStrikes, rep = TRUE)

  K = S[1:tauLen] %o% sort(runif(numStrikes, 0.5, 1.5))
  K = round(K, 1)

  
  prem = matrix(NA, tauLen, numStrikes) ## Position premia
  poff = matrix(NA, tauLen, numStrikes) ## Position payoffs
  deltas = list() ## Deltas calculated directly

  for (k in 1:tauLen)
    deltas[[k]] = matrix(NA, tauLen - k + 2, numStrikes)
  
  pfn = function(x, k, opt) { ## Payoff function

    if ("s" == opt)
      abs(x - k)
    else if ("c" == opt)
      pmax( x - k, 0)
    else if ("p" == opt)
      pmax(-x + k, 0)
    else
      NA
    
  }


  
  for (j in 1:numStrikes) { ## Loop over strikes
    
    prem[, j] = bs(tau, S[1:tauLen], K[, j], vols[1], opt = opt[j], r = r)
    
    poff[, j] = exp(-r * tau / 365) * pfn(S[tauLen + 1], K[, j], opt[j])

    for (k in 1:tauLen) { ## Loop over times to expiration

      M = tauLen - k + 2
      
      deltas[[k]][-M, j] = posn[j] *
        bs(tau[k:tauLen], S[k:tauLen], K[k, j], vols[1], r = r, ret = "d", opt = opt[j])

      deltas[[k]][M, j] = posn[j] *
        sign(S[tauLen + 1] - K[k, j]) * sign(pfn(S[tauLen + 1], K[k, j], opt[j]))

    }
    
  }

  prem = prem %*% posn # Net premium, now prem is a vector of length tauLen
  poff = poff %*% posn # Net payoff, now vector of length tauLen
  
  for (k in 1:tauLen) {
    
    deltas[[k]] = apply(deltas[[k]], 1, sum) ## Net position delta
    PL[k] = sum(-discFac[seq(tauLen - k + 1)] * dropLast(deltas[[k]]) * diff(S[seq(k, tauLen + 1)]))
    
  }
    
  hpl = bsDeltaHedge2(tau, S, K, vols,
    r = r, posn = posn, per = period,
    opt = opt, tcosts = tcosts, reltc  = reltc, retDly = TRUE)

  deltas.SSE = rep(0, tauLen)
  totalTransCosts = rep(0, tauLen)
  
  for (j in 1:tauLen) {

    hedgeErr1[j]  = (poff[j] - prem[j] + hpl$PL[j]) / S[j]
    hedgeErr2[j]  = (poff[j] - prem[j] + PL[j]) / S[j]
    
    deltas.SSE[j] = sum((hpl$deltas[[j]] - deltas[[j]])^2)

    M = tauLen - j + 1
    
    if (reltc)
      totalTransCosts[j] = (sum(discFac[1:M] * abs(diff(deltas[[j]])) * S[j:tauLen + 1]) +
               abs(deltas[[j]][1]) * S[j]) * tcosts
    else
      totalTransCosts[j] = (sum(discFac[1:M] * abs(diff(deltas[[j]]))) + abs(deltas[[j]][1])) * tcosts
    
    
  }

  
  hedgeErr1 <<- hedgeErr1
  hedgeErr2 <<- hedgeErr2
  totalTransCosts  <<- totalTransCosts

  deltas <<- deltas
  hpl  <<- hpl
  
}





test3 = function(tc = .01)
{

  L = 20
  
  vols = tau = S = K = opt = posn = r = list()
  ans1 = list()
  
  for (j in 1:L) {

    tauLen = sample(50:100, 1)

    vols[[j]] = rep(runif(1), tauLen)
    S[[j]]    = sim.gbm(vols[[j]][1], tauLen)
    tau[[j]]  = tauLen:1
    r[[j]]    = rep(runif(1, 0, .1), tauLen)
    
    numStrikes = sample(seq(5), 1)

    opt[[j]] = sample(c("c", "p", "s"), numStrikes, rep = TRUE)
    posn[[j]] = sample(c(-(5:1), 1:5),  numStrikes, rep = TRUE)
    
    K[[j]] = S[[j]][1:tauLen] %o% sort(runif(numStrikes, 0.5, 1.5))
    K[[j]] = round(K[[j]], 1)

    ans1[[j]] = bsDeltaHedge2(tau[[j]], S[[j]], K[[j]], vols[[j]],
          posn = posn[[j]], r = r[[j]],
          period = 1, opt = opt[[j]], tcosts = tc, reltc = FALSE)

  }
  
  ans2 <<- bsDeltaHedge3(tau = tau, S = S, K = K, vols = vols,
                         posn = posn, r = r, opt = opt, tcosts = tc,
                         reltc = FALSE)
  ans1 <<- ans1
}



















K = list(as.matrix(S))
S = list(S)
tau = list(tau)
posn = list(1)
vols = list(vols)
r = list(rep(0, n))
div = list(rep(0, n))
period = list(1)
option.type = list("c")
tcosts = list(rep(0, n + 1))
reltc = list(TRUE)
retdel = list(TRUE)
retdlypl = list(TRUE)


test3 = function() {

  L = 10
  n = 1e2
  v = rep(.3, n)

  tau = S = vols = r = div = K = posn =
    period = option.type = tcosts = reltc =
      retdel = retdlypl = list()
  
  for (j in 1:L) {
    
    vols[j] = list(rep(.3, n))
    tau[j]  = list(n:1)
    S[j]    = list(sim.gbm(vols[[j]][1], n)) 
    K[j]    = list(as.matrix(S[[j]]))
    posn[j] = list(1)
    r[j]    = list(rep(0, n))
    period[j] = list(1)
    option.type[j] = list("c")
    tcosts[j] = list(rep(0, n + 1))
    reltc[j] = list(TRUE)
    retdel[j] = list(TRUE)

  }

  retdlypl = TRUE
  
  div = r

  tmp <<- bsdhpl3(tau, S, K, vols, posn, r, div, period,
                  option.type, tcosts, reltc, retdel, retdlypl)
}





























test3 = function(trials = 1e2, r = 0.05, period = 1, tcosts = 0, reltc = TRUE)
{
  m   = 1e2
  vol = .3
  tau = m:1
  
  e   = rep(NA, trials)
  tc  = rep(NA, trials)
  del = rep(NA, trials)

  
  for (tr in 1:trials) {
    
    n = sample(seq(2, 8), 1)

    opt = sample(c("c", "p", "s"), n, rep = TRUE)
    
    posn = sample(c(-(5:1), 1:5), n, rep = TRUE)
    
    S  = sim.gbm(vol, m, S = 10);
    
    K = round(sort(runif(n, 0.5, 1.5)) * S[1], 1)
    
    f  = rep(NA, n)
    po = rep(NA, n)
    dels = matrix(0, m + 1, n)
    
    pfn = function(x, k, opt) {
      
      if ("s" == opt)
        abs(x - k)
      else if ("c" == opt)
        pmax( x - k, 0)
      else if ("p" == opt)
        pmax(-x + k, 0)
      else
        NA
      
    }
    
    for (j in seq(n)) {
      
      f[j]  = bs(m, S[1], K[j], vol, opt = opt[j], r = r)
      po[j] = exp(-r * m / 365) * pfn(last(S), K[j], opt[j])
      
      dels[1:m , j] = posn[j] * bs(m:1, S[1:m], K[j], vol,
            r = r, ret = "d", opt = opt[j])
      ##dels[m + 1, j] = posn[j] * sign(po[j])
    }

    f    = sum(posn * f)
    po   = sum(posn * po)
    dels = apply(dels, 1, sum)
    
    hpl = dhedge.pl3(tau, S, K, vol,
      r      = r,
      posn   = posn,
      per    = period,
      opt    = opt,
      tcosts = tcosts,
      reltc  = reltc)

    e[tr]   = (po - f + hpl$pl) / S[1]
    del[tr] = sum((hpl$del[1:m] - dels[1:m])^2)
    tc[tr]  = hpl$tc
    
  }

  e   <<- e
  del <<- del
  tc  <<- tc
}


















test1 = function(n, S, vol, period = 1, opt = "c", tcosts = 0, reltc = TRUE) {
## Test dhedge.pl1

  f = bs(n, S, S, vol, opt = opt)
  S = S * exp(diffinv(-0.5 * vol^2 / 365 + vol / sqrt(365) * rnorm(n)))


  pfn2 = function(x, opt) {
    if ("s" == opt)
      abs(last(x) - x[1])
    else if ("c" == opt)
      max(last(x) - x[1], 0)
    else if ("p" == opt)
      max(x[1] - last(x), 0)
    else
      NA
  }

  
  hpl = dhedge.pl1(n:1, S, S[1], vol, per = period,
    opt = opt, tcosts = tcosts, reltc = reltc)$pl

  (pfn2(S, opt) - f + hpl) / f
}





test2 = function(n, S, vol, period = 1, opt = "c", tcosts = 0, reltc = TRUE) {

  S = S * exp(diffinv(-0.5 * vol^2 / 365 + vol / sqrt(365) * rnorm(n)))
  f = bs(n:1, S[1:n], S[1:n], vol, opt = opt)


  pfn = function(x, opt) {

    n = length(x) - 1
  
    if ("s" == opt)
      abs(x[n + 1] - x[1:n])
    else if ("c" == opt)
      pmax( x[n + 1] - x[1:n], 0)
    else if ("p" == opt)
      pmax(-x[n + 1] + x[1:n], 0)
    else
      NA
  }

  

  hpl = dhedge.pl2(n:1, S, S, vol, per = period, opt = opt,
    tcosts = tcosts, reltc = reltc)

  (pfn(S, opt) - f + hpl) / f
}




for (i in seq(1e3)) {

  if (1 == i) {
    e = 0
    Sys.sleep(5)
  }
    
  h = test2(1e2, 100, .5)

  e = ((i - 1) * e + h) / i
  
  matplot(cbind(e, h), type = "l", ylim = c(-.1, .1)); abline(0, 0, col = 4)

  if (1 == i)
    Sys.sleep(2)
  else
    Sys.sleep(0.01)
  
}












vol = .4; m = 500; tau = m:1; S = sim.gbm(vol, m)
K = cbind(S, 1.05 * S, 1.1 * S); posn = c(1.0, -2.0, 1.0)
opts = c("c", "c", "c")

tmp4 = dhedge.pl4(tau, S, as.matrix(K), vol, posn = posn,
  opt = opts, tcost = rep(0, m + 1), retdel = TRUE)



i = 25
tmp3 = dhedge.pl3(tau[i:m], S[i:length(S)], K[i, 1], vol,
  posn = posn[1], opt = opts[1], tcost = rep(0, m - i + 2))







test4 = function(trials) {

  vol = .3
  m = 50
  
  e  = rep(NA, trials)
  de = rep(NA, trials)

  
  for (tr in 1:trials) {
    
    n = sample(seq(2, 8), 1)

    opt = sample(c("c", "p", "s"), n, rep = TRUE)
    
    posn = sample(c(seq(-5, -1), seq(1, 5)), n, rep = TRUE)
    
    S  = sim.gbm(vol, m, S = 10);
    K = round(S %o% sort(runif(n, 0.5, 1.5)), 2)
    f  = matrix(NA, m, n)
    po = matrix(NA, m, n)
    dels = list()

    
    pfn = function(x, k, opt) {

      n = length(x) - 1
      
      if ("s" == opt)
        abs(x - k)
      else if ("c" == opt)
        pmax( x - k, 0)
      else if ("p" == opt)
        pmax(-x + k, 0)
      else
        NA
      
    }

    
    for (i in seq(m))
      dels[[i]] = 0
    
    for (j in seq(n)) {
      
      f[, j]  = bs(m:1, S[1:m], K[1:m, j], vol, opt = opt[j])
      po[, j] = pfn(last(S), K[1:m, j], opt[j])

      for (i in seq(m))
        dels[[i]] = dels[[i]] +
          posn[j] * bs(i:m, S[i:m], K[i, j], vol, ret = "d", opt = opt[j])
      
    }
    
    f  = drop(f %*%  posn)
    po = drop(po %*% posn)

    
    hpl = dhedge.pl4(m:1, S, K, vol, posn = posn, per = 1,
      opt = opt, tcost = rep(0, m + 1), reltc = TRUE, retdel = TRUE)

    e[tr] = (po - f + hpl$pl) / S[1]
    de[tr] = sum((hpl$del[1:m] - dels[1:m])^2)

    
  }

  cbind(e, de)
  
}

  





######################################################################
######################################################################
######################################################################
##
##
##
######################################################################
######################################################################
######################################################################













n = 1e2
vol = .3
S = 4 * exp(diffinv(-0.5 * vol^2/365 + vol/sqrt(365) * rnorm(n)))
f = bs(n, S[1], S[1], vol, opt = opt)




tcosts = seq(0, .01, by = 1e-3)

ftmp = function(tc) dhedge.pl1(n:1, S, S[1], vol,
  tcosts = tc, per = per, reltc = FALSE, opt = opt)$pl


h = sapply(tcosts, ftmp)

hpl = (pfn(S, opt) - f + h) / f
plot(hpl)
  


efn = function(x, n = 1e2, vol = .3, per = 1, opt = "s") {
  
  S = 4 * exp(diffinv(-0.5 * vol^2/365 + vol/sqrt(365) * rnorm(n)))
  
  ftmp = function(tc) dhedge.pl1(n:1, S, S[1], vol,
    tcosts = tc, per = per, reltc = TRUE, opt = opt)$pl

  sapply(tcosts, ftmp)
  
}


n = 1e2
vol = .5
trials = 5e2
f = bs(n, 1, 1, vol, opt = opt)
per = 1

efn = function(x) {
  
  S = exp(diffinv(-0.5 * vol^2/365 + vol/sqrt(365) * rnorm(n)))

  tmp = dhedge.pl1(n:1, S, S[1], vol, opt = opt, per = per, tc = .01, rel = F)

  ans = list()
  ans$pl = pfn(S, opt) - f + tmp$pl
  ans$del = tmp$del
  
  ans
}


ctr = 0

repeat {

  ctr = ctr + 1
  
  S = exp(diffinv(-0.5 * vol^2/365 + vol/sqrt(365) * rnorm(n)))
  e = dhedge.pl1(n:1, S, S[1], vol, opt = opt, per = per, tc = .01, rel = T)
  
  if (is.na(e$pl) | ctr > 1e5 | is.na(ctr))
    break
}



n = 2e2
vol = .3
mu = 0

efn = function(x) {

  S = exp(diffinv(mu / 365 -0.5 * vol^2/365 + vol/sqrt(365) * rnorm(n)))
  
  f = bs(n:1, S[1:n], S[1:n], vol)
  
  (pmax(last(S) - S[1:n], 0) - f + dhedge.pl2(n:1, S, S, vol, per = 1)) / f
  
}

trials = 5e2
e = sapply(seq(trials), efn)







n = 2e2
vol = .5

e1 = rep(NA, n)
S  = exp(diffinv(-0.5 * vol^2 / 365 + vol / 365 * rnorm(n)))
f  = bs(n:1, S[1:n], S[1:n], vol)

for (i in 1:n)
  e1[i] = dhedge.pl((n - i + 1):1, S[i:n], S[i], vol)$pl



######################################################################
## Test vega-hedging stuff


library(uri)

n = 1e2
vols = matrix(0.4, n, 11)
volStrikes = matrix(seq(0.75, 1.25, by = 0.05), n, 11, byrow = TRUE)
tau = n:1
S = sim.gbm(vols[1], n)
names(S) = paste("t", 1:length(S) - 1, sep = "")
names(tau) = dropLast(names(S))
posn = c(1.0, -1.0)
posnStrikes = c(1.05, 1.15)
optionTypes = c("c", "c")
r = q = 0.1

tmp = dvHedge(posn, posnStrikes, optionTypes,
  1, tau, S, vols, volStrikes, r, q, 0.01, 0.05, 0, 1)


posnStrikes = matrix(posnStrikes, n, length(posnStrikes), byrow = TRUE)

tmp = dvHedge2(posn, posnStrikes, optionTypes, 1,
  tau, S, vols, volStrikes, r, q, 0.01, 0.05, 0, 0)


test.dvh = function(hedgePeriod = 1, undTC = .01, optTC = .01, isRelUndTC=FALSE,isRelOptTC=FALSE)
{

  tauLen = 1e2
  tau = tauLen:1
  dt = c(-diff(tau), last(tau))
  vol = arima.sim(model = list(ar = 0.9), n = tauLen)
  vol = exp(vol)
  vol = vol / mean(vol) * 0.4 / sqrt(365)
  S = vol * rnorm(tauLen) -0.5 * vol^2
  S = exp(diffinv(S))

  vol = rev(cum.mean(rev(sqrt(365) * vol)))

  q = r = 0.05

  strikeMatrix = matrix(seq(0.5, 2, by = 0.1), tauLen, 16, byrow = TRUE)
  volMatrix = matrix(vol, tauLen, 16)
  numVolCols = 16

  numStrikes = 3
  posn = sample(c(-(5:1),1:5), numStrikes, rep = TRUE)
  K    = sample(seq(0.8,1.2,by=.05), numStrikes)
  optionTypes = sample(c("c", "p", "s"), numStrikes, rep=TRUE)
  
  ## Get ATM strikes
  atmStrike = rep(NA, tauLen + 1)
  
  for (i in 1:tauLen) {

    j = which.min(abs(S[i] - strikeMatrix[i, ]))

    atmStrike[i] = strikeMatrix[i, j]
    
  }
  
  ## Get deltas & vegas of base posn
  deltas = vegas = numeric(tauLen + 1)

  for (j in 1:numStrikes) {
    
    tmp = bs(c(tau,0), S, K[j], c(vol,0), opt=optionTypes[j],ret=NULL,r=r)
    deltas = deltas + posn[j]*tmp[,"d"]
    vegas  = vegas + posn[j]*tmp[,"v"]
    
  }  
  
  ## Get hedges
  undPosn = strdPosn = strdStrike = strdPrice = prevStrdPrice = numeric(tauLen+1)

  dayCtr = 0
  tmp = bs(tau[1],S[1],atmStrike[1],vol[1],opt="s",ret=NULL,r=r)
  strdStrike[1] = atmStrike[1]
  strdPosn[1] = -vegas[1]/tmp["v"]
  undPosn[1] = -deltas[1] - strdPosn[1] * tmp["d"]
  strdPrice[1] = tmp["p"]
  prevStrdPrice[1]=NA
  
  for (i in 2:tauLen) {

    dayCtr = dayCtr + tau[i-1]-tau[i]
    prevStrdPrice[i]=bs(tau[i],S[i],strdStrike[i-1],vol[i],opt="s",ret="p",r=r)
    
    if (dayCtr >= hedgePeriod) {

      dayCtr = 0

      tmp = bs(tau[i],S[i],atmStrike[i],vol[i],opt="s",ret=NULL,r=r)
      strdStrike[i] = atmStrike[i]
      strdPosn[i] = -vegas[i]/tmp["v"]
      undPosn[i] = -deltas[i]-strdPosn[i]*tmp["d"]
      strdPrice[i]=tmp["p"]

    }
    else {

      strdStrike[i] = strdStrike[i-1]
      undPosn[i] = undPosn[i-1]
      strdPosn[i] = strdPosn[i-1]
      strdPrice[i]=bs(tau[i],S[i],strdStrike[i],vol[i],opt="s",ret="p",r=r)
      
    }
  }

  bs(0, last(S), strdStrike[tauLen],0.1,opt="s",ret=NULL,r=r)
  undPosn[tauLen+1] = -deltas[tauLen+1] - strdPosn[tauLen]*tmp["d"]
  prevStrdPrice[tauLen+1] = abs(last(S) - strdStrike[tauLen])
  
  ## Compute hedge P&L & TC
  discFac = exp(-r * (max(tau)-tau+1)/365)
  undPL = sum(discFac * undPosn[1:tauLen] * ((S[-1]-dropLast(S)) + dropLast(S)*expm1((q-r)*dt)))
  strdPL = sum(discFac * strdPosn[1:tauLen] * (prevStrdPrice[-1] - exp(r*dt/365)*dropLast(strdPrice)))

  if (isRelUndTC) {

    totUndTC = undTC * (sum(discFac*abs(diff(undPosn)) * S[-1]) + abs(undPosn[1])*S[1])
    
  } else {
    totUndTC = undTC * (sum(discFac*abs(diff(undPosn))) + abs(undPosn[1]))
  }

  totOptTC = abs(strdPosn[1]) * optTC
  if (isRelOptTC)
    totOptTC = totOptTC * strdPrice[1]

  for (i in 2:tauLen) {

    if (strdStrike[i-1] == strdStrike[i] & !isRelOptTC)
      optTCtoday = optTC * abs(strdPosn[i] - strdPosn[i-1])
    else if (strdStrike[i-1] == strdStrike[i] & isRelOptTC)
      optTCtoday = optTC * abs(strdPosn[i] - strdPosn[i-1]) * strdPrice[i]
    else if (!isRelOptTC)
      optTCtoday = optTC * max(abs(strdPosn[(i-1):i]))
    else
      optTCtoday = optTC * max(abs(strdPosn[(i-1):i])) * max(strdPrice[i], prevStrdPrice[i])

    totOptTC = totOptTC + discFac[i] * optTCtoday
    
  }
  

  
  tmp <<- dvHedge(posn, K, optionTypes, hedgePeriod,
                  tau, S, volMatrix, strikeMatrix, r, q,
                  undTC, optTC, isRelUndTC, isRelOptTC)
  
  undPL <<- undPL
  strdPL <<- strdPL
  undPosn <<- undPosn
  strdPosn <<- strdPosn
  strdPrice <<- strdPrice
  prevStrdPrice <<- prevStrdPrice
  deltas <<- deltas
  vegas <<- vegas
  strdStrike <<- strdStrike
  undTC <<- totUndTC
  optTC <<- totOptTC
}



