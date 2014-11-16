
library( uri )

vsf = data.frame(
    t = c( 1, 1, 1, 2, 2, 2 ),
    m = c( .5, 1, 1.5, .5, 1, 1.5 ),
    v = c( .6, .4, .2, .5, .45, .4 ))

fun1 = function( t, m ) interp.sfc( vsf$t, vsf$m, vsf$v, t, m )

t2   = 10:20/10
m2   = 5:15/10
vsf2 = matrix( NA, length( t2 ), length( m2 ))

for( i in 1:length( t2 ))
    vsf2[ i, ] = fun1( t2[ i ], m2 )

x = as.data.frame(matrix(rnorm(20),10,2))

write.table(x,"tmp.txt",quote=F,row.names=F,col.names=T)



y1 = scan("tmp.txt",nlines=1,what="")
y2 = scan("tmp.txt",skip=1)








test.tmp = function(N = 100, opt = "c",r=0,q1=r,q2=r) {
 = 60
  rho = runif(1, -1, 1)
  vol1 = runif(1)
  vol2 = runif(1)

  
  tau = n:1

  err1 = err2 = rep(NA, N)

  for (i in 1:N) {
    
    eps1 = rnorm(n)
    eps2 = rho * eps1 + sqrt(1-rho^2)*rnorm(n)
    S1 = 2.5
    S2 = runif(1, 0, 5)
    
    S1 = S1 * exp(diffinv((r-q1)/365-0.5*vol1^2/365 + vol1/sqrt(365)*eps1))
    S2 = S2 * exp(diffinv((r-q2)/365-0.5*vol2^2/365 + vol2/sqrt(365)*eps2))
    K  = S2[1] - S1[1]

    price = pearson(tau[1], S1[1], S2[1],
      K, vol1, vol2, rho, r=r,q1=q1,q2=q2,opt=opt)[,"p"]

    hedgePL1 = sprdOptDh1(tau, S1, S2,
      K, vol1, vol2, rho, r=r,q1=q1,q2=q2,opt=opt)

    
    if ("c"==opt)
      payoff = max(last(S2)-last(S1)-K, 0)
    else if ("p"==opt)
      payoff = max(K-(last(S2)-last(S1)), 0)
    else if ("s"==opt)
      payoff = abs(last(S2)-last(S1)-K)
    else
      payoff = NA

    payoff = payoff * exp(-r*tau[1]/365)
    err1[i] = (payoff + hedgePL1 - price)/price


  }

  tmp1 <<- err1
}




f1 = function(mu, sigma, n = 1e4) {

  phi = (1-2*exp(sigma^2)) / (1-exp(sigma^2))
  theta = (phi-0.5)*exp(mu+sigma^2/2)

  x = 1/rgamma(n, shape = phi, rate = theta)

  m1 = mean(x)
  v1 = var(x)
  m2 = exp(mu+sigma^2/2)
  v2 = exp(2*mu+sigma^2)*expm1(sigma^2)

  c(m1=m1,m2=m2,v1=v1,v2=v2,mdiff=abs(m1-m2),vdiff=abs(v1/v2-1))  

}


f2 = function(mu, sigma, n = 1e4) {

  phi = 1/expm1(sigma^2)
  theta = phi * exp(2*mu-sigma^2/2)
  
  x = 1/rgamma(n, shape = phi, rate = theta)

  m1 = mean(x)
  v1 = var(x)
  m2 = exp(mu+sigma^2/2)
  v2 = exp(2*mu+sigma^2)*expm1(sigma^2)
  
  c(m1=m1,m2=m2,v1=v1,v2=v2,mdiff=abs(m1-m2),vdiff=abs(v1/v2-1))
  
}



f3 = function(mu, sigma, n = 1e4) {

  phi = 1/expm1(sigma^2)
  theta = phi * exp(2*mu-sigma^2/2)
  
  x = rgamma(n, shape = phi, rate = theta)

  m1 = mean(x)
  v1 = var(x)
  m2 = exp(-mu+sigma^2/2)
  v2 = exp(-2*mu+sigma^2)*expm1(sigma^2)
  
  c(m1=m1,m2=m2,v1=v1,v2=v2,mdiff=abs(m1-m2),vdiff=abs(v1/v2-1))

}



f4 = function(y, mu, sigma, n = 1e4) {

  
  
}







f=function(S) digital(90,S,1,0.3)
f=function(vol) digital(90,1,1,vol)
f=function(S) digital(90,S,1,0.3,ret="v")
f=function(S) digital(90,S,1,0.3,ret="d")
f=function(S) {function(vol) digital(90,S,1,vol)}


f = function(b,w) 0.5*(w/(w+b) + (50-w)/(100-w-b))


tmp = outer(0:50,0:25,f)


f1 = function(pyou,popp) pyou/(pyou+popp-pyou*popp)
f2 = function(pyou,popp) pyou*(1-popp)/(pyou+popp-pyou*popp)

tmp = outer(seq(0.05,0.95,by=.05),seq(0.05,0.95,by=.05),f)



n = 1e4
P = runif(n)
H = numeric(n)
for (i in 1:n) H[i] = rbinom(1,1,P[i])

P = P[1==H]
h = function(p) {

  tmp = outer(P, p, "<=")
  apply(tmp, 2, sum)/length(P)
  
}


f=function(p) p^3/(p^3+(1-p^2)*(1-p))
