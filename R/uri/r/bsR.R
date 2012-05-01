bsR = function( tau, S, K, vol, discRate = 0,
  div = discRate, date1 = NULL, date2 = NULL, optionType="c",returnType="p" ) {
# Last modified 6/15/02 5 PM.
# tau is in calendar days
   if (  min(tau,S,K,vol,na.rm=T) < 0 ) retVal = NA
   if ( !is.null(date1) & !is.null(date2) ) {
      tau = date2-date1
      if ( tau < 0 ) retVal = NA
   }
   tau = tau/365
   optionType = substr(optionType,1,1)
   returnType = substr(returnType,1,1)
   d1 = {log(S/K)+(discRate-div+vol^2/2)*tau}/vol/sqrt(tau)
   d2 = d1 - vol*sqrt(tau)
   retVal = NA
   if ( returnType %in% c("p","P") ) {
      if (  optionType %in% c("c","C") ) {
         retVal = exp(-div*tau)*S*pnorm(d1)-exp(-discRate*tau)*K*pnorm(d2)
      } else if ( optionType %in% c("p","P") ) {
          retVal = exp(-discRate*tau)*K*pnorm(-d2) - exp(-div*tau)*S*pnorm(-d1)   
      } else if ( optionType %in% c("s","S") ) {
         retVal = exp(-div*tau)*S*{2*pnorm(d1)-1}-exp(-discRate*tau)*K*{2*pnorm(d2)-1}
      }
   } else if ( returnType %in% c("d","D") ) {
      if (  optionType %in% c("c","C") ) {
         retVal = exp(-div*tau)*pnorm(d1)
      } else if ( optionType %in% c("p","P") ) {
         retVal = exp(-div*tau)*{pnorm(d1)-1}
      } else if ( optionType %in% c("s","S") ) {
         retVal = exp(-div*tau)*{2*pnorm(d1)-1}
      }
   } else if ( returnType %in% c("g","G") ) {
      retVal = dnorm(d1) / {S * vol * sqrt(tau)}
      if ( optionType %in% c("s","S") ) {
         retVal = 2 * retVal
      }
   } else if ( returnType %in% c("v","V") ) {
      retVal =  S * sqrt(tau) * dnorm( d1 )
      if ( optionType %in% c("s","S") ) {
         retVal = 2 * retVal
      }
   }
   retVal
}




bsGamma = function(tau, S, K, vol, r=0, q=r, optionType="c") {

  optionType = tolower(substr(optionType[1],1,1))
  
  tau = tau/365
  d1 = (log(S/K) + (r-q+vol^2/2)*tau)/(vol*sqrt(tau))
  ans = exp(-q*tau)*dnorm(d1)/(S*vol*sqrt(tau))

  if (optionType == "s")
    ans = 2 * ans
  
  ans
  
}
