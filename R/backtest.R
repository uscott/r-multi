
max.drawdown = function(x, incremental = FALSE)
{

    .Call("maxDrawdown", x, incremental)
    
}


bsDeltaHedge1 = function(tau, S, K, vols, posn,
r = 0.0, q = r, hedgePeriod = 1,
optionTypes = rep("c", length(K)),
tcosts = 0, 
reltc = TRUE)
    ##
    ##  Function: bsDeltaHedgePL1
    ##
    ##  Purpose:  Compute deltas, P&L from delta-hedging
    ##    using the BS delta, and transaction costs.
    ##
    ##  Params:
    ##
    ##    tau - vector consisting of days from expiration
    ##    at which the underlying price process S is observed.
    ##    Needs to be in descending order for function to
    ##    work properly.  Assumed *not* to include 0
    ##    (i.e. expiration).
    ##    Error if length(tau) != length(S) - 1.
    ##
    ##    S - series of underlying price observations including
    ##    the price at expiration.
    ##    
    ##    K - vector of strike prices corresponding to the
    ##    options composing the position.
    ##
    ##  Returns: List object which includes the values of the
    ##    input arguments and also 'deltas', 'PL', 'dlyPL',
    ##    and 'totalTransCosts' which are respectively the deltas, total
    ##    P&L w/o transaction costs, the P&L by the day,
    ##    and the total transaction costs.
    ##
    ##  Status: .
    ##
    
{

    .Call("bsDeltaHedge1",
          tau, S, K,
          vols, posn, r,
          q, hedgePeriod,
          optionTypes, tcosts, reltc,
          PACKAGE = "uri")

}









bsDeltaHedge2 = function(tau, S, K, vols,
posn, r = 0.0, q = r,
period = 1, option.type = "c", tcosts = 0,
reltc = TRUE, retDly = FALSE)
    ##
    ##  Returns: List with elements
    ##
    ##    pl      - Vector of length = length(tau) where pl[i] is
    ##            P&L from delta-hedging (*not* including transaction costs)
    ##            the given posn from tau[i] days to expiration to
    ##            expiration.
    ##    tc      - tc[i] is total transaction costs associated with
    ##            delta-hedging beginning tau[i] days from expiration.
    ##            The pl[i] + tc[i] is the net delta-hedging P&L.
    ##    dlypl   - If retdlypl then this element is a list of
    ##            length = length(tau) and each dlypl[[i]] is a vector
    ##            whose components are the daily P&L from delta hedging
    ##            beginning tau[i] days from expiration.  Thus
    ##            sum(dlypl[[i]]) == pl[i].
    ##    deltas  - If retdel then this is a list of length = length(tau)
    ##            where deltas[[i]] is a vector giving the deltas by the
    ##            day for delta-hedging beginning at tau[i] days from exp.
    ##    
    ##
    ##  Status:  Tested & OK.
    ##
    
{
    
    .Call("bsDeltaHedge2",
          tau, S, K,
          vols, posn, r,
          q, period, option.type,
          tcosts, reltc, retDly,
          PACKAGE = "uri")

}






bsDeltaHedge3 = function(tau, S, K,
vols, posn,
r = 0, q = r, period = 1,
option.type = "c", tcosts = 0, reltc = TRUE,
retDly = TRUE)
    ##
    ##  Description:  Applies function bsdhpl2 to the *lists*
    ##    tau, S, K, vols, posn, ... where each
    ##    list element tau[[i]], S[[i]], K[[i]], ...
    ##    is of the appropriate type for bsdhpl2.
    ##
    ##  Returns:  List where the ith element is the return
    ##    value from bshpl2 evaluated on tau[[i]], S[[i]], ...
    ##    
    ##  Status: Tested & OK.
    ##

{

    .Call("bsDeltaHedge3",
          tau, S, K, vols, posn, r, q, period,
          option.type, tcosts, reltc, retDly,
          PACKAGE = "uri")
    
}





dvHedge = function(posn, posnStrikes, optionTypes,
hedgePeriod, tau, S, vols, 
volStrikes, r, q, 
undTransCosts, optTransCosts,
isRelUndCosts, isRelOptCosts)
{
    
    .Call("dvHedge",
          posn, posnStrikes, optionTypes,
          hedgePeriod, tau, S, vols, 
          volStrikes, r, q, 
          undTransCosts, 
          optTransCosts, isRelUndCosts, 
          isRelOptCosts,
          PACKAGE = "uri")
    
}




dvHedge2 = function(posn, posnStrikes, optionTypes,
hedgePeriod, tau, S, vols,
volStrikes, r, q,
undTransCosts, optTransCosts,
isRelUndCosts, isRelOptCosts)
{
    
    .Call("dvHedge2",
          posn, posnStrikes, optionTypes,
          hedgePeriod, tau, S, vols, 
          volStrikes, r, q, 
          undTransCosts, 
          optTransCosts, isRelUndCosts, 
          isRelOptCosts,
          PACKAGE = "uri")
    
}





######################################################################
######################################################################
######################################################################
##
##  Following section has functions for delta-hedging of spread
##  options.
##
######################################################################
######################################################################
######################################################################




sprdOptDh1 = function(tau, S1, S2,
K, vol1, vol2,
rho, r = 0, q1 = r,
q2 = r, optionType="c", nGrdPts=100)
{

    .Call("sprdOptDh1",
          tau, S1, S2,
          K, vol1, vol2,
          rho, r, q1,
          q2, optionType, nGrdPts, PACKAGE = "uri")


}



listSprdOptDh1 = function(tau, S1, S2,
K, vol1, vol2,
rho, r = 0, q1 = r,
q2 = r, optionType="c", nGrdPts=100)
{

    .Call("listSprdOptDh1",
          tau, S1, S2,
          K, vol1, vol2,
          rho, r, q1,
          q2, optionType, nGrdPts,
          PACKAGE = "uri")


}

sprdOptDh2 = function(tau,S1,S2,K,
vols1,vols2,volStrikes1,volStrikes2,
rho,r=0,q1=r,q2=r,
optionType="c",nGrdPts=0)
{

    .Call("sprdOptDh2",
          tau,S1,S2,K,
          vols1,vols2,volStrikes1,volStrikes2,
          rho,r,q1,q2,
          optionType,nGrdPts,
          PACKAGE="uri")

}


listSprdOptDh2 = function(tau,S1,S2,K,
vols1,vols2,volStrikes1,volStrikes2,
rho,r=0,q1=r,q2=r,
optionType="c",nGrdPts=0)
{

    .Call("listSprdOptDh2",
          tau,S1,S2,K,
          vols1,vols2,volStrikes1,volStrikes2,
          rho,r,q1,q2,
          optionType,nGrdPts,
          PACKAGE="uri")

}




sprdOptDh3 = function(tau, S1, S2,
K, vol1, vol2,
rho, r = 0, q1 = r,
q2 = r, optionType="c",
nGrdPts=100, maxdt = 7)
{

    .Call("sprdOptDh3",
          tau, S1, S2, K, vol1, vol2,
          rho, r, q1, q2, optionType,
          nGrdPts, maxdt, PACKAGE = "uri")

}



listSprdOptDh3 = function(tau, S1, S2,
K, vol1, vol2,
rho, r = 0, q1 = r,
q2 = r, optionType="c",
nGrdPts = 100, maxdt = 7)
{

    .Call("listSprdOptDh3",
          tau, S1, S2, K, vol1, vol2,
          rho, r, q1, q2, optionType,
          nGrdPts, maxdt, PACKAGE = "uri")

}

sprdOptDh4 = function(tau, S1, S2, K, vol1, vol2, rho,
shortRate = 0, longRate = 0, q1 = shortRate, q2 = shortRate,
optionType = "c", nGrdPts = 100, maxdt = 7)
{
    .Call("sprdOptDh4", tau, S1, S2, K, vol1, vol2, rho,
          shortRate, longRate, q1, q2, optionType,
          nGrdPts, maxdt, PACKAGE = "uri")
}

listSprdOptDh4 = function(tau, S1, S2, K, vol1, vol2, rho,
shortRate = 0, longRate = 0, q1 = shortRate, q2 = shortRate,
optionType = "c", nGrdPts = 100, maxdt = 7)
{
    .Call("listSprdOptDh4", tau, S1, S2, K, vol1, vol2, rho,
          shortRate, longRate, q1, q2, optionType,
          nGrdPts, maxdt, PACKAGE = "uri")
}
