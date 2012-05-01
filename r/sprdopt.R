


margrabe = function(tau, S1, S2, vol1, vol2, rho, r = 0,
  q1 = r, q2 = r, option.type = "c", return.type = NULL)
{

  ans = .Call("margrabe",
    tau, S1, S2, vol1, vol2, rho, r, q1, q2, option.type,
    PACKAGE = "uri")

  if (!is.null(return.type))
    ans[, return.type]
  else
    ans
}




margrabe2 = function(tau, S1, S2, K, vol1, vol2, rho,
  r = 0, q1 = r, q2 = r, option.type = "c", return.type = NULL,
  paths = 1e3,
  eps   = NULL)
{

  option.type = tolower(substr(option.type, 1, 1))

  ans = .Call("margrabe2",
    tau, S1, S2, K, vol1, vol2, rho,
    r, q1, q2, option.type, paths, eps,
    PACKAGE = "uri")
    

  if (!is.null(return.type))
    ans[, return.type]
  else
    ans
  
}


pearson = function(tau, S1, S2, K, vol1, vol2, rho,
r = 0, q1 = r, q2 = r, optionType = "c",
calcDeltas = FALSE, nGrdPts = 100)
{

  optionType = tolower(substr(optionType, 1, 1))

  ans = .Call("pearson",
    tau, S1, S2, K, vol1, vol2, rho,
    r, q1, q2, optionType, calcDeltas, nGrdPts,
    PACKAGE = "uri")

  ans
  
}


listPearson = function(tau, S1, S2, K, vol1, vol2, rho,
r = 0, q1 = r, q2 = r, optionType = "c",
calcDeltas = FALSE, nGrdPts = 100) {

  .Call("listPearson",
    tau, S1, S2, K, vol1, vol2, rho,
    r, q1, q2, optionType, calcDeltas, nGrdPts,
    PACKAGE = "uri")

}


  
impCorr = function(tau, S1, S2, K, vol1, vol2, marketPrice,
r = 0, q1 = r, q2 = r, optionType = "c", nGrdPts = 100,
initGuess = NULL, tol = 1e-5, maxit = 100) {

  .Call("impCorr", tau, S1, S2, K,
        vol1, vol2, marketPrice, r,
        q1, q2, optionType, nGrdPts,
        initGuess, tol, maxit,
        PACKAGE="uri")
  

}
