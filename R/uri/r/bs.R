

bs = function(tau, S, K, vol, r = 0, q = r,
  optionType = "c", returnType = "p") {
  
  ans = .Call("bs", tau, S, K, vol, r,
    q, optionType, returnType, PACKAGE = "uri")

  drop(ans)
}




bs2 = function(tau, S, K, vol, posn = 1, r = 0, q = r,
  option.type = "c", return.type = "p")
{
  n = length(K)
  
  if (1 == length(option.type))
    option.type = rep(option.type, n)
  if (1 == length(posn))
    posn = rep(posn, n)
  
  stopifnot (n == length(option.type))
  stopifnot (1 == max(length(S), length(tau)))

  bsval =
    bs(tau, S, K, vol, r = r, q = q, opt = option.type, ret = return.type)

  t(posn) %*% bsval
  
}


digital = function(tau, S, K, vol, r = 0, q = r, returnType = "p") {

  .Call("digital", tau, S, K, vol, r, q, returnType, PACKAGE = "uri")
  
}

