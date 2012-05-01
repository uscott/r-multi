cumcov = function(x, y, rev = FALSE) {
  .Call("cum_cov", x, y, rev, PACKAGE = "uri")
}
