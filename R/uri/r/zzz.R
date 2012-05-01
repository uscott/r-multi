

.First.lib = function(lib, pkg)
{
  try(library(chron))
  try(library(tseries))
  dyn.load("c:/r/library/uri/libs/uri.dll")
}


.Last.lib = function(lib, pkg)
{
  try(dyn.unload("c:/r/library/uri/libs/uri.dll"))
}
