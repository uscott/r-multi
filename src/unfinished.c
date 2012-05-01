


#ifdef U_R_GAY



SEXP optim_ga_list(SEXP fn, SEXP init, SEXP maxit, 
                   SEXP minimize, SEXP tol, SEXP env)

  /*******************************************************************
   *
   *  Description: This is a genetic algorithm optimizer each
   *    argument is a list or vector whose element types correspond to
   *    the arguments of optim_ga1, all of equal lengths with the
   *    exception of env which has length 1.
   *
   *  Status: Incomplete.
   *
   *******************************************************************/

{

  int i, n, numprot = 0;
  SEXP ans, tmp;

  if (!isNewList(fn) || !isNewList(init) || !isNumeric(maxit) ||
      !isInteger(minimize) || !isNumeric(tol) || !isEnvironment(env))
    error ("wrong argument types");
   
  n = length(fn);

  if (n != length(init) || n != length(maxit) || 
      n != length(minimize) || n != length(tol))
    error ("length mismatches in arguments");

  PROTECT(ans = NEW_LIST(n));
  numprot += 1;

  for (i = 0; i < n; i++) {

    
    PROTECT(tmp = optim_ga1(GET_ELT(fn, i), GET_ELT(init, i),
                            GET_ELT(maxit, i), GET_ELT(minimize, i),
                            GET_ELT(tol, i), env));

    numprot += 1;
    SET_ELT(ans, i, tmp);

  }

  SET_NAMES(ans, GET_NAMES(fn));

  UNPROTECT(numprot);

  return(ans);
}






#endif
