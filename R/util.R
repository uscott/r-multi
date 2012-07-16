

fillNAfwd = function(x)
{

    .Call("fillNAfwd", x, PACKAGE = "uri")

}

fillNAbkwd = function(x)
{

    .Call("fillNAbkwd", x, PACKAGE = "uri")
    
}

summary2 = function(x, maxCols = 100)
{

    if (is.matrix(x) | is.data.frame(x)) {

        if (ncol(x) > maxCols) {

            x = x[, 1:maxCols]

            warning("Only displaying some columns")
        }
        
    }

    summary(x)

}



segment = function(x, b, e = 0) {

    isMatrixLike = is.matrix(x) | is.data.frame(x)
    
    if (missing(e))
        e = ifelse(isMatrixLike, nrow(x), length(x))
    
    if (isMatrixLike)
        x[b:e, , drop = FALSE]
    else
        x[b:e]

}



lo = function(n = 10) {
    
    z = sapply(ls(".GlobalEnv"), function(x) object.size(get(x)))
    
    ans = as.matrix(rev(sort(z))[1:n])
    colnames(ans) = "bytes"
    ans[which.not.na(ans[,"bytes"]), , drop = FALSE]
    
}



permTest = function(x, y, f = mean, R = 2e3 - 1, loop = FALSE)
    ##
    ##  Description: Permutation test for difference in functions.
    ##               Gives p-value of f(x) - f(y) under the hypothesis
    ##               that x, y are drawn from same parent distribution.
    ##               f should be scalar valued.
    ##
{
    s0 = f(x) - f(y)

    stopifnot(1 == length(s0))
    
    nx = length(x)
    ny = length(y)
    x  = c(x, y)
    s  = rep(NA, R)

    if (loop) {
        for (b in seq(R)) {

            x    = sample(x)
            s[b] = f(x[1:nx]) - f(x[1:ny + nx])
            
        }
    } else {

        g = function(n) sample(x)
        x = sapply(seq(R), g)

        stopifnot(nrow(x) == nx + ny)
        
        s = apply(x[1:nx, ], 2, f) - apply(x[1:ny + nx, ], 2, f)
    }

    (1 + length(which(s >= s0))) / (1 + R)
    
}





######################################################################


newey.west.se = nwse = function(lin.model, q)
    ##
    ##  Description: Returns the vector of Newey-West
    ##    heteroskedasticity & autocorrelation-consistent
    ##    standard errors corresponding to the lags
    ##    given in argument q.
    ##

{

    n = nrow(lin.model$x)
    x = lin.model$x
    y = lin.model$y
    u = lin.model$res
    
    Qi  = solve(crossprod(x)/n)

    w = x * u

    ans = matrix(NA, nrow(Qi), length(q))
    rownames(ans) = colnames(x)
    colnames(ans) = q

    if (length(q)) {
        
        for (k in 1:length(q)) {
            
            S = 1/n * crossprod(w)

            if (q[k]) {
                
                for (i in 1:q[k]) {
                    
                    gamma = 1/n * t(w[seq(n - q[k]), ]) %*% w[seq(q[k] + 1, n), ]

                    S = S + (1 - i / (q[k] + 1)) * (gamma + t(gamma))
                    
                }
                
            }
            
            ans[, k] = sqrt(diag(1/n * Qi %*% S %*% Qi))
            
        }
        
    }

    ans
    
}




## White's heteroskedasticity-consistent standard error
wse = function(model, method = 1)
{
    if (class(model) == "numeric")
        model = lm(as.vector(model) ~ 1, x = TRUE, y = TRUE)
    stopifnot (class(model) == "lm")
    method = method[1]
    
    y = model$y
    X = model$x
    n = nrow(X)
    K = ncol(X)

    Sxx  = crossprod(X)/n
    Sxxi = solve(Sxx)
    ehat = as.double(model$res)
    
    if (method %in% 1:2) {

        denom = ifelse(1 == method, n, n - K)
        
        S = crossprod(ehat * X)/denom

    } else if (method %in% 3:4) {

        g = function(z)
            t(z) %*% Sxxi %*% z

        p = apply(X, 1, g)

        d = ifelse(3 == method, 1, 2)
        S = crossprod(ehat / (1-p)^d * X)/n
        
    }

    sqrt(diag(Sxxi %*% S %*% Sxxi) / n)
    
}



wtstat = function(model, method = 1)
{
    if (class(model)=="numeric")
        model = lm(as.vector(model) ~ 1, x = TRUE, y = TRUE)
    model$coef / wse(model, method = method)
}




lm.se.boot = function(model, R = 2e3)
{

    n    = nrow(model$x)
    x    = model$x
    yhat = model$fit
    u    = model$res
    p    = length(model$coef)
    
    coef.boot = matrix(NA, R, p)

    for (b in seq(R)) {

        coef.boot[b, ] = lm(yhat + sample(u, rep = TRUE) ~ 0 + x)$coef
        
    }

    apply(coef.boot, 2, sd)
}


######################################################################
######################################################################
######################################################################
##
##  This section contains some functions for
##  (1) reading .csv files downloaded from Yahoo!,
##  (2) reading .txt files downloaded from the Federal Reserve.
##
######################################################################
######################################################################
######################################################################



read.yahoo = function(file)
{

    ans = read.csv(file, as.is = TRUE, comment.char = "")
    
    names(ans) = tolower(names(ans))
    ans$date   = as.char(strptime(ans$date, format = "%d-%b-%y"))

    ix  = sort(posixToChron(ans$date), index = TRUE)$ix
    ans = ans[ix, ]
    
    row.names(ans) = ans$date
    ans$date       = NULL

    ans$lclose  = log(ans$close)
    ans$lreturn = c(NA, diff(ans$lclose))
    ans$t       = difftime.uri(row.names(ans), "2000-01-01")
    ans$dt      = ans$t - lag2(ans$t, 1)
    
    
    ans
}




read.fed.data = function(file, cts = TRUE)
{

    ans = read.table(file, skip = 30, as.is = TRUE, comment.char = "")
    
    names(ans) = c("date", "rate")
    
    ans$date = as.char(strptime(ans$date, format = "%m/%d/%Y"))
    ans$rate = 0.01 * ans$rate ## Convert to straight decimals

    if (cts)
        ans$rate = log((1 + ans$rate / 365)^365)
    
    row.names(ans)  = ans$date
    names(ans$rate) = ans$date
    ans$date = NULL

    ans
}



######################################################################


aic = function(ll, npar)
{
    -2 * ll + 2 * npar
}

bic = function(ll, npar, N)
{
    -2 * ll + npar * log(N)
}



######################################################################



exDatesEquities = function(mth, yr, exp.hr = 16, exp.min = 30)
{
    ##
    ##  Description: Finds the 3rd friday of the month,
    ##    which is essentially the expiration date for
    ##    equity options.
    ##    Note that the values of exp.hr and exp.min may
    ##    need to be adjusted for the time zone.
    ##

    require(chron)
    
    mth = sort(unique(mth))
    yr  = sort(unique(yr))
    
    dts =
        chron(as.vector(outer(mth, yr, function(u, v) paste(u, 1, v, sep = "/"))))

    n         = length(dts)
    friFound = rep(FALSE, n)

    
    knownMths = paste(6:9, "/01/2003", sep="")
    knownMths   = chron(knownMths)
    knownExpDts = c("6/20/2003", "7/18/2003", "8/15/2003", "9/19/2003")
    knownExpDts = chron(knownExpDts)

    ix = dts %in% knownMths
    dts[ix] = knownExpDts[knownMths %in% dts]
    friFound[ix] = TRUE
    
    repeat {

        ind = which("Fri" == as.char(weekdays(dts)) & !friFound)

        friFound[ind] = TRUE
        
        if (n == length(which(friFound)))
            break
        else
            dts[!friFound] = dts[!friFound] + 1

    }

    dts[!ix] = dts[!ix] + 14

    ans = as.char(chronToPosix(dts))
    ans.names = paste(as.char(months(dts)), as.char(years(dts)), sep = "")
    
    ans = as.char(chronToPosix(dts))
    exp.time = paste(exp.hr, exp.min, "00", sep = ":")
    ans = paste(ans, exp.time)
    names(ans) = ans.names

    ans
}





daysToExEquities = function(mth, yr,
t0 = strptime(as.char(Sys.time()), format = "%Y-%m-%d %H:%M:%S"),
exp.hr  = 16,
exp.min = 30)
{
    ##
    ##  Description:  Finds number of calender days until expiration
    ##    of the equity options contract for the given month and year.
    ##

    ex.dts = exDatesEquities(mth, yr,
    exp.hr = exp.hr, exp.min = exp.min)

    tau = as.double(as.POSIXlt(ex.dts) - t0)
    names(tau) = names(ex.dts)

    tau
}







exDatesNG = function(mth, yr, type = "f")
{

    mth  = as.integer(mth[1])
    yr   = as.integer(yr[1])
    type = as.char(tolower(substr(type[1], 1, 1)))
    
    ndays = 0
    ans   = chron(paste(mth, 1, yr, sep = "/"))

    numDaysBack = ifelse("f" == type, 3, 4)
    
    repeat {

        ans = ans - 1

        if (!is.weekend(ans))
            ndays = ndays + 1

        if (ndays >= numDaysBack)
            break
    }

    chronToPosix(ans)
    
}




######################################################################



match.names = function(x = NULL, y = NULL, list = NULL)
    ##
    ##  Description:
    ##    Matches up names of vectors/matrices/data.frames.
    ##    Returns a data.frame
    ##    of common observations or NULL if there are none.
    ##
{
    if (!is.null(x) & !is.null(y))
        list = list(x = as.data.frame(x), y = as.data.frame(y))

    n = length(list)
    
    if (is.null(names(list)) & n)
        names(list) = paste("x", seq(n), sep = "")
    
    ans = as.data.frame(list[[1]])

    
    if (2 <= n) 
        for (i in 2:n) {
            
            tmp       = as.data.frame(list[[i]])
            ind       = intersect(rownames(ans), rownames(tmp))
            x1        = as.data.frame(ans[ind, ])
            names(x1) = names(ans)
            x2        = as.data.frame(tmp[ind, ])
            names(x2) = names(tmp)
            ans       = cbind(x1, x2)
            
        }

    names(ans) = names(list)
    
    ans
}




######################################################################


rm.uri = function(..., list = character(0), envir = .GlobalEnv)
    ##
    ##   Description: Removes all objects in the given environment whose
    ##    name contains the given character strings.
    ##
{

    x = as.character(substitute(list(...)))[-1]
    x = .Primitive("c")(list, x)
    
    for (a in x)
        for (b in ls(envir = envir))
            if (length(grep(a, b))) rm(list = b, envir = envir)

}



rm3 = function() {

    rm.uri(tmp)
    rm.uri(Tmp)
    rm.uri(.warn)

}


detach.pkg = function(name)
{
    name = as.character(substitute(name))

    detach(pos = match(paste("package:", name, sep = ""), search()))
}




reattach.pkg = function(name)
{
    pkgname = as.character(substitute(name))
    
    detach(pos = match(paste("package:", pkgname, sep = ""), search()))
    library(pkgname, character.only = TRUE)
}





######################################################################
## Convenience functions



id       = function(x) x
as.char  = function(x) as.character(x)
is.char  = function(x) is.character(x)
as.posix = function(x) as.POSIXlt(x)


last = function(x, num = 1)
    x[seq(length(x) + 1 - num, length(x))]

last.row = function(x, num = 1)
    x[seq(nrow(x) + 1 - num, nrow(x)), ]


dropLast = function(x, num = 1)
{

    n = length(x)

    stopifnot (num < n)

    x[-seq(n - num + 1, n)]
    
}


randomRows  = function(x, num = 1, ordered = TRUE)
{
    ix = sample(nrow(x), num, rep = FALSE)

    if (ordered)
        ix = sort(ix)
    
    x[ix, ]
    
}



matplot.uri = function(x, y, ...)
    matplot(x, y, type = "b", lty = 1, ...)



difftime.uri = difftime2 = dtime = function(dt1, dt2)
{
    as.integer(round(difftime(as.POSIXlt(dt1), as.POSIXlt(dt2), units = "d")))
}



days.in.month = function(month, year)
{
    year  = year + trunc(month / 12)
    month = (month %% 12) + 1
    as.integer(days(chron(paste(month, 1, year, sep = "/")) - 1))
}


posixToChron = function(dts)
{
    chron(as.character(as.POSIXlt(dts)), format = "y-m-d", out.format = "m/d/y")
}

chronToPosix = function(dts)
{
    strptime(as.character(dts), format = "%m/%d/%y")
}


maxDate = function(dt1, dt2)
{
    
    dt1 = as.POSIXlt(dt1)
    dt2 = as.POSIXlt(dt2)
    ifelse(dt1 > dt2, dt1, dt2)
    
}


minDate = function(dt1, dt2)
{
    
    dt1 = as.POSIXlt(dt1)
    dt2 = as.POSIXlt(dt2)
    ifelse(dt1 < dt2, dt1, dt2)
    
}


sortDates = function(dts)
    chronToPosix(sort(posixToChron(unique(dts))))




seq.char.dts = function(date1, date2)
{

    date1 = posixToChron(date1)
    date2 = posixToChron(date2)

    as.char(chronToPosix(seq(date1, date2)))
    
}



monthnum = function(mth)
    ##
    ##  Description: Converts character string mth representing
    ##    a calender month in full or abbreviated form
    ##    to the corresponding integer 0-11.
    ##
{
    mth = substr(mth, 1, 3)
    
    strptime(paste(mth, 1, 2001, sep = "/"), format = "%b/%d/%Y")$mon
}

######################################################################

varRatio = function(x, maxOrder)
{

    .Call("varRatio", x, maxOrder, package = "util.uri")
    
}



######################################################################


cum.mean = function(x)
{
    cumsum(x) / seq(length(x))
}

cum.var = function(x, na.rm = TRUE, reverse = FALSE)
{
    if (na.rm)
        x = x[!is.na(x)]

    if (reverse)
        x = rev(x)
    
    ans = cumsum(x^2) / seq(length(x)) - (cumsum(x) / seq(length(x)))^2

    if (reverse)
        ans = rev(ans)

    ans
}

cum.moment = function(x, n)
{
    cumsum(x^n) / seq(length(x))
}

###################################################################




rollApply = function(y, fn, order, ...)
{
    
    .Call("roll_apply",
          as.double(y),
          function(x) fn(x, ...),
          as.integer(order),
          new.env(),
          PACKAGE = "uri")
    
}



movAvg = function(y, order) {
    ans = .C("m_avg", as.double(y), as.integer(order),
    as.integer(length(y)), ans = double(length(y)),
    NAOK = TRUE,
    PACKAGE = "uri")$ans

    names(ans) = names(y)

    ans
}



movMax = function(y, order) {
    ans = .C("m_max", as.double(y), as.integer(order),
    as.integer(length(y)), ans = double(length(y)),
    NAOK = TRUE,
    PACKAGE = "uri")$ans

    names(ans) = names(y)

    ans
}



movMin = function(y, order)
{
    
    ans = .C("m_min", as.double(y), as.integer(order),
    as.integer(length(y)), ans = double(length(y)),
    NAOK = TRUE,
    PACKAGE = "uri")$ans

    names(ans) = names(y)

    ans
}



seq.dates.uri = function(from, to, by = "days", length, ...) {
    
    ans = seq.dates(posixToChron(from), posixToChron(to), by, length, ...)
    strptime(as.character(ans), format = "%m/%d/%y")
    
}

subscript = function(x, index) x[[index]]

between = function(y, x) {
    ## returns i s.t. x[i] <= y < x[i + 1]
    ## NB output has same length as y.
    ifelse(max(x) <= y, length(x), apply(outer(x, y, "<="), 2, which.min) - 1)
}


#####################################################################


circle = function(f, g) function(x) f(g(x))



######################################################################
## Functions related to the vol run-up stuff I did.



localVar2 = function(t0, t1, vol, a = 67.87, alpha = 0.77) {
    ## t0, t1 are in calender days
    t1 = t1 / 365; t0 = t0 / 365; a = a / 365
    C = vol^2*(t1 + a)^alpha
    C * {(1 - alpha) * (t1 - t0) + a}*{t1 - t0 + a}^{-alpha - 1}
}



cumVar2 =  function(t0, t1, vol, a = 67.87, alpha = 0.77) {
    t1 = t1/365; t0 = t0 / 365; a = a/365
    C = vol^2 * (t1 + a)^alpha
    t1 * C * (t1 + a)^-alpha - (t1 - t0) * C * (t1 - t0 + a)^-alpha
}



vol2fn =  function(t0, t1, vol, a = 67.87, alpha = 0.77) {
    sqrt(cumVar2(t0, t1, vol, a, alpha) / t1 * 365)
}

######################################################################

which.na     = function(x) which( is.na(x))
notwhich.na  = function(x) which(!is.na(x))
which.not.na = function(x) which(!is.na(x))

######################################################################







optionPayoff = function(S, K, posn = 1, optionType = "c") {
    ##
    ##  Description: Payoff functions for various option positions.
    ##

    len = max(length(S),length(K),length(posn),length(optionType))
    if (1 == length(S)) S = rep(S, len)
    if (1 == length(K)) K = rep(K, len)
    if (1 == length(posn)) posn = rep(posn, len)
    if (1 == length(optionType)) optionType = rep(optionType, len)
    
    stopifnot (len == min(length(S),length(K),length(posn),length(optionType)))

    y = S-K

    payoffs = ifelse(optionType=="c", pmax(y,0), rep(NA,len))
    payoffs = ifelse(optionType=="p", pmax(-y,0), payoffs)
    payoffs = ifelse(optionType=="s", abs(y), payoffs)
    
    ans = posn * payoffs

    if (length(names(S))) {
        names(ans) = names(S)
    } else if (length(names(K))) {
        names(ans) = names(K)
    } else if (length(names(posn))) {
        names(ans) = names(posn)
    } else if (length(names(optionType))) {
        names(ans) = names(optionType)
    }

    ans
}




######################################################################


delta.to.strike = dtoK = function(tau, S, vol,
                  delta, optionType = "c",
                  r = 0, q = r) {

    .Call("deltaToStrike", tau, S, vol, delta, optionType, r, q)

}






###########################################################################


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

    if (0 == K) {

        ans = margrabe(tau, S1, S2, vol1, vol2, rho, option.type = option.type)
        ans = ans[, c("p", "d1", "d2")]
        
    } else {

        ans = .Call("margrabe2",
        tau, S1, S2, K, vol1, vol2, rho,
        r, q1, q2, option.type, paths, eps,
        PACKAGE = "uri")
        
    }

    if (!is.null(return.type))
        ans[, return.type]
    else
        ans
}







#############################################################################

bisectMethod = function(fn, lower, upper, tol = 1e-3, maxit = 1e2) {
    ## Last modified 3-Dec-2002.
    ## Find zeros using bisection method

    if (0 == fn(lower))
        x = lower
    else if (0 == fn(upper))
        x = upper
    else {
        stopifnot (sign(fn(lower)) == -sign(fn(upper)))

        it = 0
        
        repeat {
            x = 0.5 * (lower + upper)
            val = fn(x)
            it = it + 1

            if (abs(val) < tol | it >= maxit)
                break
            
            if (sign(val) == sign(fn(upper)))
                upper = x
            else if (sign(val) == sign(fn(lower)))
                lower = x
            else
                break
        }
    }

    list(val = fn(x), par = x)
}





implied.vol = impvol = function(tau, S, K, price,
              r = 0,
              q = r,
              option.type = "c",
              tol = 1e-6,
              maxit = 1e2)
{

    .Call("imp_vol", tau, S, K, price, r, q, option.type, tol, maxit,
          PACKAGE = "uri")
    
}



#############################################################################

sd.ml = function(x, na.rm = TRUE)
    ##
    ##  Description: Returns maximal likelihood
    ##    standard deviation.
    ##
{
    sqrt(mean((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm))
}


var.ml = function(x, na.rm = TRUE)
    ##
    ##  Description: Returns maximal likelihood
    ##    variance.
    ##
{
    mean((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm)
}




var2 = function(x, nlags = 1, na.rm = TRUE) {

    if (na.rm)
        x = x[!is.na(x)]

    n = length(x)
    
    var.ml(x) / (1 - 1/n - 2/n^2 * sum(autocor(x, 1:nlags)))
    
}



tstat = function(x, na.rm = TRUE)
{

    if (na.rm) x = x[!is.na(x)]
    
    mean(x) / sd(x) * sqrt(length(x) - 1)
    
}




sd.boot = function(f, x, R = 1e3 - 1, na.rm = TRUE, loop = FALSE)
    ##
    ##  Description: Returns the bootstrapped standard deviation
    ##    of f(x).  f may be vector-valued.  x must be a vector.
    ##    Set loop = TRUE if allocation of a vector of
    ##    length R * length(x) causes memory problems.
    ##
{

    stopifnot (is.vector(x))
    
    if (na.rm) x = x[!is.na(x)]
    
    if (loop) {

        n   = length(f(x))
        
        ans = matrix(NA, R, n)
        
        for (b in seq(R))
            ans[b, ] = f(sample(x, rep = TRUE))

        apply(ans, 2, sd)
        
    } else {

        m      = length(x)
        x      = sample(x, R * m, rep = TRUE)
        dim(x) = c(m, R)
        x      = apply(x, 2, f)

        if (is.vector(x))
            sd(x)
        else
            apply(x, 1, sd)
        
    }

}



boot.dist = function(f, x, R = 1e3 - 1, na.rm = TRUE,
summary.only = TRUE, loop = FALSE)
    ##
    ##  Description: Returns the bootstrapped distribution of f(x) if
    ##    !summary.only, else the summary of the boostrapped
    ##    distribution of f(x).  In this function f should be
    ##    scalar-valued.  Use loop = TRUE if R cannot allocate a vector
    ##    of size R * length(x).
    ##
{
    stopifnot (is.vector(x))
    
    if (na.rm) x = x[!is.na(x)]
    
    stopifnot (1 == length(f(x)))

    if (loop) {

        fstar = rep(NA, R)

        for (b in seq(R))
            fstar[b] = f(sample(x, rep = TRUE))
        
    } else {

        m      = length(x)
        x      = sample(x, m * R, rep = TRUE)
        dim(x) = c(m, R)

        fstar  = apply(x, 2, f)
        
    }

    if (summary.only) {
        
        ans    = quantile(fstar, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
        ans[8] = mean(fstar)
        ans[9] = sd(fstar)

        names(ans) =
            c("min", "0.05", "Q1", "median", "Q3", "0.95", "max", "mean", "sd")
        
    } else {

        ans = fstar
        
    }


    ans
}




skew = function(x, na.rm = TRUE) {
    ## Last modified 12-Jan-2002.

    if (na.rm)
        x = x[!is.na(x)]
    
    mean((x - mean(x))^3) / sd.ml(x)^3
    
}


skew2 = function(x) {
    ## Last modified 6/14/02.
    x = x[ !is.na(x) ]
    m = mean( x )
    s = sd( x )
    n = length(x)
    n * sum( ( (x - m)^3 ) /( s^3 ) ) / (n-1) / (n-2)
}


kurt = function(x, na.rm = TRUE, type = "excess") {
    ## Last modified 12-Jan-2002.

    if (na.rm)
        x   = x[!is.na(x)]
    
    ans = mean((x - mean(x))^4) / var.ml(x)^2
    ifelse(type == "excess", ans - 3, ans)
    
}


kurt2 = function(x, type = "excess") {
    ## Last modified 6/14/02.
    x = x[ !is.na(x) ]
    m = mean( x )
    s = sd( x )
    n = length(x)
    ans = n*(n+1)*sum( ( (x-m)^4 )/( s^4 ) ) / (n-1) / (n-2) / (n-3 )
    if ( type != "total" ) ans = ans - 3*(n-1)^2/(n-2)/(n-3)
    ans
}




modeStat = function(x, na.rm = TRUE)
{
    if (na.rm)
        x = x[!is.na(x)]
    
    ux = unique(x)
    numOccurrences = apply(outer(x, ux, "=="), 2, sum)

    maxNumOccurrences = max(numOccurrences)
    modeIndex = which(numOccurrences %in% maxNumOccurrences)

    ux[modeIndex]
}



autocor = function(x, n = 1, use = "p") {
    ## Last modified 03-Feb-2003.

    X = x
    for (j in n)
        X = cbind(X, lag2(x, j))
    
    ans = cor(X, use = use)[1, -1]
    names(ans) = paste("lag", n)
    ans
}



#########################################################################

lag2 = function(x, n = 1, MARGIN = 2)
    ## Last modified 7-Oct-2002.
    ## Returns the time-series x lagged by n.  Thus lag2(x, n)[i] = x[i - n].
    ## Missing values are placed where appropriate.
{

    if (is.matrix(x) | is.data.frame(x)) {

        y = apply(x, MARGIN = MARGIN, function(u) lag2(u, n = n))
        
        if (1 == MARGIN)
            y = t(y)

        dimnames(y) = dimnames(x)
        
    } else {
        
        M = length(x)
        y = rep(NA, M)
        
        if (n >= 0 & n < M)
            y[(n + 1):M] = x[1:(M - n)]
        else if (n < 0 & abs(n) < M)
            y[1:(M + n)] = x[(-n + 1):M]
        
        names(y) = names(x)

    }
    
    y
    
}


############################################################################


standard = function(x, na.rm = TRUE)
{
    ## Last modified 6/14/02.
    if (na.rm) x = x[!is.na(x)]
    ans = as.vector(scale(x))
    names(ans) = names(x)
    ans
}


center = function(x, na.rm = TRUE)
{
    
    if (na.rm)
        x = x[!is.na(x)]
    
    x - mean(x, na.rm = TRUE)
    
}


fisher = function(x, inverse = FALSE)
    ##
    ##  Description: Returns fisher transform of x.
    ##
{
    if (inverse)
        ans = {exp(2 * x) - 1}/{exp(2 * x)+1}
    else
        ans = log((1 + x)/(1 - x))/2
    
    names(ans) = names(x)
    ans
}



sort.uri = function(x, y) {
    ## Last modified 7/06/02.
    ## Sorts y then applies the corresponding permuation to x.
    x[sort(y, index = TRUE)$ix]
}


unif.uri = function(n, jiggle = TRUE, permute = FALSE, rep = FALSE)
    ## Last modified 2003 Feb 19.
{
    
    .Call("unifuri", n, jiggle, permute, rep, PACKAGE = "uri")

}


norm.uri = function(n = 100, jiggle = TRUE, permute = FALSE, rep = FALSE)
    ## Last modified 2003 Feb 19.
{
    .Call("normuri", n, jiggle, permute, rep, PACKAGE = "uri")
}


pseudo.inverse = function(x) {
    ## Last modified 3-Dec-2002
    stopifnot(is.matrix(x))
    ans = lm(diag(dim(x)[1]) ~ 0 + x)$coef
    dimnames(ans) = dimnames(x)
    ans
}


sample.uri = function(n, size = 1, rep = FALSE, sorted = TRUE)
{

    ans = unif.uri(size, jig = TRUE, permute = TRUE, rep = rep)
    ans = round((n - 1) * ans + 0.5)
    
    if (sorted)
        ans = sort(ans)

    ans
}



mvtunif.uri = function(n, d = 2, jiggle = TRUE, permute = TRUE,
replace = FALSE)
{

    if (d == 1) {

        unif.uri(n, jig = jiggle, rep = replace, perm = permute)    
        
    } else if (d >= 2) {

        ans = x = seq(0.5 / n, 1 - 0.5 / n, length = n)

        for (j in seq(2, d)) {
            
            ans = apply(as.matrix(ans), 2, rep, times = n)
            tmp = rep(x, times = rep(n^(j - 1), n))
            ans = cbind(ans, tmp)
            
        }

        if (jiggle)
            ans = ans + (1/n) * matrix(runif(n^d * d, -0.5, 0.5), n^d, d)

        if (permute) 
            ans = ans[sample(nrow(ans), rep = replace), ]

        colnames(ans) = NULL
        
        ans
    } else NA
}




bvtnorm.uri = function(n, rho = 0.0,
jiggle = TRUE, permute = FALSE, replace = FALSE)
{

    ans = .Call("bvt_normuri", n, rho, jiggle, PACKAGE = "uri")

    if (permute)
        ans = ans[sample(nrow(ans), rep = replace), ]

    ans
    
}


mvtnorm.uri = function(n, d = 2, jiggle = TRUE,
permute = TRUE, replace = FALSE)
{
    qnorm(mvtunif.uri(n, d, jiggle, permute, replace))
}


vector.norm = function(x) {sqrt(sum(x^2))}




sim.gbm = function(vol, nticks, freq = 365, paths = 1, S = 1,
mu = 0, terminal.only = FALSE)
    ##
    ##  Description: Returns an n x path matrix if terminal.only is
    ##    FALSE and a vector of length = paths otherwise.  nticks
    ##    is assumed to be the # of ticks from the time horizon.
    ##    freq is the frequency of ticks per time period that
    ##    mu and vol are quoted in.  The default setting basically
    ##    corresponds to one tick == one day.
    ##
{

    if (terminal.only) {
        
        tau   = nticks / freq
        e     = rnorm(paths, m = (mu - 0.5 * vol^2) * tau, s = vol * sqrt(tau))
        S * exp(e)
        
    } else {

        dt = 1.0 / freq;
        m  = (mu - 0.5 * vol^2) * dt;
        s  = vol * sqrt(dt);
        
        e      = rnorm(nticks * paths, m = m, s = s)
        dim(e) = c(nticks, paths)
        drop(S * exp(apply(e, 2, diffinv)))
        
    }
    
}



############################################################################

list.union = function(list)
{
    stopifnot(is.list(list))
    
    ans = unique(as.vector(list[[1]]))

    n = length(list)

    if (2 <= n)
        for (i in seq(2, n)) ans = union(ans, as.vector(list[[i]]))

    ans
    
}


list.intersect = function(list)
{
    
    stopifnot(is.list(list))
    
    ans = unique(list[[1]])

    n = length(list)

    if (2 <= n)
        for (i in seq(2, n)) ans = intersect(ans, list[[i]])

    ans
    
}



## This is similar to the R "stack" function.
join.list = function(list, join.names = TRUE)
{
    stopifnot(is.list(list))
    ans = NULL

    n = length(list)

    if (n > 1 & join.names & is.null(names(list)))
        names(list) = paste("L", seq(n), sep = "")
    
    if (n)
        for (i in seq(1, n)) {
            
            if (join.names) {
                
                if (is.null(names(list[[i]])))
                    names(list[[i]]) = paste("x", seq(length(list[[i]])), sep = "")
                
                names(list[[i]]) = paste(names(list)[i], names(list[[i]]), sep = ".")
            }
            
            ans = c(ans, list[[i]])
        }
    
    ans
}



## Stacks lists of matrices or data.frames
stack.matrix = function(list)
{
    n = length(list)

    ans = list[[1]]

    if (n > 1)
        for (i in 2:n) {

            ans = rbind(ans, list[[i]])
            
        }

    ans
}





######################################################################


chol.pseudo = function(X,
method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
num.init = 20)
{
    ## Last modified 7/6/02.
    ## Finds square upper-triangular matrix A that
    ## minimizes the squared distance between X and t(A)%*%A.

    stopifnot(all(is.matrix(X), dim(X)[1] == dim(X)[2]))
    options(show.error.messages = FALSE)
    A = try(chol(X))
    if (is.numeric(A))
        val = 0
    else {
        n = nrow(X)
        num.init = ifelse(num.init[1] < 1, 1, round(abs(num.init[1])))
        ind.upper = upper.tri(X, diag = TRUE)
        
        fcn0 = function(A.upper) {
            A1 = matrix(0, n, n)
            A1[ind.upper] = A.upper
            dim(A1) = dim(X)
            sum((t(A1) %*% A1 - X)^2)
        }
        
        val = Inf
        A =  matrix(0, n, n)
        
        for(k in seq(num.init)) {

            par.init = rnorm(0.5 * n * (n + 1))
            tmp = optim(par.init, fcn0, method = method[1])
            val.new = tmp$val
            
            if(val.new < val) {
                A[indices.upper] = tmp$par
                val = val.new
            }
            
        }
    }

    list(matrix = A, distance = sqrt(val))
    
}

######################################################################


## This computes a single sharpe ratio
sharpe.ratio = function(x, na.rm = TRUE)
{
    ## Last modified 03-Feb-2003.

    mean(x, na.rm = na.rm) / sd(x, na.rm = na.rm)
}

