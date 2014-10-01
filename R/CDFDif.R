

CDFDif = function(t, x, nu, z, a, Ter, eta=0, sz=0, st=0){
# purpose: compute probability Pr(T < t, X=x) in Ratcliff's diffusion model
# arguments:
#       a   - the boundary separation
#       Ter - the mean non-decisional component time
#       eta - the standard deviation of the normal drift distribution
#       z   - the mean starting point
#       sz  - the spread of the starting point distribution
#       st  - the spread of the non-decisional component time distribution
#       nu  - the mean drift rate
# output:
#       the value of the cumulative distribution function evaluated at t
#       and x for parameters Par, F_{X,T}(x,t) = Pr(X=x,T<=t)
#
# remarks:  Fortran source by Francis Tuerlinckx, University of Leuven, Belgium
# see:      Tuerlinckx (2004) The efficient computation of the cumulative
#           distribution and probability density functions in the diffusion
#           model, _Behavior Research Methods, Instruments, & Computers_,
#           vol. 34 (4), p. 702--716.
    eps = .Machine$double.eps
    sqrteps = sqrt(eps)
    if (any(c(length(z), length(a), length(Ter), length(eta),
        length(sz), length(st)) != 1))
        warning("This function is not vectorized. ", ifelse(length(t) >
            1, "Use cdfdif for a version of this function that is vectorized with respect to the parameter t",
            ""))
    if (any(t < 0 | z < 0 | a < 0 | z > a | Ter < 0 | eta < 0 |
        sz < 0 | st < 0))
        stop("Invalid input. Parameters z, a, Ter, eta, sz, and st should all be non-negative. Parameter z should be less than parameter a.")
    if (sz < eps) {
        sz = sqrteps
        warning("Parameter sz (", deparse(substitute(sz)), ") set to minimum value ",
            sqrteps)
    }
    if (st < eps) {
        st = sqrteps
        warning("Parameter st (", deparse(substitute(st)), ") set to minimum value ",
            sqrteps)
    }
    if (eta < eps) {
        eta = sqrteps
        warning("Parameter eta (", deparse(substitute(eta)),
            ") set to minimum value ", sqrteps)
    }
    .C("cdfdif", as.double(t[1]), as.double(x[1]), as.double(c(a[1],
        Ter[1], eta[1], z[1], sz[1], st[1], nu[1])), CDF = double(1))$CDF
}

cdfdif = function(t, x, nu, z, a, Ter, eta=0, sz=0, st=0){
# purpose: compute probability Pr(T < t, X=x) in Ratcliff's diffusion model
# arguments:
#       t   - array of repsonse times in seconds
#       x   - logical indicating wether error (=FALSE) or correct (=TRUE) RTs are provided
#       a   - the boundary separation
#       Ter - the mean non-decisional component time
#       eta - the standard deviation of the normal drift distribution
#       z   - the mean starting point
#       sz  - the spread of the starting point distribution
#       st  - the spread of the non-decisional component time distribution
#       nu  - the mean drift rate
# output:
#       the value of the cumulative distribution function evaluated at t
#       and x for parameters Par, F_{X,T}(x,t) = Pr(X=x,T<=t)
#
# remarks:  Fortran source by Francis Tuerlinckx, University of Leuven, Belgium
# see:      Tuerlinckx (2004) The efficient computation of the cumulative
#           distribution and probability density functions in the diffusion
#           model, _Behavior Research Methods, Instruments, & Computers_,
#           vol. 34 (4), p. 702--716.
    if (any(c(length(z), length(a), length(Ter), length(eta),
        length(sz), length(st)) != 1))
        warning("This function is not vectorized with respect to z, a, Ter, eta, sz, or st.")
    if (any(t < 0 | z < 0 | a < 0 | z > a | Ter < 0 | eta < 0 |
        sz < 0 | st < 0))
        stop("Invalid input. Parameters z, a, Ter, eta, sz, and st should all be non-negative. Parameter z should be less than parameter a.")
    eps = .Machine$double.eps
    sqrteps = sqrt(eps)
    if (sz < eps) {
        sz = sqrteps
        warning("Parameter sz (", deparse(substitute(sz)), ") set to minimum value ",
            sqrteps)
    }
    if (st < eps) {
        st = sqrteps
        warning("Parameter st (", deparse(substitute(st)), ") set to minimum value ",
            sqrteps)
    }
    if (eta < eps) {
        eta = sqrteps
        warning("Parameter eta (", deparse(substitute(eta)),
            ") set to minimum value ", sqrteps)
    }
    .C("cdfdifn", as.double(t), as.integer(length(t)), as.double(x[1]),
        as.double(c(a[1], Ter[1], eta[1], z[1], sz[1], st[1],
            nu[1])), CDF = double(length(t)))$CDF
}
