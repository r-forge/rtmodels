`EZ2.objf` <-
function(p,ObsValPairs,nms=names(p),grad=FALSE)
{
    #
    # the basic idea underlying this function is:
    #   ObsValPair = 0.3 ~ EZ2.vrt(nu=nu,z=a-z,a=a)
    #   with(list(nu=.1,z=0.09,a=0.2), eval(ObsValPair[[3]])
    #
    names(p) = nms
    p = as.list(p)
    f = sapply(ObsValPairs, function(.x) {
        y = with(p, eval(.x[[3]]))
        z = rep(0, length(nms))
        names(z) = nms
        .g = attr(y, "gradient")[1, ]
        z[names(.g)] = .g
        c(pred = y, z)
    })
    e = sapply(ObsValPairs, function(x) eval(x[[2]])) - f[1, ]
    s = 1e+06 * sum(e^2)
    if (grad) {
        df = 1e+06 * f[-1, ]
        .grad = -2 * df %*% e
        attr(s, "gradient") = drop(.grad)
    }
    return(s)
}

