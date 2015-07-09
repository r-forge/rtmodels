`EZ2batch` <-
function (pstart, ObsValPair, ..., data, nrestart = 1, method = "Nelder-Mead", 
    control = list(), hessian = FALSE) 
{
    mdl = c(ObsValPair, ...)
    t(apply(data, 1, function(x) {
        attach(as.list(x), warn.conflicts = FALSE) 
        # with(as.list(x), { # This was an attempt to remove the attach statement... but it lead to an execution ERROR in R CMD check
          fit <- EZ2(pstart, mdl, method = method, control = control, 
              hessian = hessian)
          for (i in 1:nrestart) fit <- EZ2(fit$par, mdl, method = method, 
              control = control, hessian = hessian)
          unlist(fit)
        # })
    }))
}

