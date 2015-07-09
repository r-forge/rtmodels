`EZ2` <-
function (pstart, ObsValPair, ..., method = "Nelder-Mead", control = list(), 
    hessian = FALSE) 
{
    tab = c(ObsValPair, ...)
    if (is.null(names(pstart)) | any(names(pstart) == "")) 
        stop("All elements of pstart should have names.")
    optim(pstart, EZ2.objf, method = method, control = control, 
        hessian = hessian, nms = names(pstart), ObsValPairs = tab, 
        grad = FALSE)
}

