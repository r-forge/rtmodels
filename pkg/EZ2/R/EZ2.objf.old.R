`EZ2.objf.old` <-
function(p,nms,ObsValPairs,grad=FALSE)
{
    #
    # the basic idea underlying this function is:
    #   ObsValPair = c(func="EZ2.vrt", obs=0.3,v= ~v, z=~a-z, a=~a)
    #   args = lapply(ObsValPair[-(1:2)],function(exp) eval(exp[[length(exp)]]))
    #   do.call("EZ2.vrt",args)
    #
###   eval(parse(text=paste(nms,'=',p,collapse=';'))) # define the parameters ---> changed to using 'with'
    names(p)=nms # nlm removes names 
    p=as.list(p)
    f = sapply(ObsValPairs, #1, 
		function(x) do.call(x[[1]],lapply(x[-(1:2)],function(expr) with(p,eval(expr[[length(expr)]]))))
        )
    e = unlist(sapply(ObsValPairs,function(x)x[[2]]))-unlist(f) #unlist(ObsValPairs[,2])-unlist(f)
    s = sum(e^2)
    if(grad & FALSE){ # this is currently not used...
	.grad = -2*t(df)%*%e
	attr(s,"gradient") = .grad
    }
    return(s)
}

