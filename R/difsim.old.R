
"difsim.old" <-
function (v, a, N, z = 0.5 * a, s = 0.1, var.v = 0, range.z = 0, 
    Ter = 0, range.Ter = 0, stepsize=0.0)
{

    res = .C("difsim", as.double(v), as.double(a), as.double(z), as.double(s), 
	double(1), double(1), as.double(Ter), double(1), as.integer(N), err=integer(1), 
	t.top=double(N), t.bot=double(N))
    if (res$err != 0)
        warning("An error occured.")
    return(list(t.top = res$t.top[res$t.top > 0], t.bot = res$t.bot[res$t.bot >
        0]))
}

