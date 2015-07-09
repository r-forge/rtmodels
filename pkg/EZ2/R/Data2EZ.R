`Data2EZ` <-
function (Pc,VRT,MRT,s=0.1)
# input:
#				Pc - Proportion correct
#				VRT - sample variance of the RT's
#				MRT - sample mean of the RT's
#				s - diffusion standard deviation
# returns:	Object with properties v, a, and Ter, containing EZ-estimates of these parameters
{
	#var L, x, v, a, y, MDT, Ter, s2, arr = [];
	p = Pc; 
      pow = function(a,b) a^b;
      logit = function(p) log(p/(1-p))
	if (p == 0)
	{
		stop("Oops, only errors!");
	}
	if (p == 0.5)
	{
		stop("Oops, chance performance!");
	}
	if (p == 1)
	{
		stop("Oops, only correct responses!");
	}
	s2     = s*s;
	L      = logit(p);
	x      = L*(L*p*p - L*p + p - 0.5) / VRT;
	v      = sign(p-0.5)*s*pow(x,(1/4));
	a      = s2*logit(p)/v;
	y      = -v*a/s2;
	MDT    = (a/(2*v))*(1-exp(y))/(1+exp(y));
	Ter    = ifelse(!missing(MRT),MRT - MDT,NA); # compute Ter only if MRT was provided
	return(list(v= v, a= a, Ter= Ter));
}

