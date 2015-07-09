\name{EZ2-package}
\alias{EZ2-package}
\docType{package}
\title{
   EZ2 diffusion modeling of Response Time and Accuracy
}
\description{
   EZ2 fits a simplified version of the diffusion model for human and monkey response times 
   and accuracy to the means and variances of the observables. 
}
\details{
\tabular{ll}{
Package: \tab EZ2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-08-31\cr
License: \tab GPL version 2 (or later)\cr
}
You need a data frame containing (at least) the moments (RT means, RT variances, and proportions correct) that you would like to model.

Specify for each used moment a formula of the form \code{vrt1 ~ EZ2.vrt(v1, z, a)}.

Supply the model with the data and a startingpoint for all unknown parameters to \code{EZ2} or \code{EZ2batch} to calculate method of moment estimators or least squares estimators of the unknown parameters in the model (that is, \code{v1}, \code{z}, and \code{a} in \code{vrt1 ~ EZ2.vrt(v1, z, a)}).

Your model can be complex and extensive; for instance 

\code{vrt1 ~ EZ2.vrt(v1, z, a) + p0*(maxRT-minRT)^2/12+p0*(1-p0)*(EZ2.mrt(v1, z, a)-(maxRT+minRT)/2)^2}

is equaly valid.
}

\author{
Raoul Grasman

Maintainer: Raoul Grasman <rgrasman@uva.nl>
}
\references{
Grasman, Wagenmakers & van der Maas (2007) "EZ2: An extention of the EZ-diffusion model for Response Time and Accuracy, \emph{Manuscript submitted for publication in J. Math. Psych.}
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%~~ simple examples of the most important functions ~~
## create some data (theoretical values, not simulated) Needless to say, in reality you would like to fit real data!
A = seq(.08,.13,len=6)
X2 = data.frame(A=A)
X2$vrt0 = sapply(A, function(a) EZ2.vrt(.1,.05,a))
X2$pe0 = sapply(A, function(a) EZ2.pe(.1,.05,a))
X2$vrt1 = sapply(A, function(a) EZ2.vrt(.2,a-.05,a))
X2$pe1 = sapply(A, function(a) EZ2.pe(.2,a-.05,a))

X2 = as.data.frame(X2)          # now pretend that X2 is the data frame that you may have computed from real data

## fit an EZ2 model on each row
#  method 1:
EZ2batch(c(v0=.11,v1=.21,z=.05,a=.09), 
 vrt0 ~ EZ2.vrt(v0,z,a), 
  pe0 ~ EZ2.pe(v0,z,a), 
 vrt1 ~ EZ2.vrt(v1,a-z,a), 
  pe1 ~ EZ2.pe(v1, a-z, a), data=X2)

# method 2 (eventually less typing):
mdl <- list( vrt0 ~ EZ2.vrt(v0,z,a), 
              pe0 ~ EZ2.pe(v0,z,a), 
             vrt1 ~ EZ2.vrt(v1,a-z,a), 
              pe1 ~ EZ2.pe(v1, a-z, a)
           )
EZ2batch(c(v0=.11,v1=.21,z=.05,a=.09), mdl, data=X2)
}
\keyword{ models }
\keyword{ package }