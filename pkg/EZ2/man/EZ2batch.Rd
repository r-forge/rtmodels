\name{EZ2batch}
\alias{EZ2batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Batch EZ2 model fitting }
\description{
  Fit the a simplified diffusion model for response time and accuracy in batch to each row of a data frame.
}
\usage{
EZ2batch(pstart, ObsValPair, ..., data, nrestart = 1, method = "Nelder-Mead", control = list(), hessian = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pstart}{ Numeric vector with \emph{named} elements. Starting point. }
  \item{ObsValPair}{ Formula or list of formulas relating observed value(s) to parameters. See examples below. }
  \item{\dots}{ More formulas. }
  \item{data}{ Data frame with observed values. }
  \item{nrestart}{ Number of times the fitting should restarted. Defaults to 1. }
  \item{method}{ See \code{method} parameter of \code{\link{optim}}. }
  \item{control}{ See \code{control} parameter of \code{\link{optim}}. }
  \item{hessian}{ See \code{hessian} parameter of \code{\link{optim}}. }
}
\details{
  The number of restarts needs to be increased if the model entails the same number of 
  estimated parameters as the number of modeled moments and any of the values in the
  \code{value} column in the returned data frame is not close to \code{1.0e-08}.
  
  \emph{An important bug to note: for EZbatch to fit the rows in the provided data.frame there should be no variables in the global environment with the same name as the column names of the data.frame.}
}
\value{
%  ~Describe the value returned
   Data frame with parameter estimates and convergence information.
}
\references{ Grasman, Wagenmakers, and van der Maas (2008) EZ2: An extention of the EZ-diffusion model for Response Time and Accuracy, \emph{Manuscript submitted for publication in J. Math. Psych.} }
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{EZ2}} }
\examples{
## create some data (theoretical values, not simulated)
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }