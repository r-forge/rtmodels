\name{EZ2}
\alias{EZ2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit diffusion model to observed sample moments }
\description{
  Fit the a simplified diffusion model for response time and accuracy to observed proportions of errors and response time means and variances.
}
\usage{
EZ2(pstart, ObsValPair, ..., method = "Nelder-Mead", control = list(), hessian = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pstart}{ Vector with named starting parameter values }
  \item{ObsValPair}{ Observed-predicted value pair in the form \code{0.80 ~ EZ2.pe(v1, z, a) } or \code{vrt2 ~ EZ2.pe(v2, a-z, a) if vrt2 exists in the global environment and is numeric } }
  \item{\dots}{ More Observed-predicted value pairs (as many as desired, see example) }
  \item{method}{ See \code{method} parameter of \code{\link{optim}}. }
  \item{control}{ See \code{control} parameter of \code{\link{optim}}. }
  \item{hessian}{ See \code{hessian} parameter of \code{\link{optim}}. }
}
\details{
   EZ2 fits a simplified version of the diffusion model for human and monkey response times 
   and accuracy to the means and variances of the observables. This model of information accumulation and decision is a simplified version of Ratcliff's diffusion model (1978). 
   
   Use \code{\link{EZ2batch}} for more user friendly fitting automatically each row in a \code{data.frame}.
}
\value{
%  ~Describe the value returned
  The list returned by optim. This list contains
  \item{par}{A vector containing the estimates}
  \item{value}{Sum of squared prediction errors. This should be very close to zero (order of \code{1e-8}) if there are as many observed-predicted moment value pairs as there are unknown parameters (the estimates then constitute method of moments estimators).}
  For other list members, see \code{Value} section of \code{\link{optim}} for a complete description.
}
\references{ 
Ratcliff, R. (1978). Theory of Memory Retrieval. \emph{Psychological review} vol. 85 (2) pp. 59-108 

Grasman, R. P. P., Wagenmakers, E.-J., & van der Maas, H. L. J. (2007). On the mean and variance of response times under the diffusion model with an application to parameter estimation, \emph{J. Math. Psych.} 53: 55--68.
}
\author{ Raoul P. P. P. Grasman }
\seealso{ \code{\link{EZ2-package}}, \code{\link{EZ2batch}} }
\examples{
 Vrt0 <- EZ2.vrt(0.1, 0.08, 0.12)
  Pe0 <- EZ2.pe(0.1, 0.08, 0.12)
 Vrt1 <- EZ2.vrt(0.15, 0.12-0.08, 0.12)
  Pe1 <- EZ2.pe(0.15, 0.12-0.08, 0.12)
  
 EZ2(c(v0=.156,v1=.0844,z=.0631,a=.1263), 
 Vrt0 ~ EZ2.vrt(v0,z,a), 
  Pe0 ~ EZ2.pe(v0,z,a), 
 Vrt1 ~ EZ2.vrt(v1,a-z,a), 
  Pe1 ~ EZ2.pe(v1, a-z, a))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
