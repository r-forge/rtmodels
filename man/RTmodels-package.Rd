\name{RTmodels-package}
\alias{RTmodels-package}
\docType{package}
\title{
Utility functions for simulation of reaction times from Ratcliff's diffusion model and computing associated probabilities.
}
\description{
A set of functions for simulating RTs from the diffusion model as proposed by Ratcliff,
for calculating associated probabilities, and for estimating parameters.
}
\details{
\tabular{ll}{
Package: \tab difsim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-05-19\cr
License: \tab GPL v2 or higher \cr
}
  rt = difsim(v=.2,a=.2,z=.09,Ter=.250,N=1000);
  quantile(rt[[1]],seq(0.05,0.95,l=10))
}
\author{
Raoul Grasman
Eric-Jan Wagenmakers

Maintainer: Raoul Grasman <rgrasman@uva.nl>, Eric-Jan Wagenmakers <e.j.m.wagenmakers@uva.nl>
}
\references{
Ratcliff & Tuerlinckx (2002) Estimating parameters of the diffusion model: Approaches to dealing with contaminant reaction times and parameter variability. Psychonomic Bulletin & Review, vol. 9 (3), p. 438--481.
}
\keyword{ package }
\examples{
  rt = difsim(v=.2,a=.2,z=.09,Ter=.250,N=1000);
  quantile(rt[[2]],seq(0.05,0.95,l=10))

  prt.error(.Last.value-.250,v=.2,a=.2,z=.09) # assumes Ter=0
}
