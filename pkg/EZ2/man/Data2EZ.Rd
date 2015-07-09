\name{Data2EZ}
\alias{Data2EZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert observed sample moments to diffusion parameter values }
\description{
%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Converts proportion correct (Pc), response time variance (VRT), and response 
  time mean (MRT), to EZ diffusion model parameters.
}
\usage{
Data2EZ(Pc, VRT, MRT, s = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pc}{ Proportion of correct responses.}
  \item{VRT}{ Variance of the response times in seconds.}
  \item{MRT}{ Mean of the response times in seconds.}
  \item{s}{ Scale parameter for the parameters (Ratcliff's convention is \code{s = 0.1}, the default).}
}
\details{
  Use of MRT is optional; if MRT is absent \code{Ter} will not be computed.
}
\value{
  A list with estimates
  \item{v }{Drift rate}
  \item{a }{Boundary separation}
  \item{Ter}{Non-decision time (is NA if MRT is absent)}
}
\references{ 
Grasman, R.P.P.P., Wagenmakers, E.-J.& van der Maas, H.L.J. (2009). On the mean and variance of response times under the diffusion model with an application to parameter estimation. Journal of Mathematical Psychology, 53(2), 55-68. 
}
\author{ Raoul Grasman }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{EZ2batch}} }
\examples{
    Data2EZ(Pc=0.8022, VRT=0.112035, MRT=0.7231)
}
