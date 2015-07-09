\name{drt.error}
\alias{drt.error}
\alias{prt.error}
\alias{p.error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ RT distribution of error responses}
\description{
  Density and distribution functionof the reaction times of the error responses, and probability of an error, in Ratcliff's diffusion model for two alternative choice RT's.
}
\usage{
drt.error(t, v, a, z, s = 0.1, k = 1:100)
prt.error(t, v, a, z, s = 0.1, k = 1:100)
p.error(v, a, z = a/2, s2 = 0.1^2)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ reaction time }
  \item{v}{ drift rate }
  \item{a}{ boundary separation }
  \item{z}{ starting point of the diffusion }
  \item{s}{ diffusion constant }
  \item{s2}{ square of diffusion constant }
  \item{k}{ terms evaluated in the infinite series (the more, the more accurate) }
}
\details{
  This function is vectorized in \code{t}. For small values of \code{t} you need a large range for \code{k}.
}
\value{
  Density of the error responses.
}
\references{ Ratcliff & Tuerlinckx (2002) Estimating parameters of the diffusion model: Approaches to dealing with contaminant reaction times and parameter variability. Psychonomic Bulletin & Review, 9(3), 438--481. }
\author{ Raoul Grasman }
\section{Warning }{Beta version}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (t, v, a, z, s = 0.1, k = 1:100)
{
    defect = sapply(t, function(t) pi * s^2 * (exp(-v * z/s^2)/a^2) *
        sum(k * sin(k * pi * z/a) * exp(-(v^2/s^2 + pi^2 * k^2 *
            s^2/a^2) * t/2)))
    defect/p.error(v = v, a = a, z = z, s2 = s^2)
  }
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
