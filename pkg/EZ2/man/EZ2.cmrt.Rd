\name{EZ2.cmrt}
\alias{EZ2.cmrt}
\alias{EZ2.cvrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ compute mean exit (decision) time conditional on exit point (chosen alternative) of a diffusion process}
\description{
  Given a boundary separation, a starting point, and a drift rate, this function computes the mean exit time/exit time variance of a one dimensional diffusion process under constant drift on an interval with absorbing boundaries, conditioned on the point of exit. Used as a model of information accumulation, it is gives the mean decision time/decision time variance of responses in a speeded two-alternative forced choice (2AFC) response time task, conditional on what alternative was decided upon.  
  %~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EZ2.cmrt(nu, z, a, s = 0.1)
EZ2.cvrt(nu, z, a, s = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{ Drift rate. }
  \item{z}{ Starting point. }
  \item{a}{ Boundary separation }
  \item{s}{ Scaling parameter (Ratcliff's convention is \code{s = 0.1}, the default) }
}
\details{
%  ~~ If necessary, more details than the description above ~~
  This model of information accumulation and decision is a simplified version of Ratcliff's diffusion model (1978). It can be used e.g., to compute the mean response times of the correct responses in a lexical decision time, given the drift rate, the bias (start point), and criterion (boundary separation).
}
\value{
%  ~Describe the value returned
    \code{EZ2.cmrt} returns the mean exit/decision time(s)
    \code{EZ2.cvrt} returns the exit/decision time variance(s)
    The return value has the attribute "gradient" attached: the gradient with respect to each of the parameters.
}
\references{ 
    Ratcliff. Theory of Memory Retrieval. Psychological review (1978) vol. 85 (2) pp. 59-108 
    Grasman et al. EZ2: An extension of teh EZ-diffusion model for Response Time and Accuracy. Journal of Mathematical Psychology (2007) submitted
     }
\author{ Raoul Grasman }
\note{ %~~further notes~~ 

 %~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ \code{\link{EZ2-package}}, \code{\link{EZ2.mrt}}, \code{\link{EZ2.vrt}} }
\examples{
    EZ2.mrt( 0.1, 0.08, 0.12, Ter=0.25)
    EZ2.vrt( 0.1, 0.08, 0.12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
