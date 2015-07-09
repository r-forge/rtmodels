\name{EZ2.pe}
\alias{EZ2.pe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ compute probability of exit through lower bound of a diffusion with constant drift }
\description{
  Computes the probability of exit through the lower bound of a univariate diffusion process with constant drift on an interval with absorbing boundaries. Used as a model of information accumulation, it is gives the probability of an error response in a speeded two-alternative forced choice (2AFC) response time task.
}
\usage{
EZ2.pe(nu, z, a, s = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu}{ Drift rate. }
  \item{z}{ Starting point. }
  \item{a}{ Boundary separation }
  \item{s}{ Scaling parameter (Ratcliff's convention is \code{s = 0.1}, the default) }
}
\details{
  This process as a model of information accumulation and decision is Ratcliff's diffusion model (1978). It can be used e.g., to compute the mean response times of the correct responses in a lexical decision time, given the drift rate, the bias (start point), and criterion (boundary separation).
}
\value{
%  ~Describe the value returned
    \code{EZ2.pe} returns the exit probability through lower end of the interval (0,a)
    The return value has the attribute "gradient" attached: the gradient with respect to each of the parameters.
}
\references{ 
Ratcliff, R. (1978). Theory of Memory Retrieval. \emph{Psychological review} vol. 85 (2) pp. 59-108 

Grasman, R. P. P., Wagenmakers, E.-J., & van der Maas, H. L. J. (2007). On the mean and variance of response times under the diffusion model with an application to parameter estimation, \emph{J. Math. Psych.} 53: 55--68.
}
\author{ Raoul P. P. P. Grasman }
\seealso{ \code{\link{EZ2-package}}, \code{\link{EZ2.cmrt}}, \code{\link{EZ2.cvrt}}, \code{\link{EZ2.mrt}}, \code{\link{EZ2.vrt}} }
\examples{
    EZ2.pe(.1, .08, .12)
}
