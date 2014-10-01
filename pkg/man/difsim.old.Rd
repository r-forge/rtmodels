\name{difsim.old}
\alias{difsim.old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to simulate diffusion processes model of Rattcliff for perceptual decisions }
\description{
  Calling difsim results in the reaction times of both correct and incorrect responses of \code{N} trials, as obtained from a diffusion process. All parameters of the diffusion process can be modified. The underlying process is given by the stochastic differential equation \eqn{dx = v dt + \sqrt{var.v} * dW(t)}. }
\usage{
difsim.old(v, a, N, z = 0.5 * a, s = 0.1, var.v = 0, range.z = 0, Ter = 0, range.Ter = 0, stepsize = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{ Drift rate of the diffusion }
  \item{a}{ Upper decision boundary. Lower decision boundary is 0. }
  \item{N}{ The number of trials. }
  \item{z}{ Starting point of the difusion, defaults to the center of the lower and upper boundaries.~~ }
  \item{s}{ Variance of \eqn{dW(t)} }
  \item{var.v}{ Variance of the driftrate. Defaults to zero. }
  \item{range.z}{ Radius of uniform distributed starting point around z. }
  \item{Ter}{ Fixed constant for motor respons. }
  \item{range.Ter}{ Radius of uniform distribution around Ter. }
  \item{stepsize}{ Size (msec) of the steps taken in the random walk. Default value is appropriate voor common parameter values. }
}
\details{
  For details on the meaning of the parameters, see the references.
}
\value{
  \item{t.top }{reaction times in trials hitting the boundary \code{a} ("correct" responses).}
  \item{t.bot }{reaction times in trials hitting the boundary \code{0} ("incorrect" responses).}
}
\references{ Ratcliff & Tuerlinckx, (2002), Estimating parameters of the diffusion model: Approaches to dealing with contaminant reaction times and parameter variability, \emph{Psychonomic Bulletin & Review}, vol. 9 (3): 438-481 }
\author{
Raoul Grasman
Eric-Jan Wagenmakers
}
\section{Warning}{This function is obsolete.}
\examples{
  rt = difsim(v=.2,a=.2,z=.09,Ter=.250,N=1000);
  quantile(rt[[1]],seq(0.05,0.95,l=10))
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
\keyword{ distribution }
\keyword{ models }
