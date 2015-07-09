\name{cdfdif}
\alias{cdfdif}
\alias{CDFDif}
\title{ Compute simultaneous probability of (a) response time(s) and error under Ratcliff's diffusion model. }
\description{
\code{CDFDif} Returns the probability Prob(RT <= t, X=x), where RT is a single response time, and X is an indicator
(correct or incorrect).
\code{cdfdif} does the same, except that t is allowed to be an array.
}
\usage{
CDFDif(t,x,nu,z,a,Ter,eta=0,sz=0,st=0);
cdfdif(t,x,nu,z,a,Ter,eta=0,sz=0,st=0)
}
\arguments{
\item{t}{ degrees of freedom parameter }
\item{x}{ TRUE (=correct) or FALSE (=incorrect) }
\item{a}{ the boundary separation }
\item{Ter}{ the mean non-decisional component time }
\item{eta}{ the standard deviation of the normal drift distribution }
\item{z}{ the mean starting point }
\item{sz}{ the spread of the starting point distribution }
\item{st}{ the spread of the non-decisional component time distribution }
\item{nu}{ the mean drift rate }
}
\value{
The probability/probabilities Prob(RT <= t, X=x).
}
\references{ Tuerlinckx, (2004), The efficient computation of the cumulative distribution and probability density functions in the diffusion model, \emph{Behavior Research Methods, Instruments, & Computers}, vol. 36 (4): 702-716 }
\author{ Raoul Grasman \email{rgrasman@uva.nl}}
\section{Remark}{Original Fortran sources due to Francis Tuerlinckx. C translation by Raoul Grasman.}
\examples{
##
cdfdif(c(0.3,0.5),1,nu=.2,z=.11,a=.2,Ter=.3,eta=.1,sz=0.05,st=0.1)
}
\keyword{ distribution }
