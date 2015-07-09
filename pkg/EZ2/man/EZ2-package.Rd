% Date: 2007-08-31
\name{EZ2-package}
\alias{EZ2-package}
\docType{package}
\title{
   EZ2 diffusion modeling of Response Time and Accuracy
}
\description{
   EZ2 fits a simplified version of the diffusion model for human and monkey response times 
   and accuracy to the means and variances of the observables. Ratcliff's diffusion model specifies
   the information accumulation process in the brain in a two-alternative forced choice perceptual
   decision as a drift diffusion that conforms to the stochastic differential equation
   \deqn{dZ(t) = \nu dt + sdW(t),}
   where \eqn{\nu} specifies the accumulation speed (which is a function of perceptual quality of the
   stimulus as well as the quality of the stimulus' stored representation in memory), and \eqn{s} is
   the noisiness of the accumulation process. The process is assumed to start at \eqn{Z(0) = z} (\eqn{z}
   is often called the \dQuote{starting point}). Each alternative is associated with a direction in
   which the process \eqn{Z(t)} develops. A decision can be made only if sufficient information
   is accumulated towards either of the alternatives. That is, the model assumes that there are
   \emph{thresholds} associated with each of the two alternatives, which cause the brain to
   make the decision in favor of the reached threshold. For the alternative that is associated 
   with the downward direction, the threshold is conventionally set at 0. For the alternative that 
   is associated with the upward direction the threshold is conventionally called \eqn{a} (the 
   \dQuote{boundary separation}). The starting point \eqn{z} lies in between the lower and
   upper threshold (\eqn{0 < z < a}) and is interpreted as a bias parameter: If one alternative
   is for some reason favored over the other (e.g., because of relative frequency of occurence),
   the starting point is closer to that alternative. The time it takes \eqn{Z} to reach either
   of these thresholds is considered the \emph{Decision time}. (In the Stochastic Calculus 
   literature this is also called the \emph{exit time} of the process \eqn{Z} starting at \emph{Z(0)=z}
   from the interval \emph{(0,a)}.) The Response Time for a given
   decision is the sum of the Decision time, stimulus encoding time, and motor execution time.
   The latter two times are lumped into a \emph{Non-decision time} parameter often called \eqn{T_er}.
   The upper threshold alternative is most commonly identified with the correct alternative 
   for a given stimulus, while the lower threshold alternative is identified with the incorrect 
   alternative. Hence, there are \sQuote{correct} response times which have a mean and variance, 
   and there are \sQuote{incorrect} response times which also have a mean and variance. These
   means and variances are predicted by the diffusion model. Furthermore, the diffusion model
   predicts the percentage of error-responses.
   
   In Ratcliff's diffusion model, holding all things equal across different trials, the starting
   point \eqn{z}, the information speed \eqn{\nu}, and the non-decision time \eqn{T_{er}} still
   have some variability. This variability is captured in extra parameters that specify the 
   across trials variances of these parameters. The EZ and EZ2 models assume that, holding all
   things equal across different trials, all parameters are constant. Hence, EZ and EZ2 leave
   out these extra variance parameters. The difference between EZ and EZ2 is that the EZ
   model assumes that \eqn{z = a / 2}, which is appropriate in many experiments, while EZ2
   does not make this assumption.
   
   This package contains functions to compute the response time means and variances of the
   decision (exit) times and proportions of error, given the parameter values of the EZ(2) model.
   The means and variances can be computed either conditionally on wether the response was 
   (in)correct or by pooling across correct and incorrect response. Furthermore, the package
   provides utility functions for fitting data. The original intent was to provide method
   of moment estimators, but the package also supports least squares fitting of extensive 
   models (including models that allow for contaminant RTs) to RTs from multiple experimental
   conditions.
}
\details{
\tabular{ll}{
Package: \tab EZ2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-08-31\cr
License: \tab GPL version 2 (or later)\cr
} need a data frame containing (at least) the moments (RT means, RT variances, and proportions correct) that you would like to model.

Specify for each used moment a formula of the form \code{vrt1 ~ EZ2.vrt(v1, z, a)}.

Supply the model with the data and a starting point for all unknown parameters to \code{EZ2} 
or \code{EZ2batch} to calculate method of moment estimators or least squares estimators of 
the unknown parameters in the model (that is, \code{v1}, \code{z}, and \code{a} in 
\code{vrt1 ~ EZ2.vrt(v1, z, a)}).

Your model can be complex and extensive; for instance 

\code{vrt1 ~ EZ2.vrt(v1, z, a) + p0*(maxRT-minRT)^2/12+p0*(1-p0)*(EZ2.mrt(v1, z, a)-(maxRT+minRT)/2)^2}

is equaly valid.
}

\author{
Raoul P. P. P. Grasman
Maintainer: Raoul Grasman <rgrasman@uva.nl>
}
\references{
Grasman, R. P. P., Wagenmakers, E.-J., & van der Maas, H. L. J. (2007). On the mean and variance of response times under the diffusion model with an application to parameter estimation, \emph{J. Math. Psych.} 53: 55--68.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
%~~ simple examples of the most important functions ~~
## create some data (theoretical values, not simulated) Needless to say, in reality 
## you would like to fit real data!
A = seq(.08,.13,len=6)
X2 = data.frame(A=A)
X2$vrt0 = sapply(A, function(a) EZ2.vrt(.1,.05,a))
X2$pe0 = sapply(A, function(a) EZ2.pe(.1,.05,a))
X2$vrt1 = sapply(A, function(a) EZ2.vrt(.2,a-.05,a))
X2$pe1 = sapply(A, function(a) EZ2.pe(.2,a-.05,a))

X2 = as.data.frame(X2)          # now pretend that X2 is the data frame that 
                                # you may have computed from real data

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
%\keyword{ models }
%\keyword{ package }