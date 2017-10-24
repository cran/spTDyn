\name{initials}
\alias{initials}
\title{Initial values for the spatio-temporal models.}
\description{
This command is useful to assign the initial values of the hyper-parameters of the prior distributions.
}
\usage{
initials(sig2eps=0.01, sig2eta=NULL, sig2beta=NULL, sig2delta=NULL,
   rhotp=NULL, rho=NULL, beta=NULL, phi=NULL)
}
\arguments{
\item{sig2eps}{Initial value for the parameter \eqn{\sigma}^2_\eqn{\epsilon}.}
\item{sig2eta}{Initial value for the parameter \eqn{\sigma}^2_\eqn{\eta}.}
\item{sig2beta}{Initial value for the parameter \eqn{\sigma}^2_\eqn{\beta} for spatially varying model.}
\item{sig2delta}{Initial value for the parameter \eqn{\sigma}^2_\eqn{\delta} for dynamic state-space model.}
\item{rhotp}{Value for the parameter \eqn{\rho} for dynamic state-space model. For \code{rhotp=1}, \eqn{\rho} parameters are not sampled and fixed at value 1. For \code{rhotp=0}, \eqn{\rho} parameters are sampled from the full conditional distribution via MCMC with initial value 0.}
\item{rho}{Initial value for the parameter \eqn{\rho}.}
\item{beta}{Initial value for the parameter \eqn{\beta}.}
\item{phi}{Initial value for the parameter \eqn{\phi}.}
}
\note{
Initial values are automatically given if the user does not provide these. 
}
\seealso{
\code{\link{GibbsDyn}, \link{priors}}.
}
\examples{
%\dontrun{
## 

initials<-initials(sig2eps=0.01, sig2eta=0.5, beta=NULL, phi=0.001)
initials

##
%}
}
\keyword{spTDyn}