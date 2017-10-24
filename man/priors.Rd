\name{priors}
\alias{priors}
\title{Priors for the spatio-temporal models.}
\description{This command is useful to assign the hyper-parameters of the prior distributions.
}
\usage{
priors(inv.var.prior=Gamm(a=2,b=1),beta.prior=Norm(0,10^10), 
  rho.prior=Norm(0,10^10))
}
\arguments{
\item{inv.var.prior}{The hyper-parameter for the Gamma prior distribution (with mean = a/b) of the precision (inverse variance) model parameters (e.g., 1/\eqn{\sigma}2_\eqn{\epsilon}, 1/\eqn{\sigma}2_\eqn{\eta}).}
\item{beta.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\beta} model parameters.}
\item{rho.prior}{The hyper-parameter for the Normal prior distribution of the \eqn{\rho} model parameter.}
}
\note{
If no prior information are given (assigned as NULL), then it use flat prior values of the corresponding distributions. \cr
\code{Gam} and \code{Nor} refers to Gamma and Normal distributions respectively.
}
\seealso{
\code{\link{GibbsDyn}, \link{initials}}.
}
\examples{
%\dontrun{
## 
library(spTimer)
priors<-priors(inv.var.prior=Gamm(2,1), beta.prior=Norm(0,10^4))
priors

##
%}
}
\keyword{spTDyn}
