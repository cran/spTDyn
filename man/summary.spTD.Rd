\name{summary.spTD}
\alias{summary.spTD}
\title{Summary statistics of the parameters.}
\description{
This function is used to obtain MCMC summary statistics.
}
\usage{
\method{summary}{spTD}(object, digits=4, package="spTDyn", coefficient=NULL, ...)

## 
}
\arguments{
\item{object}{Object of class inheriting from "spTD".}
\item{digits}{Rounds the specified number of decimal places (default 4).}
\item{package}{If "coda" then summary statistics are given using coda package. Defaults value is "spTDyn".}
\item{coefficient}{Takes values: "spatial", "temporal" and "rho" for summary statistics of spatial, temporal and rho coefficients respectively. If NULL then provides parameter summary without spatial and temporal coefficients.}
\item{...}{Other arguments.}
}
\value{
\item{sig2eps}{Summary statistics for \eqn{\sigma_\epsilon^2}.} 
\item{sig2eta}{Summary statistics for \eqn{\sigma_\eta^2}.} 
\item{phi}{Summary statistics for spatial decay parameter \eqn{\phi}, if estimated using \code{decay}.} 
\item{...}{Summary statistics for other parameters used in the models.} 
}
\seealso{
\code{\link{GibbsDyn}}.
}
\examples{
\dontrun{
##

summary(out) # where out is the output from spT class
summary(out, digit=2) # where out is the output from spT class
summary(out, pack="coda") # where out is the output from spT class
summary(out, coefficient="spatial") # for spatially varying coefficients
summary(out, coefficient="temporal") # for temporally varying coefficients

##
}
}
\keyword{spT}
