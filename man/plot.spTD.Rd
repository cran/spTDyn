\name{plot.spTD}
\alias{plot.spTD}
\title{Plots for spTDyn output.}
\description{
This function is used to obtain MCMC summary, residual and fitted surface plots.
}
\usage{
\method{plot}{spTD}(x, residuals=FALSE, coefficient=NULL, ...)

## 
}
\arguments{
\item{x}{Object of class inheriting from "spTD".}
\item{residuals}{If TRUE then plot residual vs. fitted and normal qqplot of the residuals. If FALSE then plot MCMC samples of the parameters using coda package. Defaults value is FALSE.}
\item{coefficient}{Takes values: "spatial", "temporal" and "rho" for summary statistics of spatial, temporal and rho coefficients respectively. If NULL then provides parameter plots without spatial and temporal coefficients.}
\item{...}{Other arguments.}
}
\seealso{
\code{\link{GibbsDyn}}.
}
\examples{
\dontrun{
##

plot(out) # where out is the output from spT class
plot(out, residuals=TRUE) # where out is the output from spT class
plot(out, coefficient="spatial") # for spatially varying coefficients

##
}
}
\keyword{spT}
