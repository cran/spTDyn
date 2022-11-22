\name{tp}
\alias{tp}
\title{Defining dynamic time-series coefficients in the formula}
\description{
This function is used to define dynamic time-series coefficients within the formula for the Gaussian process spatio-dynamic and spatio-temporal DLM. 
}
\usage{
tp(x)
}
\arguments{
\item{x}{The variable/covariate for which time varying coefficient is defined.}
}
\seealso{
\code{\link{GibbsDyn}, \link{sp}}
}
\examples{
\dontshow{
##
library(spTDyn)
# Read data 
data(NYdata); 
# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- GibbsDyn(formula=o8hrmax ~cMAXTMP+WDSP+tp(RH),   
         nItr=100,nBurn=0, 
         data=NYdata, coords=coords, scale.transform="SQRT")
print(post.gp)
##
}
\donttest{
##

###########################
## Attach library spTimer
###########################

library(spTDyn)

###########################
## The GP models:
###########################

##
## Model fitting
##

# Read data 
data(NYdata); 

# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))

# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- GibbsDyn(formula=o8hrmax ~cMAXTMP+WDSP+tp(RH),   
         data=NYdata, coords=coords, scale.transform="SQRT")
print(post.gp)

##
}
}
\keyword{spT}