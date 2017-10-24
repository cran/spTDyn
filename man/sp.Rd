\name{sp}
\alias{sp}
\title{Defining spatially varying coefficients in the formula}
\description{
This function is used to define spatially varying coefficients within the formula for the Gaussian process spatio-dynamic and spatially varying coefficient process models. 
}
\usage{
sp(x)
}
\arguments{
\item{x}{The variable/covariate for which spatially varying coefficient is defined.}
}
\seealso{
\code{\link{GibbsDyn}, \link{tp}}
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
post.gp <- GibbsDyn(formula=o8hrmax ~cMAXTMP+WDSP+sp(RH),   
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
coords<-as.matrix(unique(cbind(DataFit[,2:3])))

# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- GibbsDyn(formula=o8hrmax ~cMAXTMP+WDSP+sp(RH),   
         data=NYdata, coords=coords, scale.transform="SQRT")
print(post.gp)

}
}
\keyword{spT}