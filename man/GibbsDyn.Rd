\name{GibbsDyn}
\alias{GibbsDyn}
\title{MCMC sampling for the models.}
\description{
This function is used to draw MCMC samples using the Gibbs sampler. 
}
\usage{
GibbsDyn(formula, data=parent.frame(), time.data=NULL, coords, 
	priors=NULL, initials=NULL, nItr=5000, nBurn=1000, report=1, tol.dist=0.05, 
	distance.method="geodetic:km", cov.fnc="exponential", scale.transform="NONE", 
	spatial.decay=decay(distribution="FIXED"))
}
\arguments{
\item{formula}{The symnbolic description of the model equation of the regression part of the space-time model. The terms \link{sp} and \link{tp} are used to define spatially and temporally varying parameters for the model.} 
\item{data}{An optional data frame containing the variables in the model. If omitted, the variables are taken from environment(formula), typically the environment from which spT.Gibbs is called. The data should be ordered first by the time and then by the sites specified by the \code{coords} below. One can also supply coordinates through this argument, where coordinate names should be "Latitude" and "Longitude".}
\item{time.data}{Defining the segments of the time-series set up using the function \code{\link{def.time}}.}
\item{coords}{The n by 2 matrix or data frame defining the locations (e.g., longitude/easting, latitude/northing) of the fitting sites, where n is the number of fitting sites. One can also supply coordinates through a formula argument such as ~Longitude+Latitude.}
\item{priors}{The prior distributions for the parameters. Default distributions are specified if these are not provided. If priors=NULL a flat prior distribution will be used with large variance. See details in \code{\link{priors}}.}
\item{initials}{The preferred initial values for the parameters. If omitted, default values are provided automatically. Further  details are provided in  \code{\link{initials}}.}
\item{nItr}{Number of MCMC iterations. Default value is 5000.}
\item{nBurn}{Number of burn-in samples. This number of samples will be discarded before making any inference. Default value is 1000.}
\item{report}{Number of reports to display while running the Gibbs sampler. Defaults to number of iterations. }
\item{distance.method}{The preferred method to calculate the distance between any two locations. The available options are "geodetic:km", "geodetic:mile", "euclidean", "maximum", "manhattan", and "canberra". See details in \code{\link{dist}}. The default is "geodetic:km".}
\item{tol.dist}{Minimum separation distance between any two locations out of those  specified by coords, knots.coords and pred.coords. The default is 0.005. The programme will exit if the minimum distance is less than the non-zero specified value. This will ensure non-singularity of the covariance matrices.}
\item{cov.fnc}{Covariance function for the spatial effects. The available options  are "exponential", "gaussian", "spherical" and "matern". If "matern" is used then by default the smooth parameter (\eqn{\nu}) is estimated from (0,1) uniform distribution using discrete samples.}
\item{scale.transform}{The transformation method for the response variable. Currently implemented options are: "NONE", "SQRT", and "LOG" with  "NONE" as the deault.}
\item{spatial.decay}{Provides the prior distribution for the spatial decay parameter \eqn{\phi}. Currently implemented options  are  "FIXED", "Unif", or "Gamm".  Further details for each of these are specified by \code{\link{decay}}.}
}
\value{
\item{accept}{The acceptance rate for the \eqn{\phi} parameter if the "MH" method of sampling is chosen.}
\item{phip}{MCMC samples for the parameter \eqn{\phi}.}
\item{nup}{MCMC samples for the parameter \eqn{\nu}. Only available if "matern" covariance function is used.}
\item{sig2eps}{MCMC samples for the parameter \eqn{\sigma^2_\epsilon}.}
\item{sig2etap}{MCMC samples for the parameter \eqn{\sigma^2_\eta}.}
\item{sig2betap}{MCMC samples for the parameter \eqn{\sigma^2_\beta}, only applicable for spatially varying coefficient process model.}
\item{sig2deltap}{MCMC samples for the parameter \eqn{\sigma^2_\delta}, only applicable for spatio-temporal DLM.}
\item{sig2op}{MCMC samples for the parameter \eqn{\sigma^2_o}, only applicable for spatio-dynamic and spatio-temporal DLM.}
\item{betap}{MCMC samples for the parameter \eqn{\beta}.}
\item{rhop}{MCMC samples for \eqn{\rho}.}
\item{op}{MCMC samples for the true observations.}
\item{fitted}{MCMC summary (mean and sd) for the fitted values.}
\item{tol.dist}{Minimum tolerance distance limit between the locations.}
\item{distance.method}{Name of the distance calculation method.}
\item{cov.fnc}{Name of the covariance function used in model fitting.}
\item{scale.transform}{Name of the scale.transformation method.}
\item{sampling.sp.decay}{The method of sampling for the spatial decay parameter \eqn{\phi}.}
\item{covariate.names}{Name of the covariates used in the model.}
\item{Distance.matrix}{The distance matrix.}
\item{coords}{The coordinate values.}
\item{n}{Total number of sites.}
\item{r}{Total number of segments in time, e.g., years.}
\item{T}{Total points of time, e.g., days within each year.}
\item{p}{Total number of model coefficients, i.e., \eqn{\beta}'s including the intercept.}
\item{initials}{The initial values used in the model.}
\item{priors}{The prior distributions used in the model.}
\item{PMCC}{The predictive model choice criteria obtained by minimising the expected value of a loss function, see Gelfand and Ghosh (1998). Results for both goodness of fit and penalty are given.}
\item{iterations}{The number of samples for the MCMC chain, without burn-in.}
\item{nBurn}{The number of burn-in period for the MCMC chain.}
\item{computation.time}{The computation time required for the fitted model.}

}
\references{
 Bakar, K. S., Kokic, P. and Jin, H. (2015). A spatio-dynamic model for assessing frost risk in south-eastern Australia. Journal of the Royal Statistical Society, Series C.
 Bakar, K. S., Kokic, P. and Jin, H. (2015). Hierarchical spatially varying coefficient and temporal dynamic process models using spTDyn. Journal of Statistical Computation and Simulation.
}
%\note{
%}
\seealso{
\code{\link{priors}, \link{initials}, \link{dist}, \link{sp}, \link{tp}}.
}
\examples{
\dontrun{
##

###########################
## Attach library spTDyn
###########################

library(spTDyn)

## Read Aus data ##
data(AUSdata)
# set a side data for validation
s<-c(1,4,10)
AUSdataFit<-spT.subset(data=AUSdata, var.name=c("s.index"), s=s, reverse=TRUE)
AUSdataFit<-subset(AUSdataFit, with(AUSdataFit, !(year == 2009)))
AUSdataPred<-spT.subset(data=AUSdata, var.name=c("s.index"), s=s)
AUSdataPred<-subset(AUSdataPred, with(AUSdataPred, !(year == 2009)))
AUSdataFore<-spT.subset(data=AUSdata, var.name=c("s.index"), s=s)
AUSdataFore<-subset(AUSdataFore, with(AUSdataFore, (year == 2009)))

## Read NY data ##
data(NYdata)
# set a side data for validation
s<-c(5,8,10,15,20,22,24,26)
fday<-c(25:31)
NYdataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE)
NYdataFit<-subset(NYdataFit, with(NYdataFit, !(Day \%in\% fday & Month == 8)))
NYdataPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s)
NYdataPred<-subset(NYdataPred, with(NYdataPred, !(Day \%in\% fday & Month == 8)))
NYdataFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=s)
NYdataFore<-subset(NYdataFore, with(NYdataFore, (Day \%in\% fday & Month == 8)))

## Code for analysing temperature data in Section: 4 ##
## Model: Spatially varying coefficient process models ##

nItr<-13000
nBurn<-3000

# MCMC via Gibbs using defaults
# Spatially varying coefficient process model

library("spTDyn", warn.conflicts = FALSE)
set.seed(11)
post.sp <- GibbsDyn(tmax ~ soi+sp(soi)+grid+sp(grid),
           data=AUSdataFit, nItr=nItr, nBurn=nBurn, coords=~lon+lat,
           spatial.decay=decay(distribution=Gamm(2,1),tuning=0.06))
print(post.sp)

## Table: 3, Section: 4.1 ##
post.sp$PMCC

# parameter summary
summary(post.sp) # without spatially varying coefficients
summary(post.sp, coefficient="spatial")

#plot(post.sp, density=FALSE)  # without spatially varying coefficients
#plot(post.sp, coefficient="spatial", density=FALSE)

## Code for Figures: 3(a), 3(b) Section: 4.1 ##
Figure_3a<-function(){
  boxplot(t(post.sp$betasp[1:9,]),pch=".",main="SOI",
          xlab="Sites",ylab="Values")
}
Figure_3b<-function(){
  boxplot(t(post.sp$betasp[10:18,]),pch=".",main="Grid",
          xlab="Sites",ylab="Values")
}
Figure_3a()
Figure_3b()

## spatial prediction
set.seed(11)
pred.sp <- predict(post.sp,newcoords=~lon+lat,newdata=AUSdataPred)

## Table: 4, Section: 4.1, validations ##
spT.validation(AUSdataPred$tmax,c(pred.sp$Mean))
plot(AUSdataPred$tmax,c(pred.sp$Mean))

## temporal prediction
set.seed(11)
pred.sp.f <- predict(post.sp,type="temporal",foreStep=12,
                     newcoords=~lon+lat, newdata=AUSdataFore)

## Table: 4, Section: 4.1, validations ##
spT.validation(AUSdataFore$tmax,c(pred.sp.f$Mean))
plot(AUSdataFore$tmax,c(pred.sp.f$Mean))

## Code for analysing Ozone data in Section: 4 ##
## Model: spatio-temporal DLM ##

# MCMC via Gibbs using defaults
# spatio-temporal DLM

library("spTDyn", warn.conflicts = FALSE)
set.seed(11)
post.tp <- GibbsDyn(o8hrmax ~ tp(cMAXTMP)-1, data=NYdataFit,
           nItr=nItr, nBurn=nBurn, coords=~Longitude+Latitude,
           initials=initials(rhotp=0), scale.transform="SQRT",
           spatial.decay=decay(distribution=Gamm(2,1),tuning=0.05))
print(post.tp)
summary(post.tp)

## Table: 5, Section: 4.2 ##
post.tp$PMCC

## Figure: 5, Section: 4.2 ##
Figure_5<-function(){
  stat<-apply(post.tp$betatp[1:55,],1,quantile,prob=c(0.025,0.5,0.975))
  plot(stat[2,],type="p",lty=3,col=1,ylim=c(min(c(stat)),max(c(stat))),
       pch=19,ylab="",xlab="Days",axes=FALSE,main="cMAXTMP",cex=0.8)
  for(i in 1:55){
    segments(i, stat[2,i], i, stat[3,i])
    segments(i, stat[2,i], i, stat[1,i])
  }
  axis(1,1:55,labels=1:55);axis(2)
  abline(v=31.5,lty=2)
  text(15,0.32,"July");  text(45,0.32,"August");
}
Figure_5()

## spatial prediction
set.seed(11)
pred.tp <- predict(post.tp, newdata=NYdataPred, newcoords=~Longitude+Latitude)

## Table 6, Section: 4.2, validation ##
spT.validation(NYdataPred$o8hrmax,c(pred.tp$Mean))

## temporal prediction
set.seed(11)
pred.tp.f <- predict(post.tp, newdata=NYdataFore, newcoords=~Longitude+Latitude,
                     type="temporal", foreStep=7)

## Table 6, Section: 4.2, validation ##
spT.validation(NYdataFore$o8hrmax,c(pred.tp.f$Mean))

##############################################################################

}
}
\keyword{spTDyn}