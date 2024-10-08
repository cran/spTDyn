\name{predict.spTD}
\alias{predict.spTD}
\title{Spatial and temporal predictions for the spatio-temporal models.}
\description{
This function is used to obtain spatial predictions in the unknown locations and also to get the temporal forecasts using MCMC samples.
}
\usage{
\method{predict}{spTD}(object, newdata, newcoords, foreStep=NULL, type="spatial", 
        nBurn, tol.dist, Summary=TRUE, ...)

}

\arguments{
\item{object}{Object of class inheriting from "spTD".}
\item{newdata}{The data set providing the covariate values for spatial prediction or temporal forecasts. This data should have the same space-time structure as the original data frame.}
\item{newcoords}{The coordinates for the prediction or forecast sites. The locations are in similar format to \code{coords}, see \code{\link{spT.Gibbs}}.}
\item{foreStep}{Number of K-step (time points) ahead forecast, K=1,2, ...; Only applicable if type="temporal".}
\item{type}{If the value is "spatial" then only spatial prediction will be performed at the \code{newcoords} which must be different from the fitted sites provided by \code{coords}.  When the "temporal" option is specified then forecasting will be performed  and in this case the \code{newcoords} may also contain elements of the fitted sites in which case only temporal forecasting beyond the last fitted time point will be performed. } 
\item{nBurn}{Number of burn-in. Initial MCMC samples to discard before making inference.}
\item{tol.dist}{Minimum tolerance distance limit between fitted and predicted locations.}
\item{Summary}{To obtain summary statistics for the posterior predicted MCMC samples. Default is TRUE.}
\item{...}{Other arguments.}
}
\value{
\item{pred.samples or fore.samples}{Prediction or forecast MCMC samples.}
\item{pred.coords or fore.coords}{prediction or forecast coordinates.}
\item{Mean}{Average of the MCMC predictions}
\item{Median}{Median of the MCMC predictions}
\item{SD}{Standard deviation of the MCMC predictions}
\item{Low}{Lower limit for the 95 percent CI of the MCMC predictions}
\item{Up}{Upper limit for the 95 percent CI of the MCMC predictions}
\item{computation.time}{The computation time.}
\item{model}{The model method used for prediction.}
\item{type}{"spatial" or "temporal".}
\item{...}{Other values "obsData", "fittedData" and "residuals" are provided only for temporal prediction.}
}
\references{

 Bakar, K. S., Kokic, P. and Jin, H. (2015). A spatio-dynamic model for assessing frost risk in south-eastern Australia. Journal of the Royal Statistical Society, Series C.
 Bakar, K. S., Kokic, P. and Jin, H. (2015). Hierarchical spatially varying coefficient and temporal dynamic process models using spTDyn. Journal of Statistical Computation and Simulation.

}
\seealso{
\code{\link{GibbsDyn}}.
}
\examples{
\dontshow{
##
## Read Aus data ##
data(AUSdata)
# set a side data for validation
library(spTimer)
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

nItr<-100
nBurn<-0

# MCMC via Gibbs using defaults
# Spatially varying coefficient process model

set.seed(11)
post.sp <- GibbsDyn(tmax ~ soi+sp(soi)+grid+sp(grid),
           data=AUSdataFit, nItr=nItr, nBurn=nBurn, coords=~lon+lat,
           spatial.decay=decay(distribution=Gamm(2,1),tuning=0.06))
## spatial prediction
set.seed(11)
pred.sp <- predict(post.sp,newcoords=~lon+lat,newdata=AUSdataPred)
## temporal prediction
set.seed(11)
pred.sp.f <- predict(post.sp,type="temporal",foreStep=12,
                     newcoords=~lon+lat, newdata=AUSdataFore)

# MCMC via Gibbs using defaults
# spatio-temporal DLM

set.seed(11)
post.tp <- GibbsDyn(o8hrmax ~ tp(cMAXTMP)-1, data=NYdataFit,
           nItr=nItr, nBurn=nBurn, coords=~Longitude+Latitude,
           initials=initials(rhotp=0), scale.transform="SQRT",
           spatial.decay=decay(distribution=Gamm(2,1),tuning=0.05))
## spatial prediction
set.seed(11)
pred.tp <- predict(post.tp, newdata=NYdataPred, newcoords=~Longitude+Latitude)
## temporal prediction
set.seed(11)
pred.tp.f <- predict(post.tp, newdata=NYdataFore, newcoords=~Longitude+Latitude,
                     type="temporal", foreStep=7)
##

}
\donttest{
##

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
\keyword{spT}