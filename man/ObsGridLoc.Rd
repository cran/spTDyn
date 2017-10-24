\name{ObsGridLoc}
\alias{ObsGridLoc}
\alias{gridTodata}
\alias{ObsGridData}
\title{Combining observation and nearest grid locations and data.}
\description{
These commands combine observation and nearest grid locations, data.}
\usage{
ObsGridLoc(obsLoc, gridLoc, distance.method="geodetic:km", plot=FALSE)
gridTodata(gridData, gridLoc=NULL, gridLon=NULL, gridLat=NULL)
ObsGridData(obsData, gridData, obsLoc, gridLoc, distance.method="geodetic:km")
}
\arguments{
\item{obsLoc}{The observed/measurement locations, first column is longitude/easting/x-axis and second column is latitude/northing/y-axis.}
\item{gridLoc}{Grid locations, first column is longitude/easting/x-axis and second column is latitude/northing/y-axis.}
\item{distance.method}{The preferred method to calculate the distance between any two locations. The available options are "geodetic:km", "geodetic:mile", "euclidean", "maximum", "manhattan", and "canberra". See details in \code{\link{dist}}.}
\item{plot}{Logical argument, if TRUE then plot observed and nearest grid locations.}
\item{gridData}{Gridded data, should be in array form with dimenstions as longitude/x-axis, latitude/y-axis, day/time1, year/time2.}
\item{gridLon}{Longitude/easting/x-axis of grid locations.}
\item{gridLat}{Latitude/northing/y-axis of grid locations.}
\item{obsData}{Observation data in data frame.}
}
\examples{
%\dontrun{
##

library(spTimer)
data(NYdata)	
data(NYgrid)

obsLoc<-unique(cbind(NYdata$Longitude,NYdata$Latitude))
gridLoc<-unique(cbind(NYgrid$Longitude,NYgrid$Latitude))

# find closest observed and grid locations
dat<-ObsGridLoc(obsLoc, gridLoc)
head(dat)
# with plots
dat<-ObsGridLoc(obsLoc, gridLoc, plot=TRUE)
head(dat)

# convert array gridData to spTimer data format
gridData<-array(1:(10*10*31*2),dim=c(10,10,31,2)) # lon, lat, day, year
dat<-gridTodata(gridData, gridLoc)
head(dat)

# combine observed and grid data and locations
obsData<-NYdata
gridData<-array(1:(10*10*31*2),dim=c(10,10,31,2)) # lon, lat, day, year
dat<-ObsGridData(obsData, gridData, obsLoc, gridLoc)
head(dat)

# combine observed and more than one grid datasets
obsData<-NYdata
gridData1<-array(1:(10*10*31*2),dim=c(10,10,31,2)) # lon, lat, day, year
gridData2<-array(((10*10*31*2)+1):(2*(10*10*31*2)),dim=c(10,10,31,2)) # lon, lat, day, year
gridLoc1<-unique(cbind(NYgrid$Longitude,NYgrid$Latitude))
gridLoc2<-unique(cbind(NYgrid$Longitude,NYgrid$Latitude))
dat<-ObsGridData(obsData, gridData=list(gridData1,gridData2),
        obsLoc, gridLoc=list(gridLoc1, gridLoc2))
head(dat)

##
%}
}
\keyword{utility}