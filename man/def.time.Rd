\name{def.time}
\alias{def.time}
\title{Timer series information.}
\description{
This function defines the time series in the spatio-temporal data.
}
\usage{
def.time(t.series, segments=1)
}
\arguments{
\item{t.series}{Number of times within each segment in each series. Can take only regular time-series. }
\item{segments}{Number of segments in each time series. This should  be a constant. }
}
\seealso{
\code{\link{GibbsDyn}}.
}
\examples{
\dontrun{
## 

# regular time-series in each year
time.data<-def.time(t.series=30,segments=2)

##
}
}
\keyword{spT}