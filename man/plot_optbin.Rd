\name{plot.optbin}
\alias{plot.optbin}
\title{
Plotting Optimal Bins
}
\description{
\code{plot} method for class \code{optbin}.
}
\usage{
\method{plot}{optbin}(x, col=NULL, main="Binned Observations", \dots)
}
\arguments{
\item{x}{
an object of class \code{optbin}.
}
\item{col}{
vector of colors to apply to bins (default uses an internal set)
}
\item{main}{
title of graph
}
\item{\dots}{
other parameters passed through to the underlying plotting routines (do not
set xaxt or ann)
}
}
\details{
The plot will contain the sorted points of the data that generated the
bins.  Points are color-coded per bin, and the plot contains the average
value over the bin as a line.  x axis labels are the upper thresholds for
each bin.
}
\value{
None
}
\seealso{
\code{\link{optbin}}
}
\keyword{optbin}
