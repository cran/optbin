\name{hist.optbin}
\alias{hist.optbin}
\title{
Histogram with Optimal Bins Marked
}
\description{
Draw a histogram of the data used to build the optimal binning and mark
the extent of the bins.
}
\usage{
\method{hist}{optbin}(x, bincol=NULL, main=NULL, xlab=NULL, \dots)
}
\arguments{
\item{x}{
an object of class \code{optbin}.
}
\item{bincol}{
vector of colors for showing extent of bins (default uses an internal set)
}
\item{main}{
plot title, if not specified will modify the normal histogram title
}
\item{xlab}{
x axis label, if not specified will modify the normal histogram label
}
\item{\dots}{
other parameters passed through to hist
}
}
\details{
The points behind the binning are passed unchanged to the histogram function.
Bins are marked with colored bars under the x axis, and lines showing the
average value in each are also drawn on top.
}
\value{
None
}
\seealso{
\code{\link{optbin}},
\code{\link{hist}}
}
\keyword{optbin}
\keyword{histogram}
