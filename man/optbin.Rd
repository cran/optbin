\name{optbin}
\alias{optbin}
\title{
Optimal Binning of Continuous Variables
}
\description{
Determines break points in numeric data that minimize the difference between
each point in a bin and the average over it.
}
\usage{
optbin(x, numbin, metric=c('se', 'mse'), is.sorted=FALSE, max.cache=2^31, na.rm=FALSE)
}
\arguments{
\item{x}{
numeric data
}
\item{numbin}{
number of bins to partition vector into
}
\item{metric}{
minimize squared error (se) between values and average over bin, or mean
squared error (mse) dividing squared error by bin length
}
\item{is.sorted}{
set true if x is already in increasing order
}
\item{max.cache}{
maximum memory in bytes to use to cache bin metrics; if analysis would need
more than use slower calculation without cache
}
\item{na.rm}{
drop NA values (which may occur when converting the data to a vector),
otherwise cannot proceed with binning
}
}
\details{
Data is converted into a numeric vector and sorted if necessary.  Internally
bins are determined by positions within the vector, with the breaks inclusive
at the upper end.  The bin thresholds are the same, so bin b covers the range
thr[b-1] < x <= thr[b], where thr[0] is -Inf.  The routine finds the first
split found with the best metric, if there is more than one.

The library uses an exhaustive search over all possible breakpoints.  It
begins by finding the best splits with 2 bins for all pairs of start and
endpoints, then adds a third bin, and so on.  This rejects most alternatives
at each level, leaving an O(nbin * nval * nval) algorithm.
}
\value{
An object of class 'optbin' with components:
\item{x}{
the original data, sorted
}
\item{numbins}{
the number of bins created
}
\item{call}{
argument values when function called
}
\item{metric}{
cost function used to select best partition
}
\item{minse}{
value of SE/MSE metric for all bins
}
\item{thr}{
upper threshold of bin range, inclusive
}
\item{binavg}{
average of values in each bin
}
\item{binse}{
value of SE/MSE metric for each bin
}
\item{breaks}{
positions of endpoint (inclusive) of each bin in x
}
}
\seealso{
\code{\link{assign.optbin}}, \code{\link{print.optbin}},
\code{\link{summary.optbin}}, \code{\link{plot.optbin}}
}
\examples{
## Well separated groups
set.seed(17)
d1 <- c(rnorm(75, mean=1, sd=0.2), rnorm(75, mean=3, sd=0.2),
        rnorm(84, mean=6, sd=0.2), rnorm(75, mean=9, sd=0.2),
        rnorm(75, mean=11, sd=0.2), rnorm(150, mean=15, sd=0.2))
## Divides into groups 1+2+3, 4+5, 6, metric is 1176.3
binned3 <- optbin(d1, 3)
summary(binned3)
plot(binned3)
## Divides into groups 1, 2, 3, 4+5, and 6, metric is 169.9
binned5 <- optbin(d1, 5)
plot(binned5)
## Divides into separate groups, metric is 24.4
binned6 <- optbin(d1, 6)
summary(binned6)
plot(binned6)
## Each rnorm group divides roughly in half.
binned12 <- optbin(d1, 12)
plot(binned12)
## A grouping that overlaps, bins near but not at minima between peaks
d2 <- c(rnorm(300, mean=1, sd=0.25), rnorm(400, mean=2, sd=0.25),
        rnorm(300, mean=3, sd=0.25))
binned3b <- optbin(d2, 3)
hist(binned3b, breaks=50, col='yellow')
}
\keyword{optbin}

