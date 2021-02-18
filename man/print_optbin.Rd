\name{print.optbin}
\alias{print.optbin}
\title{
Printing Optimal Bins
}
\description{
\code{print} method for class \code{optbin}.
}
\usage{
\method{print}{optbin}(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{optbin}.
}
\item{\dots}{
generic arguments (ignored)
}
}
\details{
Shows the upper bounds of each bin, ie. bin b covers
threshold[b-1] < x <= threshold[b] where threshold[0] is -Inf.  Also prints
the total (mean) squared error sum over all bins.
}
\value{
The argument x unchanged, an object of class 'optbin' with components:
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
\code{\link{optbin}}, \code{\link{summary.optbin}}
}
\keyword{optbin}

