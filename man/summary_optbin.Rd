\name{summary.optbin}
\alias{summary.optbin}
\title{
Summarizing Optimal Bins
}
\description{
\code{summary} method for class \code{optbin}.
}
\usage{
\method{summary}{optbin}(object, show.range=FALSE, \dots)
}
\arguments{
\item{object}{
an object of class \code{optbin}
}
\item{show.range}{
if true then print the bin's range of points (endpoint inclusive) in the
sorted data
}
\item{\dots}{
generic arguments (ignored)
}
}
\details{
Prints a table with the upper threshold (inclusive), the average of the
data within the bin, and the (mean) squared error sum.  show.range also
adds a column with the start and end indices of the sorted data belonging
to the bin, although this applies to the sorted list and is less useful
in general.
}
\value{
Only called for side-effects (printing).  There is no return value.
}
\seealso{
\code{\link{optbin}}, \code{\link{print.optbin}}
}
\keyword{optbin}

