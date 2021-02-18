\name{assign.optbin}
\alias{assign.optbin}
\title{
Bin Assignment
}
\description{
\code{assign.optbin} returns an object with the same shape as the input
data and values replaced by bin numbers.
}
\usage{
assign.optbin(x, binspec, extend.upper=FALSE, by.value=FALSE)
}
\arguments{
\item{x}{
numeric data to assign
}
\item{binspec}{
an optimal binning partition
}
\item{extend.upper}{
if true then any value in x above the last bin is assigned to that bin,
otherwise its bin is set to NA
}
\item{by.value}{
if true then return average value for bin instead of bin numbers
}
}
\details{
Replaces the values in a copy of the input data by the bin number it
belongs to, or by the bin average value with by.value.  The lowest bin
always extends to -Inf.  The extend.upper argument can open the last bin
to +Inf if true.  Use this function to get in-place bin assignments for
the unsorted data that was passed to optbin.
}
\value{
An object of the same shape as the data.
}
\seealso{
\code{\link{optbin}}
}
\examples{
d <- c(rnorm(30, mean=10, sd=2), rnorm(40, mean=20, sd=2),
       rnorm(30, mean=30, sd=3))
binned <- optbin(d, 3)
assign.optbin(d, binned)
}
\keyword{optbin}
