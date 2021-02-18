
# optbin:
# Main function that splits data in x into numbin partitions by minimizing
# the (mean) squared error metric over all bins.  Set is.sorted if the data
# is already a vector in order and na.rm to drop any NA values (or NAs
# generated when converting the input to doubles).  max.cache is the maximum
# number of bytes to use before switching to a slower version.

optbin <- function(x, numbin, metric=c('se', 'mse'), is.sorted=FALSE, max.cache=2^31, na.rm=FALSE) {

  if (numbin <= 1) {
    stop("must create more than one bin")
  }
  if (max.cache < 0) {
    max.cache = 0
  }

  x <- as.vector(x, mode='double')
  if (na.rm) {
    x <- x[!is.na(x)]
  } else if (anyNA(x)) {
    stop("data must not include NA")
  }

  metric <- match.arg(metric)
  if (is.na(metric)) {
    stop(sprintf("unknown optimization by %s", metric))
  }

  if (!is.sorted) {
    x <- sort(x)
  }

  res <- .Call('C_optbin', x, numbin, metric=='mse', max.cache, PACKAGE='optbin')
  res$call <- match.call()
  return(res)
}


# assign.optbin:
# Apply the bins defined in binspec to data in x, which need not be sorted.
# The result will have the same shape as x except the values will be the bin
# numbers.  Set extend.upper to put any values above the last bin's threshold
# in that bin, otherwise assign them to NA.  Set by.value to return the
# bin average instead of the bin number.

assign.optbin <- function(x, binspec, extend.upper=FALSE, by.value=FALSE) {

  if (!inherits(binspec, 'optbin')) {
    stop('Bin specification must be from the optbin class')
  }

  binID <- x
  for (b in binspec$numbins:1) {
    if (by.value) {
      binID[x<=binspec$thr[b]] <- binspec$binavg[b]
    } else {
      binID[x<=binspec$thr[b]] <- b
    }
  }

  if (extend.upper) {
    if (by.value) {
      binID[binspec$thr[binspec$numbins]<x] <- binspec$binavg[binspec$numbins]
    } else {
      binID[binspec$thr[binspec$numbins]<x] <- binspec$numbins
    }
  } else {
    binID[binspec$thr[binspec$numbins]<x] <- NA
  }

  binID
}


# print.optbin:
# Show thresholds and (mean) squared error of best partition.

print.optbin <- function(x, ...) {
  if (!inherits(x, 'optbin')) {
    stop('Object must be from the optbin class')
  }

  cat('\n')
  cat('Upper Thresholds (inclusive)\n')
  thr <- x$thr
  names(thr) <- sprintf('bin %d', 1:x$numbins)
  print(thr, print.gap=3)
  cat('\n')
  cat(sprintf('Best %s:', toupper(x$metric)), format(x$minse), '\n')
  cat('\n')

  invisible(x)
}


# summary.optbin:
# Print a table with the bin thresholds, averages, and (mean) squared
# error metric.  Set show.range to show the start/end points of the bin
# in the sorted data.

summary.optbin <- function(object, show.range=FALSE, ...) {
  if (!inherits(object, 'optbin')) {
    stop('Object must be from the optbin class')
  }

  if (show.range) {
    rng <- rep("", object$numbins)
    wpt <- floor(log10(length(object$x))) + 1
    for (b in 1:object$numbins) {
      if (1 == b) {
        stpt <- 1
      } else {
        stpt <- object$breaks[b-1] + 1
      }
      rng[b] <- sprintf('%*d:%*d', wpt, stpt, wpt, object$breaks[b])
    }
    df <- data.frame(bin=1:object$numbins, threshold=object$thr,
                     mean=object$binavg, metric=object$binse, range=rng)
  } else {
    df <- data.frame(bin=1:object$numbins, threshold=object$thr,
                     mean=object$binavg, metric=object$binse)
  }
  colnames(df)[4] <- object$metric

  cat('\n')
  cat(sprintf('Optimal %s Binning of Data with %d Elements\n',
              toupper(object$metric), length(object$x)))
  print(df, print.gap=3, row.names=FALSE)
  cat('\n')
}


# plot.optbin:
# Draw a plot with the points color-coded according to bin (col here is a
# vector of color names to use, or an internal list if left NULL), as well
# as the bin averages.  The x axis labels will be the bin thresholds.  Do
# not pass ann or xaxt as arguments.

plot.optbin <- function(x, col=NULL, main="Binned Observations", ...) {
  if (!inherits(x, 'optbin')) {
    stop('Object must be from the optbin class')
  }

  binID <- rep(1, length(x$x))
  for (b in 2:x$numbins) {
    binID[(x$breaks[b-1]+1):x$breaks[b]] <- b
  }

  if (is.null(col)) {
    # Should be color-blind friendly.  Same as list in hist.optbin.
    clrs <- c('#DF6BBD', '#AA0A3C', '#FF825F', '#0A9B4B',
              '#87D0CB', '#00A0FA', '#005AC8', '#8214A0')
  } else {
    clrs <- col
  }
  while (length(clrs) < x$numbins) {
    clrs <- c(clrs, clrs)
  }

  plot(x$x, col=clrs[binID], xaxt='n', ann=FALSE, ...)
  title(main, xlab='thresholds', ylab='observations')
  axis(1, at=x$breaks, labels=format(signif(x$thr, digits=3)))

  stpt <- 1
  for (b in 1:x$numbins) {
    lines(c(stpt,x$breaks[b]),c(x$binavg[b],x$binavg[b]), col=clrs[b])
    stpt <- x$breaks[b] + 1
  }

  if (x$numbins <= 10) {
    for (b in 1:x$numbins) {
      abline(v=x$breaks[b], col='grey75', lty=3)
    }
  }
}


# hist.optbin:
# Draw a histogram with bars under the x axis marking the bins, and also
# show the bin averages.  bincols # is a vector of color names to use to
# distinguish bins, or we use an internal list if left NULL.  Other
# arguments are passed through normally to the histogram function.

hist.optbin <- function(x, bincol=NULL, main=NULL, xlab=NULL, ...) {
  if (!inherits(x, 'optbin')) {
    stop('Object must be from the optbin class')
  }

  if (is.null(bincol)) {
    # Same as list in plot.optbin.
    clrs <- c('#DF6BBD', '#AA0A3C', '#FF825F', '#0A9B4B',
              '#87D0CB', '#00A0FA', '#005AC8', '#8214A0')
  } else {
    clrs <- bincol
  }
  while (length(clrs) < x$numbins) {
    clrs <- c(clrs, clrs)
  }

  srcvar <- as.character(x$call['x'])
  if (is.null(main)) {
    main <- paste("Histogram and Binning of", srcvar)
  }
  if (is.null(xlab)) {
    xlab <- srcvar
  }

  h <- hist(x$x, main=main, xlab=xlab, ...)
  stpt <- min(x$x)
  # Try to put the bars a little under the axis.  This might not always
  # give a clean separation.
  yoff <- -strheight("Test") / 4
  for (i in 1:x$numbins) {
    lines(c(stpt, x$thr[i]), c(yoff,yoff), col=clrs[i], lwd=3)
    abline(v=x$binavg[i], col=clrs[i], lwd=2)
    if (i < x$numbins) {
      abline(v=x$thr[i], col='grey75', lwd=1, lty=3)
    }
    stpt <- x$thr[i]
  }

  invisible(h)
}
