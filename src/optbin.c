
/*****
      optbin.c -
      Optimized binning with an efficient (nbin * nval^2) search to minimize
      the (mean) squared error of the data against the bin average.

      Shared library for R, not an executable.

      c 2019-2020 Primordial Machine Vision Systems
*****/


#include "optbin.h"


/**** R Registration ****/

static const R_CallMethodDef callMethods[] = {
	{"C_optbin", (DL_FUNC) &C_optbin, 4},
	{NULL, NULL, 0}
};

void R_init_optbin(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}


/**** Macros ****/

/* triangular indexing - see calc_se for a description */

/* matrix index for bin between stpt and endpt (incl., > stpt) of npt total */
#define BIN_TOIND(stpt, endpt, npt)                                \
	((((npt) - 1 - (stpt)) * ((npt) - (stpt)) / 2) - ((npt) - 1 - (endpt)) - 1)

/* number of points in (M)SE cache for npt data points
   cast to double to handle squaring although this may lose precision */
#define SE_SIZE(npt)                                               \
	((double) (npt) * ((double) (npt) - 1.0) / 2.0)



/**** Prototypes ****/

static void optbin_cached(int, double *, size_t, int, size_t **, double **, double *);
static void optbin_raw(int, double *, size_t, int, size_t **, double **, double *);

static void cache_se(double *, size_t, int, double **);
static void calc_sums(double *, size_t, double **, double **);
static inline double calc_se(double *, double *, size_t, size_t, int);
static void allocate_optbin(int, size_t, size_t **, double **,
                            size_t ***, double **, double **);
static SEXP set_results(SEXP, SEXP, int, size_t *, double *, double, SEXP);



/**** Public Interface ****/

/***
    C_optbin:  Break a vector of values data, already sorted, into a number
               of bins numbins, minimizing either the squared difference 
               between each point in the bin and the average over the bin, 
               or that difference divided by the number of points in the bin
               (ie. the squared error or mean squared error, which is also 
               the variance).  The routine will use up to cachemem bytes if
               possible to pre-calculate bins and save substantial run time.
    args:      data - sorted vector of double values to bin
               numbins - number of bins to create
               usemse - true to minimize mean squared error, false just SE
               cachemem - use cache if would be less than this many bytes
    returns:   R list with elements $breaks, $se
***/
SEXP C_optbin(SEXP data, SEXP numbins, SEXP usemse, SEXP cachemem) {
	SEXP res;                             /* our calculation */
	double minse;                         /* best (M)SE metric at chosen bins */
	double *binse;                        /* (M)SE per bin */
	size_t *endpts;                       /* last element in each bin */
	size_t nval;                          /* number of values in data */
	int nbin;                             /* number of bins to create */
	int mse;                              /* true to use MSE, 0 SE */

	endpts = NULL;
	binse = NULL;

	if (length(data) < (2 * asReal(numbins))) {
		error("too few data points for bins (need at least 2 per)");
	}
	mse = Rf_asLogical(usemse);

	nval = (size_t) length(data);
	nbin = (int) asReal(numbins);

	if (asReal(cachemem) <= (SE_SIZE(nval) * sizeof(double))) {
		optbin_raw(nbin, REAL(data), nval, mse, &endpts, &binse, &minse);
	} else {
		optbin_cached(nbin, REAL(data), nval, mse, &endpts, &binse, &minse);
	}

	res = PROTECT(allocVector(VECSXP, 8));
	res = set_results(numbins, data, mse, endpts, binse, minse, res);
	/* Bugfix 13 Sep 22: Originally unprotected in set_results, but automated
	   checks complain about mis-match in this function.  Decrease UNPROTECT
	   in set_results. */
	UNPROTECT(1);
	return res;
}


/**** Implementation ****/

/***
    optbin_cached:  Pick the bin partitions that minimizes the (M)SE within
                    each bin.  Allocates a large cache holding the cost of
                    all potential bins to reduce run time.  endpt is an
                    array that stores the last point in each bin (incl.); 
                    the last will always be nval-1.  Will call R's error()
                    to abort if there's a problem.
    args:           nbin - number of bins to split data into
                    data - sorted values
                    nval - number of values
                    mse - true to optimize mean squared error, 0 SE
                    binends - bin endpoints (incl.)
                    binse - (M)SE for each bin
                    minse - best (M)SE for chosen partition
    modifies:  binends, minse, binse
***/
static void optbin_cached(int nbin, double *data, size_t nval, int mse,
                          size_t **binends, double **binse, double *minse) {
	size_t **endpt;                       /* endpoints of bins */
	double *currse;                       /* new (M)SE with extra bin, per s */
	double *basese;                       /* (M)SE previous bin, per st pt */
	double *se;                           /* cached (M)SE */
	size_t s;                             /* start of left bin */
	size_t e;                             /* end of left bin */
	size_t slo;                           /* first start point to check */
	size_t shi;                           /* last start point (excl.) to check */
	size_t seind;                         /* index into se array */
	int b;                                /* bin counter */

	basese = NULL;
	currse = NULL;
	se = NULL;
	endpt = NULL;

	/* Simplify an exhaustive check of all possible splits with the recursion
	   for bin b from start/end points s and e
	     binning(b,s,e) = min se(s,e) + binning(b-1,e+1,nval-1)
	   where we check all possible endpoints e.  The seed is a single bin
	   from each starting point to the last
	     binning(1,s,nval-1) = se(s,nval-1)
	   The best 3-way split is a new bin on the left and the best 2-way split
	   on the right, where we have to check all possible middle points for
	   the best.  Final complexity is O(nbin * nval * nval).  Note we've
	   inverted b below - we count down from nbin-1.
	*/

	cache_se(data, nval, mse, &se);

	allocate_optbin(nbin, nval, binends, binse, &endpt, &basese, &currse);

	/* Seed with a single bin from each starting point to the end. */
	for (s=0; s<(nval-1); s++) {
		seind = BIN_TOIND(s, nval-1, nval);
		basese[s] = se[seind];
		endpt[s][nbin-1] = nval - 1;
	}

	/* Add bins to the left by scanning over all intermediate endpoints
	   and picking the best.  We write into endpt with a different bin
	   index so there's no collision, but would overwrite the best (M)SE
	   so need a current version.  You could use a scalar since basese/currse
	   only depend on s, but that blows out the runtime - there's some
	   parallelization in the CPU that the arrays allow? */
	for (b=nbin-2; 0<=b; b--) {
		R_CheckUserInterrupt();

		slo = 2 * b;
		/* Lose 2 per previous bin (nbin - 1 - b) and - 1 for the second point. */
		shi = nval - 2 * (nbin - 1 - b) - 1;
		for (s=slo; s<shi; s++) {
			currse[s] = DBL_MAX;
			for (e=s+1; e<(shi+1); e++) {
				seind = BIN_TOIND(s, e, nval);
				if ((se[seind] + basese[e+1]) < currse[s]) {
					endpt[s][b] = e;
					currse[s] = se[seind] + basese[e+1];
				}
			}

			/* We only need the first point for the final bin.  Break by changing
			   shi, which also cuts short the following copy loop. */
			if ((0 == b) && (slo == s)) {
				shi = 1;
			}

			/* A PELT test would go here, but on test data it screened out
			   very few points and the extra checks increased the run time. */
		}

		/* Copy current pass to base.  Could also set up a ping-pong array
		   to avoid the copy. */
		for (s=slo; s<shi; s++) {
			basese[s] = currse[s];
		}
	}

	*minse = basese[0];
	/* Untangle the endpoints and copy them into the result.  We simply
	   take the best endpoint for a given starting point, and set the new
	   start to the next point. */
	s = 0;
	for (b=0; b<nbin; b++) {
		(*binends)[b] = endpt[s][b];
		seind = BIN_TOIND(s, endpt[s][b], nval);
		(*binse)[b] = se[seind];
		s = endpt[s][b] + 1;
	}
}

/***
    optbin_raw:  Variant of optbin_cached that calculates the (M)SE each
                 time rather than using a cache.  This would be needed for
                 long vectors whose cache won't fit in memory.  endpt is
                 an array that stores the last point in each bin (incl.);
                 the last value will always be nval-1.  Will call R's error()
                 to abort if there's a problem.
    args:        nbin - number of bins to split data into
                 data - sorted values
                 nval - number of values
                 mse - true to optimize mean squared error, 0 SE
                 binends - bin endpoints (incl.)
                 binse - (M)SE for each bin
                 minse - best (M)SE for chosen partition
    modifies:  binends, minse
***/
static void optbin_raw(int nbin, double *data, size_t nval, int mse,
                       size_t **binends, double **binse, double *minse) {
	size_t **endpt;                       /* endpoints of bins */
	double *currse;                       /* new (M)SE with extra bin, per s */
	double *basese;                       /* (M)SE previous bin, per st pt */
	double *sum;                          /* cumulative sum data (1 based) */
	double *sumsq;                        /* cum sum data squared (1 based) */
	double se;                            /* bin metric */
	size_t s;                             /* start of left bin */
	size_t e;                             /* end of left bin */
	size_t slo;                           /* first start point to check */
	size_t shi;                           /* last start point (excl.) to check */
	int b;                                /* bin counter */

	basese = NULL;
	currse = NULL;
	sum = NULL;
	sumsq = NULL;
	endpt = NULL;

	/* This is a copy of optbin_cached except we make sum and sumsq and replace
	   se cache indices with calls to binse. */

	calc_sums(data, nval, &sum, &sumsq);

	allocate_optbin(nbin, nval, binends, binse, &endpt, &basese, &currse);

	for (s=0; s<nval; s++) {
		basese[s] = calc_se(sum, sumsq, s, nval-1, mse);
		endpt[s][nbin-1] = nval - 1;
	}

	for (b=nbin-2; 0<=b; b--) {
		R_CheckUserInterrupt();

		slo = 2 * b;
		/* Lose 2 per previous bin (nbin - 1 - b) and - 1 for the second point. */
		shi = nval - 2 * (nbin - 1 - b) - 1;
		for (s=slo; s<shi; s++) {
			currse[s] = DBL_MAX;
			for (e=s+1; e<(shi+1); e++) {
				se = calc_se(sum, sumsq, s, e, mse);
				if ((se + basese[e+1]) < currse[s]) {
					endpt[s][b] = e;
					currse[s] = se + basese[e+1];
				}
			}

			if ((0 == b) && (slo == s)) {
				shi = 1;
			}
		}

		/* Copy current pass to base.  Could also set up a ping-pong array
		   to avoid the copy. */
		for (s=slo; s<shi; s++) {
			basese[s] = currse[s];
		}
	}

	*minse = basese[0];

	s = 0;
	for (b=0; b<nbin; b++) {
		(*binends)[b] = endpt[s][b];
		(*binse)[b] = calc_se(sum, sumsq, s, endpt[s][b], mse);
		s = endpt[s][b] + 1;
	}
}


/**** Internal Functions ****/

/***
    cache_se:  Calculate the (mean) squared error for each possible bin (s,e)
               where the start and endpoints are inclusive.  se uses 
               triangular indexing to reduce memory; use the BIN_TOIND macro
               to convert the endpoints to an index into the array.  Will
               call R's error() to abort if there's a problem.
    args:      data - sorted values to bin
               nval - number of data values
               mse - true to calculate mean squared error, 0 SE
               se - cost metric for bin (re-alloc if non-NULL)
    modifies:  se
***/
static void cache_se(double *data, size_t nval, int mse, double **se) {
	double *sum;                          /* cumulative sum data (1 based) */
	double *sumsq;                        /* cum sum data squared (1 based) */
	size_t s, e;                          /* start, endpoint (incl.) of bin */

	sum = NULL;
	sumsq = NULL;

	Free(*se);

	/* Allocating a square array for se is fastest - simple indexing - but
	   we only fill the upper triangle so half the memory is wasted.  The 
	   diagonal is skipped because we require two points in a bin; or, the
	   (M)SE for a single point bin would always be 0.

	   Upper triangle indexing uses Gauss' formula to sum a series 1..N,
	   or N * (N+1) / 2.  These are also called triangular numbers.  Counting
	   starts in the lower right corner:
               c 0  1  2  3  4  5  6  7  8  9
            r 0  . 37 38 39 40 41 42 43 44 45     9 r' = (nval - 1) - r
              1     . 29 30 31 32 33 34 35 36     8
              2        . 22 23 24 25 26 27 28     7
              3           . 16 17 18 19 20 21     6
              4              . 11 12 13 14 15     5
              5                 .  7  8  9 10     4
              6                    .  4  5  6     3
              7                       .  2  3     2
              8                          .  1     1
              9                             .     0
	   The rows count backwards from the number of values (less 1 for 0 index)
	   and r * (r' + 1) / 2 is the rightmost index on the row.  The column
	   adjustment subtracts (nval - 1 -c ), and we lose another 1 for 0 offset,
	   giving the equation in the BIN_TOIND macro:
	     index = ((nval - 1 - r) * (nval - r) / 2) - (nval - 1 - c) - 1
	   You could reverse the order per row (ex. r=5: 10 9 8 7) by using a
	   column adjustment of (c - r - 1), but this goes against the column loop
	   and costs about 10% run time.
	*/

	*se = (double *) R_alloc(SE_SIZE(nval), sizeof(**se));
	if (NULL == *se) {
		error("cache allocation failed (out of memory)");
	}

	calc_sums(data, nval, &sum, &sumsq);

	for (s=0; s<nval-1; s++) {
		for (e=s+1; e<nval; e++) {
			(*se)[BIN_TOIND(s, e, nval)] = calc_se(sum, sumsq, s, e, mse);
		}
	}
}

/***
    calc_sums:  Generate the cumulative sum and sum-of-squares of the data
                array, ie. sum[j] = data[0] + data[1] + ... + data[j] and
                sumsq[j] = data[0]*data[0] + ... + data[j]*data[j].  Note
                arrays are 1 based and sum[0] == sumsq[0] == 0, so arrays
                are allocated to hold nval+1 elements.  Will call R's error()
                to abort if there's a problem (allocation).
    args:       data - analysis data
                nval - length of data
                sum - cumulative sum of data (re-alloc if non-NULL)
                sumsq - cumulative sum of data squared (re-alloc if non-NULL)
    modifies:  sum, sumsq
***/
static void calc_sums(double *data, size_t nval, double **sum, double **sumsq) {
	size_t i;

	Free(*sum);
	Free(*sumsq);

	*sum = (double *) R_alloc(nval+1, sizeof(**sum));
	if (NULL == *sum) {
		error("allocation of cumulative sum failed");
	}

	*sumsq = (double *) R_alloc(nval+1, sizeof(**sumsq));
	if (NULL == *sumsq) {
		error("allocation of cumulative sum of squares failed");
	}

	(*sum)[0] = 0.0;
	(*sumsq)[0] = 0.0;
	for (i=0; i<nval; i++) {
		(*sum)[i+1] = (*sum)[i] + data[i];
		if ((0 < i) && (data[i] < data[i-1])) {
			error("data has not been sorted");
		}
		(*sumsq)[i+1] = (*sumsq)[i] + (data[i] * data[i]);
	}
}

/***
    calc_se:  Calculate the (mean) squared error over the bin.  Assumes valid
              points (ie. stpt < endpt).
    args:     sum - cumulative sum of data values (index 1 based)
              sumsq - cumulative sum of data values squared (index 1 based)
              stpt - first point in bin
              endpt - last point in bin (incl.)
              mse - true to calculate MSE, 0 for SE
    returns:   (mean) squared error
***/
static inline double calc_se(double *sum, double *sumsq,
                             size_t stpt, size_t endpt, int mse) {
	double se;                            /* bin metric */
	double mu;                            /* mean over bin */
	size_t bsize;                         /* number of points in bin */

	/* Precision might cause errors with this one-pass approach, but a
	   two pass calculating first the average and then summing the differences
	   would be too expensive.

	   se = sum((x - mu) * (x - mu)) = sum(x*x - 2*mu*x + mu*mu)
	   se = sum(x*x) - 2 * mu * sum(x) + N * mu * mu
	   se = sum(x*x) - 2 * mu * N * sum(x) / N + N * mu * mu
	   se = sum(x*x) - 2 * mu * N * mu + N * mu * mu
	   se = sum(x*x) - N * mu * mu
	*/

	bsize = endpt - stpt + 1;
	mu = (sum[endpt+1] - sum[stpt]) / bsize;
	se = (sumsq[endpt+1] - sumsq[stpt]) - bsize * mu * mu;
	if (mse) {
		se /= bsize;
	}

	return se;
}

/***
    allocate_optbin:  Common memory operations for the two optbin variants.
                      Will call R's error() to abort if there's a problem.
                      Arrays are initialized if needed (binends, binse, 
                      currse need not be).
    args:             nbin - number of bins we're creating
                      nval - number of data values
                      binends - last bin point (incl.)
                      binse - bin (M)SE
                      endpt - bin endpoints for given starting point
                      basese - previous pass (M)SE per start point
                      currse - this pass (M)SE per start point
    returns:   0 if successful
               < 0 on failure (value depends on error)
    modifies:  binends, binse, endpt, basese, currse
***/
static void allocate_optbin(int nbin, size_t nval, size_t **binends,
                            double **binse, size_t ***endpt,
                            double **basese, double **currse) {
	int b;                                /* bin counter */
	size_t i;

	*binends = (size_t *) R_alloc(nbin, sizeof(**binends));
	if (NULL == *binends) {
		error("allocation of bin endpoints failed");
	}

	*binse = (double *) R_alloc(nbin, sizeof(**binse));
	if (NULL == *binse) {
		error("allocation of bin metrics failed");
	}

	*endpt = (size_t **) R_alloc(nval, sizeof(**endpt));
	if (NULL == *endpt) {
		error("allocation of bin endpoints failed");
	}
	for (i=0; i<nval; i++) {
		(*endpt)[i] = (size_t *) R_alloc(nbin, sizeof(***endpt));
		if (NULL == (*endpt)[i]) {
			error("allocation of endpoints for bin %ld failed", (long) i);
		}
	}

	*basese = (double *) R_alloc(nval, sizeof(**basese));
	if (NULL == *basese) {
		error("allocation of base-pass SE array failed");
	}
	*currse = (double *) R_alloc(nval, sizeof(**currse));
	if (NULL == *currse) {
		error("allocation of current-pass SE array failed");
	}

	for (i=0; i<nval; i++) {
		(*basese)[i] = DBL_MAX;
		for (b=0; b<nbin; b++) {
			(*endpt)[i][b] = 0;
		}
	}
}

/***
    set_results:  Store our results plus some derived values in the results.
    args:         numbins - number of bins we've created
                  data - sorted data we've used
                  mse - true if using MSE to optimize bins, 0 for SE
                  endpts - upper index in data of bins (incl.)
                  binse - (M)SE per bin
                  minse - (M)SE summed over all bins
                  res - results of operation
    returns:   0 if successful
               < 0 on failure (value depends on error)
    returns:   modified res
    modifies:  res
***/
static SEXP set_results(SEXP numbins, SEXP data, int mse, size_t *endpts,
                        double *binse, double minse, SEXP res) {
	SEXP names;                           /* identifier strings res entries */
	SEXP class;                           /* class name assigned to res */
	SEXP vends;                           /* bin endpoints */
	SEXP vse;                             /* (M)SE per bin */
	SEXP vthr;                            /* bin thresholds */
	SEXP vavg;                            /* average value over bin */
	double *data_ends;                    /* contents of vends */
	double *data_se;                      /* contents of vse */
	double *data_thr;                     /* contents of vthr */
	double *data_avg;                     /* contents of vavg */
	int nbin;                             /* number of bins to create */
	size_t s;                             /* first point of bin */
	int b;                                /* bin counter */
	size_t i;

	names = PROTECT(allocVector(STRSXP, 8));
	
	nbin = (int) asReal(numbins);

	vends = PROTECT(allocVector(REALSXP, nbin));
	data_ends = REAL(vends);
	vse = PROTECT(allocVector(REALSXP, nbin));
	data_se = REAL(vse);
	vthr = PROTECT(allocVector(REALSXP, nbin));
	data_thr = REAL(vthr);
	vavg = PROTECT(allocVector(REALSXP, nbin));
	data_avg = REAL(vavg);

	for (b=0; b<nbin; b++) {
		/* + 1 because R is 1 based. */
		data_ends[b] = endpts[b] + 1;
		data_se[b] = binse[b];
		data_thr[b] = REAL(data)[endpts[b]];
		data_avg[b] = 0.0;
		if (0 == b) {
			s = 0;
		} else {
			s = endpts[b-1] + 1;
		}
		for (i=s; i<=endpts[b]; i++) {
			data_avg[b] += REAL(data)[i];
		}
		data_avg[b] /= (endpts[b] - s + 1);
	}

	i = 0;
	SET_VECTOR_ELT(res, i, data);
	SET_STRING_ELT(names, i, mkChar("x"));

	i += 1;
	SET_VECTOR_ELT(res, i, numbins);
	SET_STRING_ELT(names, i, mkChar("numbins"));

	i += 1;
	if (mse) {
		SET_VECTOR_ELT(res, i, mkString("mse"));
	} else {
		SET_VECTOR_ELT(res, i, mkString("se"));
	}
	SET_STRING_ELT(names, i, mkChar("metric"));

	i += 1;
	SET_VECTOR_ELT(res, i, ScalarReal(minse));
	SET_STRING_ELT(names, i, mkChar("minse"));

	i += 1;
	SET_VECTOR_ELT(res, i, vthr);
	SET_STRING_ELT(names, i, mkChar("thr"));

	i += 1;
	SET_VECTOR_ELT(res, i, vavg);
	SET_STRING_ELT(names, i, mkChar("binavg"));

	i += 1;
	SET_VECTOR_ELT(res, i, vse);
	SET_STRING_ELT(names, i, mkChar("binse"));

	i += 1;
	SET_VECTOR_ELT(res, i, vends);
	SET_STRING_ELT(names, i, mkChar("breaks"));

	setAttrib(res, R_NamesSymbol, names);

	class = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("optbin"));
	classgets(res, class);

	UNPROTECT(6);
	return res;
}


