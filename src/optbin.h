
#include <R.h>
#include <Rinternals.h>

/* STRICT_R_HEADERS removes this in R_ext/Rconstants.h, but doesn't have a
   replacement for DBL_MAX. */
#include "float.h"

extern SEXP C_optbin(SEXP, SEXP, SEXP, SEXP);
