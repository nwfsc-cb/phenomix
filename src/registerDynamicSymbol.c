// RegisteringDynamic Symbols
// from Rcpp https://github.com/RcppCore/Rcpp/issues/636
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_phenomix(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
