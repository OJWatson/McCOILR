// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// McCOIL_categorical_cpp
void McCOIL_categorical_cpp(Rcpp::List paramList);
RcppExport SEXP McCOILR_McCOIL_categorical_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    McCOIL_categorical_cpp(paramList);
    return R_NilValue;
END_RCPP
}
// McCOIL_proportional_cpp
void McCOIL_proportional_cpp(Rcpp::List paramList);
RcppExport SEXP McCOILR_McCOIL_proportional_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    McCOIL_proportional_cpp(paramList);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"McCOILR_McCOIL_categorical_cpp", (DL_FUNC) &McCOILR_McCOIL_categorical_cpp, 1},
    {"McCOILR_McCOIL_proportional_cpp", (DL_FUNC) &McCOILR_McCOIL_proportional_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_McCOILR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
