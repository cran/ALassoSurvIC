// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fun_hcross
NumericMatrix fun_hcross(NumericMatrix x);
RcppExport SEXP _ALassoSurvIC_fun_hcross(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fun_hcross(x));
    return rcpp_result_gen;
END_RCPP
}
// fun_subless
LogicalMatrix fun_subless(NumericVector u, NumericVector lessthan);
RcppExport SEXP _ALassoSurvIC_fun_subless(SEXP uSEXP, SEXP lessthanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lessthan(lessthanSEXP);
    rcpp_result_gen = Rcpp::wrap(fun_subless(u, lessthan));
    return rcpp_result_gen;
END_RCPP
}
// fun_sublr
LogicalMatrix fun_sublr(NumericVector u, NumericVector l, NumericVector r);
RcppExport SEXP _ALassoSurvIC_fun_sublr(SEXP uSEXP, SEXP lSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(fun_sublr(u, l, r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ALassoSurvIC_fun_hcross", (DL_FUNC) &_ALassoSurvIC_fun_hcross, 1},
    {"_ALassoSurvIC_fun_subless", (DL_FUNC) &_ALassoSurvIC_fun_subless, 2},
    {"_ALassoSurvIC_fun_sublr", (DL_FUNC) &_ALassoSurvIC_fun_sublr, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ALassoSurvIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
