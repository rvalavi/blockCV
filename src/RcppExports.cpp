// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// nndm_cpp
Rcpp::NumericMatrix nndm_cpp(Rcpp::NumericMatrix X, Rcpp::NumericVector Gjstar, Rcpp::NumericVector Gij, double rmin, double phi, int imin, int jmin, double min_train);
RcppExport SEXP _blockCV_nndm_cpp(SEXP XSEXP, SEXP GjstarSEXP, SEXP GijSEXP, SEXP rminSEXP, SEXP phiSEXP, SEXP iminSEXP, SEXP jminSEXP, SEXP min_trainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Gjstar(GjstarSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Gij(GijSEXP);
    Rcpp::traits::input_parameter< double >::type rmin(rminSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type imin(iminSEXP);
    Rcpp::traits::input_parameter< int >::type jmin(jminSEXP);
    Rcpp::traits::input_parameter< double >::type min_train(min_trainSEXP);
    rcpp_result_gen = Rcpp::wrap(nndm_cpp(X, Gjstar, Gij, rmin, phi, imin, jmin, min_train));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_blockCV_nndm_cpp", (DL_FUNC) &_blockCV_nndm_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_blockCV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
