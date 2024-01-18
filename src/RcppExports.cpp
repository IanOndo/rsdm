// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_mahalanobis
arma::vec arma_mahalanobis(arma::mat const& x, arma::vec const& center, arma::mat const& cov);
RcppExport SEXP _rsdm_arma_mahalanobis(SEXP xSEXP, SEXP centerSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_mahalanobis(x, center, cov));
    return rcpp_result_gen;
END_RCPP
}
// arma_calc_nt1
arma::vec arma_calc_nt1(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_calc_nt1(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_calc_nt1(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_calc_mic_nt1
Rcpp::IntegerVector arma_calc_mic_nt1(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_calc_mic_nt1(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_calc_mic_nt1(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_calc_nt2
arma::vec arma_calc_nt2(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_calc_nt2(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_calc_nt2(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_calc_mic_nt2
arma::uvec arma_calc_mic_nt2(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_calc_mic_nt2(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_calc_mic_nt2(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_micdet
SEXP arma_micdet(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_micdet(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_micdet(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_micdet_raster
Rcpp::S4 arma_micdet_raster(Rcpp::S4& trg, Rcpp::S4& ref);
RcppExport SEXP _rsdm_arma_micdet_raster(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_micdet_raster(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_exdet
arma::vec arma_exdet(arma::mat& trg, arma::mat& ref);
RcppExport SEXP _rsdm_arma_exdet(SEXP trgSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_exdet(trg, ref));
    return rcpp_result_gen;
END_RCPP
}
// arma_exdet_raster
Rcpp::S4 arma_exdet_raster(Rcpp::S4& trg, Rcpp::S4& ref, bool compute_mic);
RcppExport SEXP _rsdm_arma_exdet_raster(SEXP trgSEXP, SEXP refSEXP, SEXP compute_micSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type trg(trgSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type ref(refSEXP);
    Rcpp::traits::input_parameter< bool >::type compute_mic(compute_micSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_exdet_raster(trg, ref, compute_mic));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rsdm_arma_mahalanobis", (DL_FUNC) &_rsdm_arma_mahalanobis, 3},
    {"_rsdm_arma_calc_nt1", (DL_FUNC) &_rsdm_arma_calc_nt1, 2},
    {"_rsdm_arma_calc_mic_nt1", (DL_FUNC) &_rsdm_arma_calc_mic_nt1, 2},
    {"_rsdm_arma_calc_nt2", (DL_FUNC) &_rsdm_arma_calc_nt2, 2},
    {"_rsdm_arma_calc_mic_nt2", (DL_FUNC) &_rsdm_arma_calc_mic_nt2, 2},
    {"_rsdm_arma_micdet", (DL_FUNC) &_rsdm_arma_micdet, 2},
    {"_rsdm_arma_micdet_raster", (DL_FUNC) &_rsdm_arma_micdet_raster, 2},
    {"_rsdm_arma_exdet", (DL_FUNC) &_rsdm_arma_exdet, 2},
    {"_rsdm_arma_exdet_raster", (DL_FUNC) &_rsdm_arma_exdet_raster, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_rsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}