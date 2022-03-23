// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// permutedFC_RCPP
List permutedFC_RCPP(arma::mat logCPM, int NB, int sEachp, arma::vec weight);
RcppExport SEXP _SSPT_permutedFC_RCPP(SEXP logCPMSEXP, SEXP NBSEXP, SEXP sEachpSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type logCPM(logCPMSEXP);
    Rcpp::traits::input_parameter< int >::type NB(NBSEXP);
    Rcpp::traits::input_parameter< int >::type sEachp(sEachpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(permutedFC_RCPP(logCPM, NB, sEachp, weight));
    return rcpp_result_gen;
END_RCPP
}
// permutedPertScore_RCPP
List permutedPertScore_RCPP(const List& BminsI, const CharacterVector& expressedG, arma::mat LogCPM, int NB, int sEachp, arma::vec weight);
RcppExport SEXP _SSPT_permutedPertScore_RCPP(SEXP BminsISEXP, SEXP expressedGSEXP, SEXP LogCPMSEXP, SEXP NBSEXP, SEXP sEachpSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type BminsI(BminsISEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type expressedG(expressedGSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LogCPM(LogCPMSEXP);
    Rcpp::traits::input_parameter< int >::type NB(NBSEXP);
    Rcpp::traits::input_parameter< int >::type sEachp(sEachpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(permutedPertScore_RCPP(BminsI, expressedG, LogCPM, NB, sEachp, weight));
    return rcpp_result_gen;
END_RCPP
}
// permutedPertScore_RCPP_indiPathway_Alt
List permutedPertScore_RCPP_indiPathway_Alt(const arma::mat& X, const CharacterVector& pathwayG, const CharacterVector& expressedG, const List& permutedFC, int newS);
RcppExport SEXP _SSPT_permutedPertScore_RCPP_indiPathway_Alt(SEXP XSEXP, SEXP pathwayGSEXP, SEXP expressedGSEXP, SEXP permutedFCSEXP, SEXP newSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type pathwayG(pathwayGSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type expressedG(expressedGSEXP);
    Rcpp::traits::input_parameter< const List& >::type permutedFC(permutedFCSEXP);
    Rcpp::traits::input_parameter< int >::type newS(newSSEXP);
    rcpp_result_gen = Rcpp::wrap(permutedPertScore_RCPP_indiPathway_Alt(X, pathwayG, expressedG, permutedFC, newS));
    return rcpp_result_gen;
END_RCPP
}
// permutedPertScore_RCPP_indiPathway
List permutedPertScore_RCPP_indiPathway(arma::mat X, const CharacterVector& pathwayG, const CharacterVector& expressedG, const List& permutedFC, int newS);
RcppExport SEXP _SSPT_permutedPertScore_RCPP_indiPathway(SEXP XSEXP, SEXP pathwayGSEXP, SEXP expressedGSEXP, SEXP permutedFCSEXP, SEXP newSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type pathwayG(pathwayGSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type expressedG(expressedGSEXP);
    Rcpp::traits::input_parameter< const List& >::type permutedFC(permutedFCSEXP);
    Rcpp::traits::input_parameter< int >::type newS(newSSEXP);
    rcpp_result_gen = Rcpp::wrap(permutedPertScore_RCPP_indiPathway(X, pathwayG, expressedG, permutedFC, newS));
    return rcpp_result_gen;
END_RCPP
}
// ssPertScore_RCPP
List ssPertScore_RCPP(const List& BminsI, arma::mat weightedFC, const CharacterVector& expressedG, const CharacterVector& sample);
RcppExport SEXP _SSPT_ssPertScore_RCPP(SEXP BminsISEXP, SEXP weightedFCSEXP, SEXP expressedGSEXP, SEXP sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type BminsI(BminsISEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weightedFC(weightedFCSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type expressedG(expressedGSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type sample(sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(ssPertScore_RCPP(BminsI, weightedFC, expressedG, sample));
    return rcpp_result_gen;
END_RCPP
}
// ssPertScore_RCPP_oneP
NumericVector ssPertScore_RCPP_oneP(arma::mat adjMatrix, const CharacterVector& pathwayG, arma::mat weightedFC, const CharacterVector& expressedG, const CharacterVector& sample);
RcppExport SEXP _SSPT_ssPertScore_RCPP_oneP(SEXP adjMatrixSEXP, SEXP pathwayGSEXP, SEXP weightedFCSEXP, SEXP expressedGSEXP, SEXP sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type adjMatrix(adjMatrixSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type pathwayG(pathwayGSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weightedFC(weightedFCSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type expressedG(expressedGSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type sample(sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(ssPertScore_RCPP_oneP(adjMatrix, pathwayG, weightedFC, expressedG, sample));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSPT_permutedFC_RCPP", (DL_FUNC) &_SSPT_permutedFC_RCPP, 4},
    {"_SSPT_permutedPertScore_RCPP", (DL_FUNC) &_SSPT_permutedPertScore_RCPP, 6},
    {"_SSPT_permutedPertScore_RCPP_indiPathway_Alt", (DL_FUNC) &_SSPT_permutedPertScore_RCPP_indiPathway_Alt, 5},
    {"_SSPT_permutedPertScore_RCPP_indiPathway", (DL_FUNC) &_SSPT_permutedPertScore_RCPP_indiPathway, 5},
    {"_SSPT_ssPertScore_RCPP", (DL_FUNC) &_SSPT_ssPertScore_RCPP, 4},
    {"_SSPT_ssPertScore_RCPP_oneP", (DL_FUNC) &_SSPT_ssPertScore_RCPP_oneP, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSPT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
