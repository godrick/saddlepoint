// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/saddlepoint_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// makeAdaptorUsingRfunctions
Rcpp::XPtr<Adaptor> makeAdaptorUsingRfunctions(Rcpp::Function r_function);
RcppExport SEXP _saddlepoint_makeAdaptorUsingRfunctions(SEXP r_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type r_function(r_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(makeAdaptorUsingRfunctions(r_function));
    return rcpp_result_gen;
END_RCPP
}
// makeVectorSubsetByIndicesAdaptor
Rcpp::XPtr<Adaptor> makeVectorSubsetByIndicesAdaptor(Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices);
RcppExport SEXP _saddlepoint_makeVectorSubsetByIndicesAdaptor(SEXP indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> >::type indices(indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(makeVectorSubsetByIndicesAdaptor(indices));
    return rcpp_result_gen;
END_RCPP
}
// makeSavedVectorAdaptor
Rcpp::XPtr<Adaptor> makeSavedVectorAdaptor(vec fixed_parameter_values);
RcppExport SEXP _saddlepoint_makeSavedVectorAdaptor(SEXP fixed_parameter_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type fixed_parameter_values(fixed_parameter_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSavedVectorAdaptor(fixed_parameter_values));
    return rcpp_result_gen;
END_RCPP
}
// adapt_CGF
Rcpp::XPtr<CGF_with_AD> adapt_CGF(Rcpp::XPtr<CGF_with_AD> cgf, Rcpp::XPtr<Adaptor> adaptor);
RcppExport SEXP _saddlepoint_adapt_CGF(SEXP cgfSEXP, SEXP adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type adaptor(adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(adapt_CGF(cgf, adaptor));
    return rcpp_result_gen;
END_RCPP
}
// makeADFunNegll
Rcpp::XPtr<TMBad::ADFun<>> makeADFunNegll(const vec& tvec, const vec& theta, Rcpp::XPtr<CGF_with_AD> cgf, bool optimize);
RcppExport SEXP _saddlepoint_makeADFunNegll(SEXP tvecSEXP, SEXP thetaSEXP, SEXP cgfSEXP, SEXP optimizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< bool >::type optimize(optimizeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeADFunNegll(tvec, theta, cgf, optimize));
    return rcpp_result_gen;
END_RCPP
}
// makeADFunK1
Rcpp::XPtr<TMBad::ADFun<>> makeADFunK1(const vec& tvec, const vec& theta, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_makeADFunK1(SEXP tvecSEXP, SEXP thetaSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(makeADFunK1(tvec, theta, cgf));
    return rcpp_result_gen;
END_RCPP
}
// makeADFunIneqConstraint
Rcpp::XPtr<TMBad::ADFun<>> makeADFunIneqConstraint(const vec& tvec, const vec& theta, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_makeADFunIneqConstraint(SEXP tvecSEXP, SEXP thetaSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(makeADFunIneqConstraint(tvec, theta, cgf));
    return rcpp_result_gen;
END_RCPP
}
// makeADFunZerothNegll
Rcpp::XPtr<TMBad::ADFun<>> makeADFunZerothNegll(const vec& tvec, const vec& theta, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_makeADFunZerothNegll(SEXP tvecSEXP, SEXP thetaSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(makeADFunZerothNegll(tvec, theta, cgf));
    return rcpp_result_gen;
END_RCPP
}
// makeADFunCustom1Negll
Rcpp::XPtr<TMBad::ADFun<>> makeADFunCustom1Negll(const vec& tvec, const vec& theta, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_makeADFunCustom1Negll(SEXP tvecSEXP, SEXP thetaSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(makeADFunCustom1Negll(tvec, theta, cgf));
    return rcpp_result_gen;
END_RCPP
}
// computeCombinedGradient
Rcpp::List computeCombinedGradient(const vec& combined_vector, Rcpp::XPtr<TMBad::ADFun<>> adf);
RcppExport SEXP _saddlepoint_computeCombinedGradient(SEXP combined_vectorSEXP, SEXP adfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type combined_vector(combined_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<TMBad::ADFun<>> >::type adf(adfSEXP);
    rcpp_result_gen = Rcpp::wrap(computeCombinedGradient(combined_vector, adf));
    return rcpp_result_gen;
END_RCPP
}
// computeFuncT
Rcpp::List computeFuncT(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> cgf, bool optimize);
RcppExport SEXP _saddlepoint_computeFuncT(SEXP tvecSEXP, SEXP thetaSEXP, SEXP observationsSEXP, SEXP cgfSEXP, SEXP optimizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const vec& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< bool >::type optimize(optimizeSEXP);
    rcpp_result_gen = Rcpp::wrap(computeFuncT(tvec, theta, observations, cgf, optimize));
    return rcpp_result_gen;
END_RCPP
}
// computeZerothFuncT
Rcpp::List computeZerothFuncT(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_computeZerothFuncT(SEXP tvecSEXP, SEXP thetaSEXP, SEXP observationsSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const vec& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(computeZerothFuncT(tvec, theta, observations, cgf));
    return rcpp_result_gen;
END_RCPP
}
// computeNegll
Rcpp::List computeNegll(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_computeNegll(SEXP tvecSEXP, SEXP thetaSEXP, SEXP observationsSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const vec& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(computeNegll(tvec, theta, observations, cgf));
    return rcpp_result_gen;
END_RCPP
}
// computeZerothNegll
Rcpp::List computeZerothNegll(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> cgf);
RcppExport SEXP _saddlepoint_computeZerothNegll(SEXP tvecSEXP, SEXP thetaSEXP, SEXP observationsSEXP, SEXP cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const vec& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(computeZerothNegll(tvec, theta, observations, cgf));
    return rcpp_result_gen;
END_RCPP
}
// make_BinomialCGF
Rcpp::XPtr<CGF_with_AD> make_BinomialCGF();
RcppExport SEXP _saddlepoint_make_BinomialCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_BinomialCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_BinomialModelCGF
Rcpp::XPtr<CGF_with_AD> make_BinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_adaptor);
RcppExport SEXP _saddlepoint_make_BinomialModelCGF(SEXP n_adaptorSEXP, SEXP prob_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type n_adaptor(n_adaptorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type prob_adaptor(prob_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_BinomialModelCGF(n_adaptor, prob_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_PoissonCGF
Rcpp::XPtr<CGF_with_AD> make_PoissonCGF();
RcppExport SEXP _saddlepoint_make_PoissonCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_PoissonCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_PoissonModelCGF
Rcpp::XPtr<CGF_with_AD> make_PoissonModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor);
RcppExport SEXP _saddlepoint_make_PoissonModelCGF(SEXP lambda_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type lambda_adaptor(lambda_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_PoissonModelCGF(lambda_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_ExponentialCGF
Rcpp::XPtr<CGF_with_AD> make_ExponentialCGF();
RcppExport SEXP _saddlepoint_make_ExponentialCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_ExponentialCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_ExponentialModelCGF
Rcpp::XPtr<CGF_with_AD> make_ExponentialModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor);
RcppExport SEXP _saddlepoint_make_ExponentialModelCGF(SEXP lambda_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type lambda_adaptor(lambda_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_ExponentialModelCGF(lambda_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_GeometricCGF
Rcpp::XPtr<CGF_with_AD> make_GeometricCGF();
RcppExport SEXP _saddlepoint_make_GeometricCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_GeometricCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_GeometricModelCGF
Rcpp::XPtr<CGF_with_AD> make_GeometricModelCGF(Rcpp::XPtr<Adaptor> p_adaptor);
RcppExport SEXP _saddlepoint_make_GeometricModelCGF(SEXP p_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type p_adaptor(p_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_GeometricModelCGF(p_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_GammaCGF
Rcpp::XPtr<CGF_with_AD> make_GammaCGF();
RcppExport SEXP _saddlepoint_make_GammaCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_GammaCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_GammaModelCGF
Rcpp::XPtr<CGF_with_AD> make_GammaModelCGF(Rcpp::XPtr<Adaptor> shape_adaptor, Rcpp::XPtr<Adaptor> rate_adaptor);
RcppExport SEXP _saddlepoint_make_GammaModelCGF(SEXP shape_adaptorSEXP, SEXP rate_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type shape_adaptor(shape_adaptorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type rate_adaptor(rate_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_GammaModelCGF(shape_adaptor, rate_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_MultinomialCGF
Rcpp::XPtr<CGF_with_AD> make_MultinomialCGF();
RcppExport SEXP _saddlepoint_make_MultinomialCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_MultinomialCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_MultinomialModelCGF
Rcpp::XPtr<CGF_with_AD> make_MultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_vector_adaptor);
RcppExport SEXP _saddlepoint_make_MultinomialModelCGF(SEXP n_adaptorSEXP, SEXP prob_vector_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type n_adaptor(n_adaptorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type prob_vector_adaptor(prob_vector_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_MultinomialModelCGF(n_adaptor, prob_vector_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_SubunitaryMultinomialModelCGF
Rcpp::XPtr<CGF_with_AD> make_SubunitaryMultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_vector_adaptor);
RcppExport SEXP _saddlepoint_make_SubunitaryMultinomialModelCGF(SEXP n_adaptorSEXP, SEXP prob_vector_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type n_adaptor(n_adaptorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type prob_vector_adaptor(prob_vector_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_SubunitaryMultinomialModelCGF(n_adaptor, prob_vector_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_NegBinCGF
Rcpp::XPtr<CGF_with_AD> make_NegBinCGF();
RcppExport SEXP _saddlepoint_make_NegBinCGF() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_NegBinCGF());
    return rcpp_result_gen;
END_RCPP
}
// make_NegBinModelCGF
Rcpp::XPtr<CGF_with_AD> make_NegBinModelCGF(Rcpp::XPtr<Adaptor> r_adaptor, Rcpp::XPtr<Adaptor> p_adaptor);
RcppExport SEXP _saddlepoint_make_NegBinModelCGF(SEXP r_adaptorSEXP, SEXP p_adaptorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type r_adaptor(r_adaptorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Adaptor> >::type p_adaptor(p_adaptorSEXP);
    rcpp_result_gen = Rcpp::wrap(make_NegBinModelCGF(r_adaptor, p_adaptor));
    return rcpp_result_gen;
END_RCPP
}
// make_SumOfIIDCGF
Rcpp::XPtr<CGF_with_AD> make_SumOfIIDCGF(Rcpp::XPtr<CGF_with_AD> cgf, double n);
RcppExport SEXP _saddlepoint_make_SumOfIIDCGF(SEXP cgfSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(make_SumOfIIDCGF(cgf, n));
    return rcpp_result_gen;
END_RCPP
}
// make_SumOfIndependentCGF
Rcpp::XPtr<CGF_with_AD> make_SumOfIndependentCGF(Rcpp::List cgf_list);
RcppExport SEXP _saddlepoint_make_SumOfIndependentCGF(SEXP cgf_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type cgf_list(cgf_listSEXP);
    rcpp_result_gen = Rcpp::wrap(make_SumOfIndependentCGF(cgf_list));
    return rcpp_result_gen;
END_RCPP
}
// make_ConcatenationCGF
Rcpp::XPtr<CGF_with_AD> make_ConcatenationCGF(Rcpp::List cgf_length_list);
RcppExport SEXP _saddlepoint_make_ConcatenationCGF(SEXP cgf_length_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type cgf_length_list(cgf_length_listSEXP);
    rcpp_result_gen = Rcpp::wrap(make_ConcatenationCGF(cgf_length_list));
    return rcpp_result_gen;
END_RCPP
}
// make_LinearlyMappedCGF
Rcpp::XPtr<CGF_with_AD> make_LinearlyMappedCGF(Rcpp::XPtr<CGF_with_AD> cgf, mat Amat);
RcppExport SEXP _saddlepoint_make_LinearlyMappedCGF(SEXP cgfSEXP, SEXP AmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< mat >::type Amat(AmatSEXP);
    rcpp_result_gen = Rcpp::wrap(make_LinearlyMappedCGF(cgf, Amat));
    return rcpp_result_gen;
END_RCPP
}
// make_RandomlyStoppedSumCGF
Rcpp::XPtr<CGF_with_AD> make_RandomlyStoppedSumCGF(Rcpp::XPtr<CGF_with_AD> count_cgf, Rcpp::XPtr<CGF_with_AD> summand_cgf);
RcppExport SEXP _saddlepoint_make_RandomlyStoppedSumCGF(SEXP count_cgfSEXP, SEXP summand_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type count_cgf(count_cgfSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type summand_cgf(summand_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(make_RandomlyStoppedSumCGF(count_cgf, summand_cgf));
    return rcpp_result_gen;
END_RCPP
}
// make_IIDReplicatesCGF
Rcpp::XPtr<CGF_with_AD> make_IIDReplicatesCGF(Rcpp::XPtr<CGF_with_AD> cgf, double block_size);
RcppExport SEXP _saddlepoint_make_IIDReplicatesCGF(SEXP cgfSEXP, SEXP block_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type cgf(cgfSEXP);
    Rcpp::traits::input_parameter< double >::type block_size(block_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(make_IIDReplicatesCGF(cgf, block_size));
    return rcpp_result_gen;
END_RCPP
}
// K_impl
double K_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_K_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(K_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// K1_impl
vec K1_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_K1_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(K1_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// K2_impl
mat K2_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_K2_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(K2_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// ineq_constraint_impl
vec ineq_constraint_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_ineq_constraint_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(ineq_constraint_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// neg_ll_impl
double neg_ll_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_neg_ll_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_ll_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// func_T_impl
double func_T_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_func_T_impl(SEXP tvecSEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(func_T_impl(tvec, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// K4operatorAABB_impl
double K4operatorAABB_impl(vec tvec, mat a1, mat a2, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_K4operatorAABB_impl(SEXP tvecSEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< mat >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< mat >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(K4operatorAABB_impl(tvec, a1, a2, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// K3K3operatorAABBCC_impl
double K3K3operatorAABBCC_impl(vec tvec, mat a1, mat a2, mat a3, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf);
RcppExport SEXP _saddlepoint_K3K3operatorAABBCC_impl(SEXP tvecSEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP a3SEXP, SEXP parameter_vectorSEXP, SEXP base_cgfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< mat >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< mat >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< mat >::type a3(a3SEXP);
    Rcpp::traits::input_parameter< vec >::type parameter_vector(parameter_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<CGF_with_AD> >::type base_cgf(base_cgfSEXP);
    rcpp_result_gen = Rcpp::wrap(K3K3operatorAABBCC_impl(tvec, a1, a2, a3, parameter_vector, base_cgf));
    return rcpp_result_gen;
END_RCPP
}
// make_CustomVectorizedScalarCGF
Rcpp::XPtr<CGF_with_AD> make_CustomVectorizedScalarCGF(Rcpp::Function Kfunc, Rcpp::Function K1func, Rcpp::Function K2func, Rcpp::Function K3func, Rcpp::Function K4func, Rcpp::Function ineq_constraint_func);
RcppExport SEXP _saddlepoint_make_CustomVectorizedScalarCGF(SEXP KfuncSEXP, SEXP K1funcSEXP, SEXP K2funcSEXP, SEXP K3funcSEXP, SEXP K4funcSEXP, SEXP ineq_constraint_funcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type Kfunc(KfuncSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type K1func(K1funcSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type K2func(K2funcSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type K3func(K3funcSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type K4func(K4funcSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type ineq_constraint_func(ineq_constraint_funcSEXP);
    rcpp_result_gen = Rcpp::wrap(make_CustomVectorizedScalarCGF(Kfunc, K1func, K2func, K3func, K4func, ineq_constraint_func));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_saddlepoint_makeAdaptorUsingRfunctions", (DL_FUNC) &_saddlepoint_makeAdaptorUsingRfunctions, 1},
    {"_saddlepoint_makeVectorSubsetByIndicesAdaptor", (DL_FUNC) &_saddlepoint_makeVectorSubsetByIndicesAdaptor, 1},
    {"_saddlepoint_makeSavedVectorAdaptor", (DL_FUNC) &_saddlepoint_makeSavedVectorAdaptor, 1},
    {"_saddlepoint_adapt_CGF", (DL_FUNC) &_saddlepoint_adapt_CGF, 2},
    {"_saddlepoint_makeADFunNegll", (DL_FUNC) &_saddlepoint_makeADFunNegll, 4},
    {"_saddlepoint_makeADFunK1", (DL_FUNC) &_saddlepoint_makeADFunK1, 3},
    {"_saddlepoint_makeADFunIneqConstraint", (DL_FUNC) &_saddlepoint_makeADFunIneqConstraint, 3},
    {"_saddlepoint_makeADFunZerothNegll", (DL_FUNC) &_saddlepoint_makeADFunZerothNegll, 3},
    {"_saddlepoint_makeADFunCustom1Negll", (DL_FUNC) &_saddlepoint_makeADFunCustom1Negll, 3},
    {"_saddlepoint_computeCombinedGradient", (DL_FUNC) &_saddlepoint_computeCombinedGradient, 2},
    {"_saddlepoint_computeFuncT", (DL_FUNC) &_saddlepoint_computeFuncT, 5},
    {"_saddlepoint_computeZerothFuncT", (DL_FUNC) &_saddlepoint_computeZerothFuncT, 4},
    {"_saddlepoint_computeNegll", (DL_FUNC) &_saddlepoint_computeNegll, 4},
    {"_saddlepoint_computeZerothNegll", (DL_FUNC) &_saddlepoint_computeZerothNegll, 4},
    {"_saddlepoint_make_BinomialCGF", (DL_FUNC) &_saddlepoint_make_BinomialCGF, 0},
    {"_saddlepoint_make_BinomialModelCGF", (DL_FUNC) &_saddlepoint_make_BinomialModelCGF, 2},
    {"_saddlepoint_make_PoissonCGF", (DL_FUNC) &_saddlepoint_make_PoissonCGF, 0},
    {"_saddlepoint_make_PoissonModelCGF", (DL_FUNC) &_saddlepoint_make_PoissonModelCGF, 1},
    {"_saddlepoint_make_ExponentialCGF", (DL_FUNC) &_saddlepoint_make_ExponentialCGF, 0},
    {"_saddlepoint_make_ExponentialModelCGF", (DL_FUNC) &_saddlepoint_make_ExponentialModelCGF, 1},
    {"_saddlepoint_make_GeometricCGF", (DL_FUNC) &_saddlepoint_make_GeometricCGF, 0},
    {"_saddlepoint_make_GeometricModelCGF", (DL_FUNC) &_saddlepoint_make_GeometricModelCGF, 1},
    {"_saddlepoint_make_GammaCGF", (DL_FUNC) &_saddlepoint_make_GammaCGF, 0},
    {"_saddlepoint_make_GammaModelCGF", (DL_FUNC) &_saddlepoint_make_GammaModelCGF, 2},
    {"_saddlepoint_make_MultinomialCGF", (DL_FUNC) &_saddlepoint_make_MultinomialCGF, 0},
    {"_saddlepoint_make_MultinomialModelCGF", (DL_FUNC) &_saddlepoint_make_MultinomialModelCGF, 2},
    {"_saddlepoint_make_SubunitaryMultinomialModelCGF", (DL_FUNC) &_saddlepoint_make_SubunitaryMultinomialModelCGF, 2},
    {"_saddlepoint_make_NegBinCGF", (DL_FUNC) &_saddlepoint_make_NegBinCGF, 0},
    {"_saddlepoint_make_NegBinModelCGF", (DL_FUNC) &_saddlepoint_make_NegBinModelCGF, 2},
    {"_saddlepoint_make_SumOfIIDCGF", (DL_FUNC) &_saddlepoint_make_SumOfIIDCGF, 2},
    {"_saddlepoint_make_SumOfIndependentCGF", (DL_FUNC) &_saddlepoint_make_SumOfIndependentCGF, 1},
    {"_saddlepoint_make_ConcatenationCGF", (DL_FUNC) &_saddlepoint_make_ConcatenationCGF, 1},
    {"_saddlepoint_make_LinearlyMappedCGF", (DL_FUNC) &_saddlepoint_make_LinearlyMappedCGF, 2},
    {"_saddlepoint_make_RandomlyStoppedSumCGF", (DL_FUNC) &_saddlepoint_make_RandomlyStoppedSumCGF, 2},
    {"_saddlepoint_make_IIDReplicatesCGF", (DL_FUNC) &_saddlepoint_make_IIDReplicatesCGF, 2},
    {"_saddlepoint_K_impl", (DL_FUNC) &_saddlepoint_K_impl, 3},
    {"_saddlepoint_K1_impl", (DL_FUNC) &_saddlepoint_K1_impl, 3},
    {"_saddlepoint_K2_impl", (DL_FUNC) &_saddlepoint_K2_impl, 3},
    {"_saddlepoint_ineq_constraint_impl", (DL_FUNC) &_saddlepoint_ineq_constraint_impl, 3},
    {"_saddlepoint_neg_ll_impl", (DL_FUNC) &_saddlepoint_neg_ll_impl, 3},
    {"_saddlepoint_func_T_impl", (DL_FUNC) &_saddlepoint_func_T_impl, 3},
    {"_saddlepoint_K4operatorAABB_impl", (DL_FUNC) &_saddlepoint_K4operatorAABB_impl, 5},
    {"_saddlepoint_K3K3operatorAABBCC_impl", (DL_FUNC) &_saddlepoint_K3K3operatorAABBCC_impl, 6},
    {"_saddlepoint_make_CustomVectorizedScalarCGF", (DL_FUNC) &_saddlepoint_make_CustomVectorizedScalarCGF, 6},
    {NULL, NULL, 0}
};

void rtmb_set_shared_pointers();
RcppExport void R_init_saddlepoint(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  rtmb_set_shared_pointers();
}
