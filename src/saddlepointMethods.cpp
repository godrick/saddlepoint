# include "saddlepoint_types.h"
# include "parametric_submodelCGF.h"
# include "adaptor_using_r.h"
# include "tvec_hat_implicit_function.h"

 
# include "binomialCGF.h"
# include "poissonCGF.h"
# include "exponentialCGF.h"
# include "geometricCGF.h"
# include "gammaCGF.h"
# include "multinomialCGF.h"
# include "subunitaryMultinomialCGF.h"
 
 
# include "sumOfIID_CGF.h"
# include "sumOfIndependentCGF.h"
# include "concatenationCGF.h"
# include "linearly_mappedCGF.h"
# include "randomlyStoppedSumCGF.h"
# include "vectorToIIDCGF.h"  // replicates
 
 
# include "cgf_from_rfunctions.h"



using namespace saddlepoint::CGFs_with_AD;
using namespace saddlepoint::CGFs_with_Rcpp;


// [[Rcpp::export]]
Rcpp::XPtr<Adaptor> makeAdaptorUsingRfunctions(Rcpp::Function r_function) {
  Adaptor* adaptor_ptr = new AdaptorUsingRFunctions(r_function);
  Rcpp::XPtr<Adaptor> ptr(adaptor_ptr);
  attach_attributes(ptr, r_function);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<Adaptor> makeVectorSubsetByIndicesAdaptor(Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices) {
  return Rcpp::XPtr<Adaptor>(new SubsetVectorByIndicesAdaptor(indices.array() - 1));
}
// [[Rcpp::export]]
Rcpp::XPtr<Adaptor> makeSavedVectorAdaptor(vec fixed_parameter_values) {
  return Rcpp::XPtr<Adaptor>(new SavedVectorAdaptor(fixed_parameter_values));
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> adapt_CGF(Rcpp::XPtr<CGF_with_AD> cgf, Rcpp::XPtr<Adaptor> adaptor)
{
  auto modelCGF_ptr = new parametric_submodelCGF(cgf, adaptor);
  Rcpp::XPtr<CGF_with_AD> ptr(modelCGF_ptr);
  attach_attributes(ptr, cgf, adaptor);
  return ptr;
}









// [[Rcpp::export]]
Rcpp::XPtr< TMBad::ADFun<> > makeADFunK1(const vec& tvec,
                                         const vec& theta,
                                         Rcpp::XPtr<CGF_with_AD> cgf){
  const CGF_with_AD* cgf_ptr = cgf.get();
  auto func = [cgf_ptr, tvec_size=tvec.size(), theta_size=theta.size()](const std::vector<a_scalar>& x) {
    a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), tvec_size);
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + tvec_size, theta_size);
    a_vector result = cgf_ptr->K1(tvec_ad, theta_ad);
    std::vector<a_scalar> result_std(result.data(), result.data() + result.size());
    return result_std;
  };
  
  std::vector<double> combined_input(tvec.size() + theta.size());
  std::copy(tvec.begin(), tvec.end(), combined_input.begin());
  std::copy(theta.begin(), theta.end(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr< TMBad::ADFun<> > makeADFunNegll(const vec& tvec,
                                            const vec& theta,
                                            Rcpp::XPtr<CGF_with_AD> cgf){
  const CGF_with_AD* cgf_ptr = cgf.get();
  auto func = [cgf_ptr, tvec_size=tvec.size(), theta_size=theta.size()](const std::vector<a_scalar>& x) {
    a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), tvec_size);
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + tvec_size, theta_size);
    a_scalar result = cgf_ptr->neg_ll(tvec_ad, theta_ad);
    return std::vector<a_scalar>{result};
  };
  
  std::vector<double> combined_input(tvec.size() + theta.size());
  std::copy(tvec.begin(), tvec.end(), combined_input.begin());
  std::copy(theta.begin(), theta.end(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr< TMBad::ADFun<> > makeADFunIneqConstraint(const vec& tvec,
                                                     const vec& theta,
                                                     Rcpp::XPtr<CGF_with_AD> cgf){
  const CGF_with_AD* cgf_ptr = cgf.get();
  auto func = [cgf_ptr, tvec_size=tvec.size(), theta_size=theta.size()](const std::vector<a_scalar>& x) {
    a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), tvec_size);
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + tvec_size, theta_size);
    a_vector result = cgf_ptr->ineq_constraint(tvec_ad, theta_ad);
    std::vector<a_scalar> result_std(result.data(), result.data() + result.size());
    return result_std;
  };
  
  std::vector<double> combined_input(tvec.size() + theta.size());
  std::copy(tvec.begin(), tvec.end(), combined_input.begin());
  std::copy(theta.begin(), theta.end(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List computeCombinedGradient(const vec& combined_vector, Rcpp::XPtr<TMBad::ADFun<>> adf){
  std::vector<double> combined_input(combined_vector.size());
  std::copy(combined_vector.begin(), combined_vector.end(), combined_input.begin());
  
  return Rcpp::List::create(Rcpp::Named("objective") = (*adf)(combined_input),
                            Rcpp::Named("gradient") = adf->Jacobian(combined_input));
}










// [[Rcpp::export]]
Rcpp::List computeFuncT(const vec& tvec,
                        const vec& theta,
                        const vec& observations,
                        Rcpp::XPtr<CGF_with_AD> cgf){
  a_vector tvec_ad = tvec.cast<a_scalar>();
  a_vector observations_ad = observations.cast<a_scalar>();
  const CGF_with_AD* cgf_ptr = cgf.get();

  auto func = [cgf_ptr, tvec_ad, observations_ad](const std::vector<a_scalar>& x) {
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data(), x.size());
    a_vector tvec_ad_hat = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_ptr);
    a_scalar result = cgf_ptr->func_T(tvec_ad_hat, theta_ad);
    return std::vector<a_scalar>{result};
  };

  std::vector<double> theta_std(theta.data(), theta.data() + theta.size());
  TMBad::ADFun<> adf(func, theta_std);
  TMBad::ADFun<> jac_adf = adf.JacFun();

  return Rcpp::List::create(Rcpp::Named("value") = adf(theta_std),
                            Rcpp::Named("gradient") = adf.Jacobian(theta_std),
                            Rcpp::Named("hessian") = jac_adf.Jacobian(theta_std)
                            );

}


// [[Rcpp::export]]
Rcpp::List computeNegll(const vec& tvec,
                        const vec& theta,
                        const vec& observations,
                        Rcpp::XPtr<CGF_with_AD> cgf){
  a_vector tvec_ad = tvec.cast<a_scalar>();
  a_vector observations_ad = observations.cast<a_scalar>();
  const CGF_with_AD* cgf_ptr = cgf.get();

  auto func = [cgf_ptr, tvec_ad, observations_ad](const std::vector<a_scalar>& x) {
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data(), x.size());
    a_vector tvec_ad_hat = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_ptr);
    a_scalar result = cgf_ptr->neg_ll(tvec_ad_hat, theta_ad);
    return std::vector<a_scalar>{result};
  };

  std::vector<double> theta_std(theta.data(), theta.data() + theta.size());
  TMBad::ADFun<> adf(func, theta_std);
  TMBad::ADFun<> jac_adf = adf.JacFun();

  return Rcpp::List::create(Rcpp::Named("value") = adf(theta_std),
                            Rcpp::Named("gradient") = adf.Jacobian(theta_std),
                            Rcpp::Named("hessian") = jac_adf.Jacobian(theta_std)
  );

}

// RCPP_MODULE(thetaGradientFunctions) {
//   using namespace Rcpp;
// 
//   function("computeFuncT", &computeFuncT,
//            "Computes function T and its gradient with respect to theta.");
//   function("computeNegll", &computeNegll,
//            "Computes saddlepoint negative log-likelihood and its gradient with respect theta.");
// }





// // Template to handle common logic
// template<typename Func>
// std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
// computeGeneral(const vec& tvec, const vec& theta, const vec& observations, CGF_with_AD* cgf_ptr, Func funcMember) {
//   a_vector tvec_ad = tvec.cast<a_scalar>();
//   a_vector observations_ad = observations.cast<a_scalar>();
//   
//   auto func = [cgf_ptr, tvec_ad, observations_ad, funcMember](const std::vector<a_scalar>& x) {
//     a_vector theta_ad = Eigen::Map<const a_vector>(x.data(), x.size());
//     a_vector tvec_ad_hat = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_ptr);
//     a_scalar result = funcMember(tvec_ad_hat, theta_ad);
//     return std::vector<a_scalar>{result};
//   };
//   
//   std::vector<double> theta_std(theta.data(), theta.data() + theta.size());
//   TMBad::ADFun<> adf(func, theta_std);
//   TMBad::ADFun<> jac_adf = adf.JacFun();
//   
//   return std::make_tuple(
//     adf(theta_std),
//     adf.Jacobian(theta_std),
//     jac_adf.Jacobian(theta_std)
//   );
// }
// 
// // [[Rcpp::export]]
// Rcpp::List computeFuncT(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> modelCGF) {
//   auto tempfuncT = [cgf = modelCGF.get()](const a_vector& arg1, const a_vector& arg2) {
//     return cgf->func_T(arg1, arg2);
//   };
//   auto result = computeGeneral(tvec, theta, observations, modelCGF.get(), tempfuncT);
//   return Rcpp::List::create(
//     Rcpp::Named("value") = std::get<0>(result),
//     Rcpp::Named("gradient") = std::get<1>(result),
//     Rcpp::Named("hessian") = std::get<2>(result)
//   );
// }
// 
// // [[Rcpp::export]]
// Rcpp::List computeNegll(const vec& tvec, const vec& theta, const vec& observations, Rcpp::XPtr<CGF_with_AD> modelCGF) {
//   auto tempneg_ll = [cgf = modelCGF.get()](const a_vector& arg1, const a_vector& arg2) {
//     return cgf->neg_ll(arg1, arg2);
//   };
//   auto result = computeGeneral(tvec, theta, observations, modelCGF.get(), tempneg_ll);
//   return Rcpp::List::create(
//     Rcpp::Named("value") = std::get<0>(result),
//     Rcpp::Named("gradient") = std::get<1>(result),
//     Rcpp::Named("hessian") = std::get<2>(result)
//   );
// }










// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_BinomialCGF(){
  CGF_with_AD* CGF_base_ptr = new BinomialCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_BinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_adaptor){
  CGF_with_AD* binomial_model_cgf = new BinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), new ScalarAdaptorFromVectorAdaptor(prob_adaptor));
  Rcpp::XPtr<CGF_with_AD> ptr(binomial_model_cgf);
  attach_attributes(ptr, n_adaptor, prob_adaptor);
  return ptr;
}


// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_PoissonCGF(){
  CGF_with_AD* CGF_base_ptr = new PoissonCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_PoissonModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor){
  CGF_with_AD* poisson_model_cgf = new PoissonModelCGF( new ScalarAdaptorFromVectorAdaptor(lambda_adaptor) );
  Rcpp::XPtr<CGF_with_AD> ptr(poisson_model_cgf);
  attach_attributes(ptr, lambda_adaptor);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_ExponentialCGF(){
  CGF_with_AD* CGF_base_ptr = new ExponentialCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_ExponentialModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor){
  CGF_with_AD* exponential_model_cgf = new ExponentialModelCGF( new ScalarAdaptorFromVectorAdaptor(lambda_adaptor));
  Rcpp::XPtr<CGF_with_AD> ptr(exponential_model_cgf);
  attach_attributes(ptr, lambda_adaptor);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_GeometricCGF(){
  CGF_with_AD* CGF_base_ptr = new GeometricCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_GeometricModelCGF(Rcpp::XPtr<Adaptor> p_adaptor){
  CGF_with_AD* geometric_model_cgf = new GeometricModelCGF( new ScalarAdaptorFromVectorAdaptor(p_adaptor) );
  Rcpp::XPtr<CGF_with_AD> ptr(geometric_model_cgf);
  attach_attributes(ptr, p_adaptor);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_GammaCGF(){
  CGF_with_AD* CGF_base_ptr = new GammaCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_GammaModelCGF(Rcpp::XPtr<Adaptor> shape_adaptor, Rcpp::XPtr<Adaptor> rate_adaptor){
  CGF_with_AD* gamma_model_cgf = new GammaModelCGF( new ScalarAdaptorFromVectorAdaptor(shape_adaptor), new ScalarAdaptorFromVectorAdaptor(rate_adaptor) );
  Rcpp::XPtr<CGF_with_AD> ptr(gamma_model_cgf);
  attach_attributes(ptr, shape_adaptor, rate_adaptor);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_MultinomialCGF(){
  CGF_with_AD* CGF_base_ptr = new MultinomialCGF();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_MultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_vector_adaptor){
  CGF_with_AD* multinomial_model_cgf = new MultinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), prob_vector_adaptor);
  Rcpp::XPtr<CGF_with_AD> ptr(multinomial_model_cgf);
  attach_attributes(ptr, n_adaptor, prob_vector_adaptor);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_SubunitaryMultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_vector_adaptor){
  CGF_with_AD* cgf = new SubunitaryMultinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), prob_vector_adaptor);
  Rcpp::XPtr<CGF_with_AD> ptr(cgf);
  attach_attributes(ptr, n_adaptor, prob_vector_adaptor);
  return ptr;
}



























// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_SumOfIIDCGF(Rcpp::XPtr<CGF_with_AD> cgf, double n){
  CGF_with_AD* cgf_ptr = new SumOfIID_CGF(cgf, n);
  Rcpp::XPtr<CGF_with_AD> ptr(cgf_ptr);
  attach_attributes(ptr, cgf, n);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_SumOfIndependentCGF(Rcpp::List cgf_list) {
  // Populate a container of cgfs from cgf_list
  std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
  for (auto i = cgf_list.begin(); i != cgf_list.end(); ++i) {
    Rcpp::XPtr<CGF_with_AD> xp_cgf = *i;
    cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(xp_cgf.get()));
  }
  // Pass this container into SumOfIndependentCGF
  CGF_with_AD* CGF_base_ptr = new SumOfIndependentCGF(cgfs);
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  attach_attributes(ptr, cgf_list);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_ConcatenationCGF(Rcpp::List cgf_length_list) {
  std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
  Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> lengths(cgf_length_list.size());
  //// Populate the containers from cgf_length_list

  for (int i = 0; i < cgf_length_list.size(); ++i) {
    Rcpp::List currentList = Rcpp::as<Rcpp::List>(cgf_length_list[i]);
    Rcpp::XPtr<CGF_with_AD> xp_cgf = Rcpp::as<Rcpp::XPtr<CGF_with_AD>>(currentList[0]);
    Eigen::Index length = Rcpp::as<Eigen::Index>(currentList[1]);
    cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(xp_cgf.get()));
    lengths[i] = length;
  }

  //// Pass these vectors into ConcatenationCGF
  CGF_with_AD* CGF_base_ptr = new ConcatenationCGF(cgfs, lengths);
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);

  attach_attributes(ptr, cgf_length_list);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_LinearlyMappedCGF(Rcpp::XPtr<CGF_with_AD> cgf, mat Amat){
  CGF_with_AD* CGF_ptr = new LinearlyMappedCGF(cgf, Amat);
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_ptr);
  attach_attributes(ptr, cgf, Amat);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_RandomlyStoppedSumCGF(Rcpp::XPtr<CGF_with_AD> count_cgf, Rcpp::XPtr<CGF_with_AD> summand_cgf){
  CGF_with_AD* rss_cgf = new RandomlyStoppedSumCGF(count_cgf, summand_cgf);
  Rcpp::XPtr<CGF_with_AD> ptr(rss_cgf);
  attach_attributes(ptr, count_cgf, summand_cgf);
  return ptr;
}
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_IIDReplicatesCGF(Rcpp::XPtr<CGF_with_AD> cgf, double block_size){
  CGF_with_AD* vector_to_iid_cgf = new IIDReplicates_CGF(cgf, block_size);
  Rcpp::XPtr<CGF_with_AD> ptr(vector_to_iid_cgf);
  attach_attributes(ptr, cgf);
  return ptr;
}






















// TO DO: The following functions of class CGF_with_AD should probably 
// be exposed as a collection by RCPP_MODULE. Here, we define 
// individual wrapper functions and export them separately. 
// The issue with RCPP_MODULE, we will need to find a way that 
// the module definition correctly handles the overloads and expose only the double versions of these functions.
// [[Rcpp::export]]
double K_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
  try {
    return base_cgf->K(tvec, parameter_vector);
  } catch (...) {
    Rcpp::stop("An error occurred in the function 'K'. Check the arguments.");
  }
}
// [[Rcpp::export]]
vec K1_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){ return base_cgf->K1(tvec, parameter_vector);}
// [[Rcpp::export]]
mat K2_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){return base_cgf->K2(tvec, parameter_vector);}
// [[Rcpp::export]]
vec ineq_constraint_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){return base_cgf->ineq_constraint(tvec, parameter_vector);}
// [[Rcpp::export]]
double neg_ll_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){ return base_cgf->neg_ll(tvec, parameter_vector);}
// [[Rcpp::export]]
double func_T_impl(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){ return base_cgf->func_T(tvec, parameter_vector);}
// [[Rcpp::export]]
double K4operatorAABB_impl(vec tvec, mat a1, mat a2, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){return base_cgf->K4operatorAABB(tvec, a1, a2, parameter_vector);}
// [[Rcpp::export]]
double K3K3operatorAABBCC_impl(vec tvec, mat a1, mat a2, mat a3, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){ return base_cgf->K3K3operatorAABBCC(tvec, a1, a2, a3, parameter_vector);}






// Create a CGF object from R functions
// The function takes in 6 R functions for K, K1, K2, K3, K4, and ineq_constraint
// and returns a pointer to a CGF object
// [[Rcpp::export]]
Rcpp::XPtr<CGF_with_AD> make_CustomVectorizedScalarCGF(Rcpp::Function Kfunc, Rcpp::Function K1func, Rcpp::Function K2func,
                                                       Rcpp::Function K3func, Rcpp::Function K4func, Rcpp::Function ineq_constraint_func) {
  CGF_with_AD_from_r_functions* cgf = new CGF_with_AD_from_r_functions(Kfunc, K1func, K2func, K3func, K4func, ineq_constraint_func);
  Rcpp::XPtr<CGF_with_AD> ptr(cgf);
  attach_attributes(ptr, Kfunc, K1func, K2func, K3func, K4func, ineq_constraint_func);
  return ptr;
}






// # include "saddlepoint_types.h"
// # include "atomic_funcs.h"
// 
// # include "atomic_implicit_function.h"
// 
// # include "suppliedFunctions.h"
// # include "parametric_submodelCGF.h"
// # include "AdaptorsWithR.h"
// 
// # include "binomialCGF.h"
// # include "poissonCGF.h"
// # include "exponentialCGF.h"
// # include "geometricCGF.h"
// # include "gammaCGF.h"
// # include "multinomialCGF.h"
// # include "subunitaryMultinomialCGF.h"
// 
// # include "sumOfIID_CGF.h"
// # include "sumOfIndependentCGF.h"
// # include "linearly_mappedCGF.h"
// # include "randomlyStoppedSumCGF.h"
// # include "scalarToIIDCGF.h"
// # include "vectorToIIDCGF.h"
// # include "concatenationCGF.h"
// 
// 
// 
// using namespace saddlepoint::CGFs_with_AD;
// using namespace saddlepoint::CGFs_with_Rcpp;
// using namespace saddlepoint::atomic_funcs;
// 
// /*
//  // ' Create an Adaptor Using R Functions
//  // '
//  // ' This function creates an `adaptor` object using two R functions.
//  // ' The `adaptor` object is returned as an external pointer (`XPtr`).
//  // '
//  // ' This function is primarily intended for internal use.
//  // '
//  // ' @param fn An R function for computing `h(theta)`. This function should return a vector.
//  // ' @param grad_fn An R function for computing the gradient of `h(theta)` with respect to `theta`. This function should return a matrix with `length(h(theta))` rows and `length(theta)` columns.
//  // '
//  // ' @return An `XPtr` to an `Adaptor` object.
//  // ' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<Adaptor> makeAdaptorUsingRfunctions(Rcpp::Function fn, Rcpp::Function grad_fn)
// {
//   SuppliedFunctions supplied_functions(fn, grad_fn);
//   auto atomic_ptr = new atomic_class_for_supplied_functions<SuppliedFunctions> ("supplied_functions_atomic", supplied_functions);
//   Rcpp::XPtr<atomic_class_for_supplied_functions<SuppliedFunctions>> atomic_Xptr(atomic_ptr);
//   Adaptor* adaptor_ptr = new AdaptorFromAtomicFunctionWithRFunctions(fn, atomic_ptr);
//   Rcpp::XPtr<Adaptor> ptr(adaptor_ptr);
//   attach_attributes(ptr, atomic_Xptr, fn, grad_fn);
//   return ptr;
// }
// /*
//  //' Create a Subvector Adaptor
//  //'
//  //' This function creates an `adaptor` object that extracts a subvector of a given length, starting at a specified position.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param pos The starting position of the subvector (1-indexed).
//  //' @param len The length of the subvector.
//  //'
//  //' @return An `XPtr` to an `adaptor` object.
//  //' @keywords internal
//  //'
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<Adaptor> makeSubvectorAdaptor(Eigen::Index pos, Eigen::Index len) {
//   return Rcpp::XPtr<Adaptor>(new SubvectorAdaptor(pos - 1, len));
// }
// /*
//  //' Create a Vector Subset By Indices Adaptor
//  //'
//  //' This function creates an `adaptor` object that extracts a subset of a vector based on given indices.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param indices A vector of indices for the subset.
//  //'
//  //' @return An `XPtr` to an `adaptor` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<Adaptor> makeVectorSubsetByIndicesAdaptor(Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices) {
//   return Rcpp::XPtr<Adaptor>(new SubsetVectorByIndicesAdaptor(indices.array() - 1));
// }
// /*
//  //' Create a Saved Vector Adaptor
//  //'
//  //' This function creates an `adaptor` object that holds fixed parameter values.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param fixed_parameter_values A vector of fixed parameter values.
//  //'
//  //' @return An `XPtr` to an `adaptor` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<Adaptor> makeSavedVectorAdaptor(vec fixed_parameter_values) {
//   return Rcpp::XPtr<Adaptor>(new SavedVectorAdaptor(fixed_parameter_values));
// }
// /*
//  //' Adapt a CGF Object
//  //'
//  //' This function attaches an `Adaptor` to a `CGF_with_AD` object.
//  //' The resulting object can be used as any other `CGF_with_AD` object (e.g., to create `CppAD::ADFun` objects).
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param base_cgf An `XPtr` to a `CGF_with_AD` object.
//  //' @param adaptor An `XPtr` to an `Adaptor` object.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> adapt_CGF(Rcpp::XPtr<CGF_with_AD> base_cgf, Rcpp::XPtr<Adaptor> adaptor)
// {
//   auto modelCGF_ptr = new parametric_submodelCGF(base_cgf, adaptor);
//   Rcpp::XPtr<CGF_with_AD> ptr(modelCGF_ptr);
//   attach_attributes(ptr, base_cgf, adaptor);
//   return ptr;
// }
// /*
//  //' Create an ADFun Object for K1
//  //'
//  //' This function creates a `CppAD::ADFun<double>` object for the K1 function of a `CGF` object.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param tvec A numeric vector.
//  //' @param theta A numeric vector.
//  //' @param modelCGF An `XPtr` to a `CGF_with_AD` object or of class `CGF`.
//  //'
//  //' @return An `XPtr` to a `CppAD::ADFun<double>` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr< CppAD::ADFun<double> > makeADFunK1(vec tvec,
//                                                vec theta,
//                                                Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   vec combined_vector(tvec.size() + theta.size());
//   combined_vector << tvec, theta;
//   //----------------------
//   a_vector a_combined_vector =  combined_vector.cast<CppAD::AD<double>>();
//   //----------------------
//   CppAD::Independent(a_combined_vector);
//   a_vector a_tvec = a_combined_vector.head(tvec.size());
//   a_vector a_distributional_pars = a_combined_vector.tail(theta.size());
//   //----------------------
//   a_vector a_K1 = modelCGF->K1(a_tvec, a_distributional_pars);
// 
//   CppAD::ADFun<double>* ADFun_ptr = new CppAD::ADFun<double>;
//   ADFun_ptr->Dependent(a_combined_vector, a_K1);
// 
//   Rcpp::XPtr< CppAD::ADFun<double> > ptr(ADFun_ptr);
//   attach_attributes(ptr, modelCGF);
//   return ptr;
// }
// /*
//  //' Compute K1 and its gradient
//  //'
//  //' This function computes the K1 function and its gradient for a given combined vector of `tvec` and `theta`.
//  //' The K1 function is computed using a `CppAD::ADFun<double>` object created by the `make.ADFunK1()`-R function.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param combined_vector A numeric vector of `tvec` and `theta` in that order.
//  //' @param ADfunK1 An `XPtr` to a `CppAD::ADFun<double>` object, typically the result of the `makeADFunK1()` function.
//  //'
//  //' @return A list with two elements: `fn`, the value of the K1 function, and `gr`, the gradient of the K1 function.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::List K1_with_gradient(vec combined_vector, Rcpp::XPtr< CppAD::ADFun<double> > ADfunK1){
//   return Rcpp::List::create(Rcpp::Named("fn") = ADfunK1->Forward(0, combined_vector),
//                             Rcpp::Named("gr") = ADfunK1->Jacobian(combined_vector));
// }
// /*
//  //' Create an ADFun Object for Saddlepoint negative log-likelihood
//  //'
//  //' This function creates a `CppAD::ADFun<double>` object for the negative log-likelihood function of a `CGF_with_AD` object.
//  //'
//  //' The function is intended for internal use.
//  //'
//  //' @param tvec A numeric vector.
//  //' @param theta A numeric vector.
//  //' @param modelCGF An `XPtr` to a `CGF_with_AD` object.
//  //'
//  //' @return An `XPtr` to a `CppAD::ADFun<double>` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr< CppAD::ADFun<double> > makeADFunNegll(vec tvec,
//                                                   vec theta,
//                                                   Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   //---------
//   vec combined_vector(tvec.size() + theta.size());
//   combined_vector << tvec, theta;
//   a_vector a_combined_vector =  combined_vector.cast<CppAD::AD<double>>();
//   //---------
//   CppAD::Independent(a_combined_vector);
//   a_vector a_tvec = a_combined_vector.head(tvec.size());
//   a_vector a_distributional_pars = a_combined_vector.tail(theta.size());
//   //---------
//   a_vector a_neg_ll(1);
//   a_neg_ll(0) = modelCGF->neg_ll(a_tvec, a_distributional_pars);
// 
//   CppAD::ADFun<double>* ADFun_ptr = new CppAD::ADFun<double>;
//   ADFun_ptr->Dependent(a_combined_vector, a_neg_ll);
// 
//   Rcpp::XPtr< CppAD::ADFun<double> > ptr(ADFun_ptr);
//   attach_attributes(ptr, modelCGF);
//   return ptr;
// }
// 
// // Zeroth order saddlepoint log-likelihood
// // [[Rcpp::export]]
// Rcpp::XPtr< CppAD::ADFun<double> > makeADFunZerothLL(vec tvec,
//                                                      vec theta,
//                                                      Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   //---------
//   vec combined_vector(tvec.size() + theta.size());
//   combined_vector << tvec, theta;
//   a_vector a_combined_vector =  combined_vector.cast<CppAD::AD<double>>();
//   //---------
//   CppAD::Independent(a_combined_vector);
//   a_vector a_tvec = a_combined_vector.head(tvec.size());
//   a_vector a_distributional_pars = a_combined_vector.tail(theta.size());
//   //---------
//   a_vector a_zeroth_order_ll(1);
//   a_zeroth_order_ll(0) = modelCGF->tilting_exponent(a_tvec, a_distributional_pars);
//   
//   CppAD::ADFun<double>* ADFun_ptr = new CppAD::ADFun<double>;
//   ADFun_ptr->Dependent(a_combined_vector, a_zeroth_order_ll);
//   
//   Rcpp::XPtr< CppAD::ADFun<double> > ptr(ADFun_ptr);
//   attach_attributes(ptr, modelCGF);
//   return ptr;
// }
// 
// /*
//  //' Compute log-likelihood and its gradient
//  //'
//  //' This function computes the log-likelihood function and its gradient for a given combined vector of `tvec` and `theta`.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param combined_vector A numeric vector of `tvec` and `theta` in that order.
//  //' @param ADfun_ll An `XPtr` to a `CppAD::ADFun<double>` object, typically the result of the `make.ADFunNegll()`-R function.
//  //'
//  //' @return A list with two elements: `objective`, the value of the log-likelihood function, and `gradient`, the gradient of the function.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::List ll_with_gradient(vec combined_vector, Rcpp::XPtr< CppAD::ADFun<double> > ADfun_ll)
// {
//   return Rcpp::List::create(Rcpp::Named("objective") = ADfun_ll->Forward(0, combined_vector),
//                             Rcpp::Named("gradient") = ADfun_ll->Jacobian(combined_vector));
// }
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// /*
//  //' Create an ADFun Object for inequality constraint
//  //'
//  //' This function creates a `CppAD::ADFun<double>` object for the inequality constraint function of a `CGF_with_AD` object.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param tvec A numeric vector.
//  //' @param theta A numeric vector.
//  //' @param modelCGF An `XPtr` to a `CGF_with_AD` object.
//  //'
//  //' @return An `XPtr` to a `CppAD::ADFun<double>` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr< CppAD::ADFun<double> > makeADFunIneqConstraint(vec tvec,
//                                                            vec theta,
//                                                            Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   vec combined_vector(tvec.size() + theta.size());
//   combined_vector << tvec, theta;
//   //----------------------
//   a_vector a_combined_vector =  combined_vector.cast<CppAD::AD<double>>();
//   //----------------------
//   CppAD::Independent(a_combined_vector);
//   a_vector a_tvec = a_combined_vector.head(tvec.size());
//   a_vector a_distributional_pars = a_combined_vector.tail(theta.size());
//   //----------------------
//   a_vector a_ineq_constraint = modelCGF->ineq_constraint(a_tvec, a_distributional_pars);
// 
//   CppAD::ADFun<double>* ADFun_ptr = new CppAD::ADFun<double>;
//   ADFun_ptr->Dependent(a_combined_vector, a_ineq_constraint);
// 
//   Rcpp::XPtr< CppAD::ADFun<double> > ptr(ADFun_ptr);
//   attach_attributes(ptr, modelCGF);
//   return ptr;
// }
// /*
//  //' Compute inequality constraint and its gradient
//  //'
//  //' This function computes the inequality constraint function and its gradient for a given combined vector of `tvec` and `theta`.
//  //'
//  //' This function is primarily intended for internal use.
//  //'
//  //' @param combined_vector A numeric vector of `tvec` and `theta` in that order.
//  //' @param ADfun_ineqConstraint An `XPtr` to a `CppAD::ADFun<double>` object, typically the result of the `make.ADFunIneqConstraint()`-R function.
//  //'
//  //' @return A list with two elements: `fn`, the value of the inequality constraint function, and `gr`, the gradient of the function.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::List ineqConstraint_with_gradient(vec combined_vector, Rcpp::XPtr< CppAD::ADFun<double> > ADfun_ineqConstraint)
// {
//   return Rcpp::List::create(Rcpp::Named("fn") = ADfun_ineqConstraint->Forward(0, combined_vector),
//                             Rcpp::Named("gr") = ADfun_ineqConstraint->Jacobian(combined_vector));
// }
// /*
//  //' Create a Binomial CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the binomial CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_BinomialCGF(){
//   CGF_with_AD* CGF_base_ptr = new BinomialCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Binomial Non-Identical CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the binomial non-identical CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_BinomialNonIdenticalCGF(){
//   CGF_with_AD* CGF_base_ptr = new BinomialNonIdenticalCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Binomial Model CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the binomial model CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @param n_adaptor An `XPtr` to an `Adaptor` object.
//  //' @param prob_adaptor An `XPtr` to an `Adaptor` object.
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_BinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> prob_adaptor){
//   CGF_with_AD* binomial_model_cgf = new BinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), new ScalarAdaptorFromVectorAdaptor(prob_adaptor));
//   Rcpp::XPtr<CGF_with_AD> ptr(binomial_model_cgf);
//   attach_attributes(ptr, n_adaptor, prob_adaptor);
//   return ptr;
// }
// /*
//  //' Create a Poisson CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Poisson CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_PoissonCGF(){
//   CGF_with_AD* CGF_base_ptr = new PoissonCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Poisson Non-Identical CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Poisson non-identical CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_PoissonNonIdenticalCGF(){
//   CGF_with_AD* CGF_base_ptr = new PoissonNonIdenticalCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Poisson Model CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Poisson model CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @param lambda_adaptor An `XPtr` to an `Adaptor` object.
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_PoissonModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor){
//   CGF_with_AD* poisson_model_cgf = new PoissonModelCGF( new ScalarAdaptorFromVectorAdaptor(lambda_adaptor) );
//   Rcpp::XPtr<CGF_with_AD> ptr(poisson_model_cgf);
//   attach_attributes(ptr, lambda_adaptor);
//   return ptr;
// }
// /*
//  //' Create a Exponential CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Exponential CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_ExponentialCGF(){
//   CGF_with_AD* CGF_base_ptr = new ExponentialCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Exponential Non-Identical CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Exponential non-identical CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_ExponentialNonIdenticalCGF(){
//   CGF_with_AD* CGF_base_ptr = new ExponentialNonIdenticalCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Exponential Model CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the Exponential model CGF.
//  //'
//  //' This function is intended for internal use.
//  //'
//  //' @param lambda_adaptor An `XPtr` to an `Adaptor` object.
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_ExponentialModelCGF(Rcpp::XPtr<Adaptor> lambda_adaptor){
//   CGF_with_AD* exponential_model_cgf = new ExponentialModelCGF( new ScalarAdaptorFromVectorAdaptor(lambda_adaptor));
//   Rcpp::XPtr<CGF_with_AD> ptr(exponential_model_cgf);
//   attach_attributes(ptr, lambda_adaptor);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GeometricCGF(){
//   CGF_with_AD* CGF_base_ptr = new GeometricCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GeometricNonIdenticalCGF(){
//   CGF_with_AD* CGF_base_ptr = new GeometricNonIdenticalCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GeometricModelCGF(Rcpp::XPtr<Adaptor> p_adaptor){
//   CGF_with_AD* geometric_model_cgf = new GeometricModelCGF( new ScalarAdaptorFromVectorAdaptor(p_adaptor) );
//   Rcpp::XPtr<CGF_with_AD> ptr(geometric_model_cgf);
//   attach_attributes(ptr, p_adaptor);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GeometricNonIdenticalModelCGF(Rcpp::XPtr<Adaptor> p_adaptor){
//   CGF_with_AD* CGF_base_ptr = new GeometricNonIdenticalModelCGF(p_adaptor);
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   attach_attributes(ptr, p_adaptor);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GammaCGF(){
//   CGF_with_AD* CGF_base_ptr = new GammaCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_GammaModelCGF(Rcpp::XPtr<Adaptor> shape_adaptor, Rcpp::XPtr<Adaptor> rate_adaptor){
//   CGF_with_AD* gamma_model_cgf = new GammaModelCGF( new ScalarAdaptorFromVectorAdaptor(shape_adaptor), new ScalarAdaptorFromVectorAdaptor(rate_adaptor) );
//   Rcpp::XPtr<CGF_with_AD> ptr(gamma_model_cgf);
//   attach_attributes(ptr, shape_adaptor, rate_adaptor);
//   return ptr;
// }
// /*
//  //' Create a Multinomial CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the multinomial (CGF).
//  //'
//  //' The function is intended for internal use.
//  //'
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_MultinomialCGF(){
//   CGF_with_AD* CGF_base_ptr = new MultinomialCGF();
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   return ptr;
// }
// /*
//  //' Create a Multinomial Model CGF Object
//  //'
//  //' This function creates a `CGF_with_AD` object for the multinomial model (CGF).
//  //'
//  //' The function is intended for internal use.
//  //'
//  //' @param n_adaptor An `XPtr` to an `Adaptor` object.
//  //' @param probVector_adaptor An `XPtr` to an `Adaptor` object.
//  //' @return An `XPtr` to a `CGF_with_AD` object.
//  //' @keywords internal
//  */
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_MultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> probVector_adaptor){
//   CGF_with_AD* multinomial_model_cgf = new MultinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), probVector_adaptor);
//   Rcpp::XPtr<CGF_with_AD> ptr(multinomial_model_cgf);
//   attach_attributes(ptr, n_adaptor, probVector_adaptor);
//   return ptr;
// }
// 
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_SubunitaryMultinomialModelCGF(Rcpp::XPtr<Adaptor> n_adaptor, Rcpp::XPtr<Adaptor> probVector_adaptor){
//   CGF_with_AD* model_cgf = new SubunitaryMultinomialModelCGF( new ScalarAdaptorFromVectorAdaptor(n_adaptor), probVector_adaptor);
//   Rcpp::XPtr<CGF_with_AD> ptr(model_cgf);
//   attach_attributes(ptr, n_adaptor, probVector_adaptor);
//   return ptr;
// }
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// double K(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   try {
//     return base_cgf->K(tvec, parameter_vector);
//   } catch (...) {
//     Rcpp::stop("An error occurred in the function 'K'. Check the arguments.");
//   }
// }
// 
// // [[Rcpp::export]]
// vec K1(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->K1(tvec, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// mat K2(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->K2(tvec, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// vec ineq_constraint(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->ineq_constraint(tvec, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// double neg_ll(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->neg_ll(tvec, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// double func_T(vec tvec, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->func_T(tvec, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// double K4operatorAABB(vec tvec, mat a1, mat a2, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->K4operatorAABB(tvec, a1, a2, parameter_vector);
// }
// 
// // [[Rcpp::export]]
// double K3K3operatorAABBCC(vec tvec, mat a1, mat a2, mat a3, vec parameter_vector, Rcpp::XPtr<CGF_with_AD> base_cgf){
//   return base_cgf->K3K3operatorAABBCC(tvec, a1, a2, a3, parameter_vector);
// }
// 
// 
// 
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_sum_of_iidCGF(Rcpp::XPtr<CGF_with_AD> base_cgf, double n){
//   CGF_with_AD* CGF_base_ptr = new SumOfIID_CGF(base_cgf, n);
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   attach_attributes(ptr, base_cgf, n);
//   return ptr;
// }
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_sum_of_independentCGF(Rcpp::List cgf_list) {
//   // Populate a container of cgfs from cgf_list
//   std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
//   for (auto i = cgf_list.begin(); i != cgf_list.end(); ++i) {
//     Rcpp::XPtr<CGF_with_AD> xp_cgf = *i;
//     cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(xp_cgf.get()));
//   }
//   // Pass this container into SumOfIndependentCGF
//   CGF_with_AD* CGF_base_ptr = new SumOfIndependentCGF(cgfs);
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
//   // return Rcpp::XPtr<CGF_with_AD>(new SumOfIndependentCGF(cgfs));
//   attach_attributes(ptr, cgf_list);
//   return ptr;
// }
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_ConcatenationCGF(Rcpp::List cgf_length_list) {
//   std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
//   Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> lengths(cgf_length_list.size());
//   //// Populate the containers from cgf_length_list
// 
//   for (int i = 0; i < cgf_length_list.size(); ++i) {
//     Rcpp::List currentList = Rcpp::as<Rcpp::List>(cgf_length_list[i]);
//     Rcpp::XPtr<CGF_with_AD> xp_cgf = Rcpp::as<Rcpp::XPtr<CGF_with_AD>>(currentList[0]);
//     Eigen::Index length = Rcpp::as<Eigen::Index>(currentList[1]);
//     cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(xp_cgf.get()));
//     lengths[i] = length;
//   }
// 
//   //// Pass these vectors into ConcatenationCGF
//   CGF_with_AD* CGF_base_ptr = new ConcatenationCGF(cgfs, lengths);
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
// 
//   attach_attributes(ptr, cgf_length_list);
//   return ptr;
// }
// 
// 
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_linearly_mappedCGF(Rcpp::XPtr<CGF_with_AD> base_cgf, mat Amat){
//   CGF_with_AD* CGF_ptr = new LinearlyMappedCGF(base_cgf, Amat);
//   Rcpp::XPtr<CGF_with_AD> ptr(CGF_ptr);
//   attach_attributes(ptr, base_cgf, Amat);
//   return ptr;
// }
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_randomly_stopped_sumCGF(Rcpp::XPtr<CGF_with_AD> count_base_cgf, Rcpp::XPtr<CGF_with_AD> summand_base_cgf){
//   CGF_with_AD* RSS_CGF = new RandomlyStoppedSumCGF(count_base_cgf, summand_base_cgf);
//   Rcpp::XPtr<CGF_with_AD> ptr(RSS_CGF);
//   attach_attributes(ptr, count_base_cgf, summand_base_cgf);
//   return ptr;
// }
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_scalar_to_iidCGF(Rcpp::XPtr<CGF_with_AD> scalar_base_cgf){
//   CGF_with_AD* Scalar_To_IID = new ScalarToIIDCGF(scalar_base_cgf);
//   Rcpp::XPtr<CGF_with_AD> ptr(Scalar_To_IID);
//   attach_attributes(ptr, scalar_base_cgf);
//   return ptr;
// }
// 
// // [[Rcpp::export]]
// Rcpp::XPtr<CGF_with_AD> make_iidReplicatesCGF(Rcpp::XPtr<CGF_with_AD> base_cgf, double block_size){
//   // base_cgf: A vector based CGF
//   // block_size: defines the size of the vector that is replicated. this is based on the observed data
//   CGF_with_AD* vector_to_iid_cgf = new IIDReplicates_CGF(base_cgf, block_size);
//   Rcpp::XPtr<CGF_with_AD> ptr(vector_to_iid_cgf);
//   attach_attributes(ptr, base_cgf);
//   return ptr;
// }
// 
// 
// 
// 
// 
// // Concatenation of Eigen - ADVectors
// // For a_vector vec1, vec2, vec3; // concatADVector({vec1, vec2, vec3}) returns a combined a_vector
// a_vector concatADVector(const std::vector<a_vector>& vecs) {
//   int totalSize = 0;
//   for (const auto& vec : vecs) { totalSize += vec.size();}
//   a_vector result(totalSize);
//   int currentIndex = 0;
//   for (const auto& vec : vecs) {
//     for (int i = 0; i < vec.size(); ++i) {
//       result[currentIndex++] = vec[i];
//     }
//   }
//   return result;
// }
// 
// template <class ModelCGFType>
// class K1ImplicitFunctionObject {
// private:
//   ModelCGFType* modelCGFInstance;
// 
// public:
//   K1ImplicitFunctionObject(ModelCGFType* instance) : modelCGFInstance(instance) {}
// 
//   // the method solve() is not directly in Matrix class but rather on the decomposition classes in Eigen.
//   // For dense-K2 something like K2.householderQr().solve(w) should work or even K2.llt().solve(w) for
//   // positive definite K2
//   template <class tvec_type, class theta_type, class w_type>
//   auto dfdu_solve(const tvec_type& tvec, const theta_type& theta_and_x, const w_type& w) {
//     Eigen::VectorXd tvec_converted = tvec;
//     Eigen::VectorXd theta_and_x_converted = theta_and_x;
//     Eigen::VectorXd theta_converted = theta_and_x_converted.head(theta_and_x.size()-tvec.size());
//     return modelCGFInstance->K2(tvec_converted, theta_converted).llt().solve(w).eval();
//   }
//   template <class tvec_type, class theta_type, class wT_type>
//   auto dfdu_solve_row(const tvec_type& tvec, const theta_type& theta_and_x, const wT_type& wT) {
//     // Since K2 is symmetric dfdu_solve and dfdu_solve_row are the same
//     return dfdu_solve(tvec, theta_and_x, wT);
//   }
// };
// 
// 
// 
// // This class acts as a callable interface to the provided CGF and the atomic_implicit_function class
// // It captures the CGF and provides a function-like interface for executing operations in the CGF
// template <class ModelCGFType>
// class ImplicitAtomicFunctionCGFWrapper {
// private:
//   ModelCGFType* modelCGFInstance;
//   atomic_implicit_function<K1ImplicitFunctionObject<CGF_with_AD>> aif;
// 
// public:
//   ImplicitAtomicFunctionCGFWrapper(ModelCGFType* instance) : modelCGFInstance(instance), aif("K1Implicit_Saddlepoint", instance) {}
// 
//   a_vector operator()(const a_vector& a_tvec, const a_vector& a_x, const a_vector& a_theta) {
//     a_vector a_c = modelCGFInstance->K1(a_tvec, a_theta) - a_x;
// 
//     a_vector a_input = concatADVector({a_tvec, a_theta, a_x, a_c});
//     a_vector a_tvec_hat(a_tvec.size());
//     aif(a_input, a_tvec_hat);
// 
//     return a_tvec_hat;
//   }
// };
// 
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// Rcpp::List computeFuncTGradient(vec tvec,
//                                 vec theta,
//                                 vec observations,
//                                 Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   a_vector a_theta = theta.cast<CppAD::AD<double>>();
//   a_vector a_tvec = tvec.cast<CppAD::AD<double>>();
//   a_vector a_x = observations.cast<CppAD::AD<double>>();
// 
//   // Create an instance of the ImplicitAtomicFunctionCGFWrapper using the provided modelCGF
//   ImplicitAtomicFunctionCGFWrapper<CGF_with_AD> implicit_func_cgfWrapper(modelCGF);
// 
//   CppAD::Independent(a_theta);
// 
//   // Use the cgfWrapper object to compute a_tvec_hat
//   a_vector a_tvec_hat = implicit_func_cgfWrapper(a_tvec, a_x, a_theta);
// 
//   a_vector a_funcT(1);
//   a_funcT(0) = modelCGF->func_T(a_tvec_hat, a_theta);
// 
//   CppAD::ADFun<double> ad_fun;
//   ad_fun.Dependent(a_theta, a_funcT);
// 
//   Rcpp::List result = Rcpp::List::create(Rcpp::Named("funcT") = ad_fun.Forward(0, theta),
//                                          Rcpp::Named("grad.theta.funcT") = ad_fun.Jacobian(theta)
//                                            // Rcpp::Named("grad.theta.t_hat") = negQB,
//   );
// 
//   return result;
// }
// 
// // [[Rcpp::export]]
// Rcpp::List computeZerothFuncTGradient(vec tvec,
//                                       vec theta,
//                                       vec observations,
//                                       Rcpp::XPtr<CGF_with_AD> modelCGF)
// {
//   a_vector a_theta = theta.cast<CppAD::AD<double>>();
//   a_vector a_tvec = tvec.cast<CppAD::AD<double>>();
//   a_vector a_x = observations.cast<CppAD::AD<double>>();
//   
//   // Create an instance of the ImplicitAtomicFunctionCGFWrapper using the provided modelCGF
//   ImplicitAtomicFunctionCGFWrapper<CGF_with_AD> implicit_func_cgfWrapper(modelCGF);
//   
//   CppAD::Independent(a_theta);
//   
//   // Use the cgfWrapper object to compute a_tvec_hat
//   a_vector a_tvec_hat = implicit_func_cgfWrapper(a_tvec, a_x, a_theta);
//   
//   a_vector a_funcT(1);
//   a_funcT(0) = -0.5*saddlepoint::atomic_funcs::logdet(modelCGF->K2(a_tvec_hat, a_theta));
//   
//   
//   CppAD::ADFun<double> ad_fun;
//   ad_fun.Dependent(a_theta, a_funcT);
//   
//   Rcpp::List result = Rcpp::List::create(Rcpp::Named("funcT") = ad_fun.Forward(0, theta),
//                                          Rcpp::Named("grad.theta.funcT") = ad_fun.Jacobian(theta)
//                                            // Rcpp::Named("grad.theta.t_hat") = negQB,
//   );
//   
//   return result;
// }


