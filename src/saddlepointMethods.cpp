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
  std::copy(tvec.data(), tvec.data()+tvec.size(), combined_input.begin());
  std::copy(theta.data(), theta.data()+theta.size(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr< TMBad::ADFun<> > makeADFunNegll(const vec& tvec,
                                            const vec& theta,
                                            Rcpp::XPtr<CGF_with_AD> cgf){
  try{
      const CGF_with_AD* cgf_ptr = cgf.get();
      auto func = [cgf_ptr, tvec_size=tvec.size(), theta_size=theta.size()](const std::vector<a_scalar>& x) {
        a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), tvec_size);
        a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + tvec_size, theta_size);
        a_scalar result = cgf_ptr->neg_ll(tvec_ad, theta_ad);
        return std::vector<a_scalar>{result};
      };
      
      std::vector<double> combined_input(tvec.size() + theta.size());
      std::copy(tvec.data(), tvec.data()+tvec.size(), combined_input.begin());
      std::copy(theta.data(), theta.data()+theta.size(), combined_input.begin() + tvec.size());
      
      Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
      attach_attributes(ptr, cgf);
      return ptr;
  } catch (...) {
    Rcpp::stop("An unknown error occurred while creating ADFun object.");
  }
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
  std::copy(tvec.data(), tvec.data()+tvec.size(), combined_input.begin());
  std::copy(theta.data(), theta.data()+theta.size(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}


// [[Rcpp::export]]
Rcpp::XPtr< TMBad::ADFun<> > makeADFunZerothNegll(const vec& tvec,
                                                  const vec& theta,
                                                  Rcpp::XPtr<CGF_with_AD> cgf){
  const CGF_with_AD* cgf_ptr = cgf.get();
  auto func = [cgf_ptr, tvec_size=tvec.size(), theta_size=theta.size()](const std::vector<a_scalar>& x) {
    a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), tvec_size);
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + tvec_size, theta_size);
    a_scalar result = -cgf_ptr->tilting_exponent(tvec_ad, theta_ad); // return negative log-likelihood
    return std::vector<a_scalar>{result};
  };
  
  std::vector<double> combined_input(tvec.size() + theta.size());
  std::copy(tvec.data(), tvec.data()+tvec.size(), combined_input.begin());
  std::copy(theta.data(), theta.data()+theta.size(), combined_input.begin() + tvec.size());
  
  Rcpp::XPtr<TMBad::ADFun<>> ptr(new TMBad::ADFun<>(func, combined_input), true);
  attach_attributes(ptr, cgf);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List computeCombinedGradient(const vec& combined_vector, Rcpp::XPtr<TMBad::ADFun<>> adf){
  std::vector<double> combined_input(combined_vector.size());
  std::copy(combined_vector.data(), combined_vector.data()+combined_vector.size(), combined_input.begin());
  
  // Rcpp::Named("gradient") = adf->Jacobian(combined_input)
  
  return Rcpp::List::create(Rcpp::Named("objective") = (*adf)(combined_input),
                            Rcpp::Named("gradient") = (*adf).Jacobian(combined_input)
                            );
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
  // TMBad::ADFun<> jac_adf = adf.JacFun();

  return Rcpp::List::create(Rcpp::Named("value") = adf(theta_std),
                            Rcpp::Named("gradient") = adf.Jacobian(theta_std)
                            // Rcpp::Named("hessian") = jac_adf.Jacobian(theta_std)
                            );

}

// [[Rcpp::export]]
Rcpp::List computeZerothFuncT(const vec& tvec,
                              const vec& theta,
                              const vec& observations,
                              Rcpp::XPtr<CGF_with_AD> cgf)
{
  a_vector tvec_ad = tvec.cast<a_scalar>();
  a_vector observations_ad = observations.cast<a_scalar>();
  const CGF_with_AD* cgf_ptr = cgf.get();
  
  auto func = [cgf_ptr, tvec_ad, observations_ad](const std::vector<a_scalar>& x) {
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data(), x.size());
    a_vector tvec_ad_hat = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_ptr);
    matrix<a_scalar> K2_val = cgf_ptr->K2(tvec_ad_hat, theta_ad);
    return std::vector<a_scalar>{-0.5*atomic::logdet(K2_val)};
  };
  
  std::vector<double> theta_std(theta.data(), theta.data() + theta.size());
  TMBad::ADFun<> adf(func, theta_std);
  // TMBad::ADFun<> jac_adf = adf.JacFun();
  
  return Rcpp::List::create(Rcpp::Named("value") = adf(theta_std),
                            Rcpp::Named("gradient") = adf.Jacobian(theta_std)
                            // Rcpp::Named("hessian") = jac_adf.Jacobian(theta_std)
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


// [[Rcpp::export]]
Rcpp::List computeZerothNegll(const vec& tvec,
                              const vec& theta,
                              const vec& observations,
                              Rcpp::XPtr<CGF_with_AD> cgf){
  a_vector tvec_ad = tvec.cast<a_scalar>();
  a_vector observations_ad = observations.cast<a_scalar>();
  const CGF_with_AD* cgf_ptr = cgf.get();
  
  auto func = [cgf_ptr, tvec_ad, observations_ad](const std::vector<a_scalar>& x) {
    a_vector theta_ad = Eigen::Map<const a_vector>(x.data(), x.size());
    a_vector tvec_ad_hat = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_ptr);
    a_scalar result = -cgf_ptr->tilting_exponent(tvec_ad_hat, theta_ad);
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









