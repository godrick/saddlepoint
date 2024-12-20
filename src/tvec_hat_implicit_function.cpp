
#include "atomic_implicit_function.h"


typedef Rcpp::ComplexVector RADvector;
typedef Rcpp::ComplexMatrix RADmatrix; 

RADvector send_a_vector_to_advector(const a_vector ax) {
  RADvector rax(ax.size());
  for(int i = 0; i < ax.size(); i++) rax[i] = ad2cplx(ax[i]);
  return as_advector(rax);
}
a_vector get_a_vector_from_advector(const RADvector rax) {
  a_vector avec_x(rax.size());
  for (int i = 0; i < rax.size(); ++i) avec_x[i] = cplx2ad(rax[i]);
  return avec_x;
}




// This class represents your specific implicit function object for tvec_hat.
// It inherits from ImplicitFunctionFO and uses an R6 CGF object with methods K1 and K2.
struct TvecHatImplicitFunctionFO : public ImplicitFunctionFO {
private:
  // Store a reference to the R6 CGF object environment.
  // Assume you have something like: my_r6_obj = new.env(); with methods my_r6_obj$K1 and my_r6_obj$K2 in R.
  Rcpp::Environment r6_obj_env;
  Rcpp::Function K1_r;  
  Rcpp::Function K2_r;  
  
public:
  // Constructor takes the R environment (or however you access the R6 object)
  TvecHatImplicitFunctionFO(const Rcpp::Environment& env) 
  : r6_obj_env(env), 
    K1_r(env["K1"]),  
    K2_r(env["K2"]) {}
  
  // f(u,v) calling K1
  a_vector f(const a_vector& tvec, const a_vector& theta) const override {
    ADrep res_r(K1_r(send_a_vector_to_advector(tvec), 
                     send_a_vector_to_advector(theta)));
    RADvector advector_result = Rcpp::as<RADvector>(res_r);
    return get_a_vector_from_advector(advector_result);
  }
  
  
  
  
  
  matrix<double> dfdu_solve(const ConstDoubleVecMap& tvec,
                            const ConstDoubleVecMap& theta,
                            const ConstDoubleVecMap& w) const override {

    Eigen::Matrix<double, Eigen::Dynamic, 1> tvec_converted = tvec;
    Eigen::Matrix<double, Eigen::Dynamic, 1> theta_converted = theta;

    SEXP K2_val_r = K2_r(tvec_converted, theta_converted);
    Rcpp::NumericMatrix K2_val(K2_val_r);

    matrix<double> K2_val_converted(K2_val.nrow(), K2_val.ncol());
    for (int i=0; i<K2_val.nrow(); i++)
      for (int j=0; j<K2_val.ncol(); j++)
        K2_val_converted(i,j) = K2_val(i,j);

    matrix<double> w_converted = w;

    return atomic::matmul(atomic::matinv(K2_val_converted), w_converted);
  }


  matrix<ad> dfdu_solve(const ConstADVecMap& tvec,
                        const ConstADVecMap& theta,
                        const ConstADVecMap& w) const override {

    Eigen::Matrix<ad, Eigen::Dynamic, 1> tvec_converted = tvec;
    Eigen::Matrix<ad, Eigen::Dynamic, 1> theta_converted = theta;

    ADrep k2_val_ad(K2_r(send_a_vector_to_advector(tvec_converted),
                         send_a_vector_to_advector(theta_converted)) );

    if (k2_val_ad.size() == 1) {
      // Ensuring single-element ADrep is handled as a 1x1 matrix:
      // If K2_r function returns a single-element ADrep (i.e., size() == 1),
      // it is initially treated as a scalar (vector without explicit dimensions).
      // Rcpp::ComplexMatrix/MatrixInput function requires dimension attributes to treat the data as a matrix.
      // Thus, we manually set the dimensions for single-element ADrep to be treated as a 1x1 matrix.
      Rcpp::IntegerVector new_dim = Rcpp::IntegerVector::create(1, 1);
      Rf_setAttrib(k2_val_ad, R_DimSymbol, new_dim);
    }

    ConstMapMatrix K2_val_admat = MatrixInput(k2_val_ad);
    matrix<ad> w_converted = w;

    return atomic::matmul(atomic::matinv(matrix<ad>(K2_val_admat)), w_converted);
  }
  
};
















template<class Type>
vector<Type> tvec_hat(const vector<Type>& tvec, 
                      const vector<Type>& theta, 
                      const vector<Type>& observations, 
                      const Rcpp::Environment& cgf_env){
  
  // Create an instance of the function object
  TvecHatImplicitFunctionFO* fobj = new TvecHatImplicitFunctionFO(cgf_env);
  
  CppAD::vector<Type> arg(tvec.size()+theta.size()+observations.size());
  for(int i=0;i<tvec.size();++i){arg[i]=tvec(i);}
  for(int i=0;i<theta.size();++i){arg[i+tvec.size()]=theta(i);}
  for(int i=0;i<observations.size();++i){arg[i+tvec.size()+theta.size()]=observations(i);}
  
  CppAD::vector<Type> res(tvec.size());
  implicit_function_internal(arg, res, fobj);
  vector<Type> result(res.size());
  for(size_t i=0;i<res.size();++i){result[i]=res[i];}
  return result;
}


//**** At the time of writing this ADrep class is still not available in the master branch of RTMB.
// //**** Remember to revert to Rcpp::ComplexVector if ADrep remains unavailable. 
// [[Rcpp::export]]
ADrep tvec_hat(ADrep tvec, ADrep theta, ADrep observations, SEXP cgf_env){
  
  // Extract raw pointers to the underlying array of 'ad' elements.
  // 'adptr(x)' returns a pointer (ad*) to the first element of x.
  // We can then treat these pointers as arrays and index them.
  const a_scalar* tvec_ptr = adptr(tvec);
  const a_scalar* theta_ptr = adptr(theta);
  const a_scalar* observations_ptr = adptr(observations);
  
  
  // Convert ADrep (R objects) to vector<a_scalar> (C++ objects).
  // We copy the elements from the pointers into our C++ vectors.
  vector<a_scalar> tvec_ad(tvec.size()), theta_ad(theta.size()), observations_ad(observations.size());
  for (size_t i=0; i<tvec.size(); i++) tvec_ad(i)=tvec_ptr[i];
  for (size_t i=0; i<theta.size(); i++) theta_ad(i)=theta_ptr[i];
  for (size_t i=0; i<observations.size(); i++) observations_ad(i)=observations_ptr[i];
  
  // Call the main tvec_hat function that works with vector<a_scalar>.
  vector<a_scalar> result = tvec_hat(tvec_ad, theta_ad, observations_ad, cgf_env);
  
  ADrep out_res(result.size());
  a_scalar* out_ptr = adptr(out_res);
  for (size_t i=0; i<result.size(); i++) out_ptr[i] = result[i];
  return out_res;
}

