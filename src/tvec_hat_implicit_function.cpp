
// To streamline interactions between R and C++, all matrix operations are handled in R.
// Previously, sparse matrices caused issues, but with RTMB's current support for them,
// it's advantageous to manage all functions on the R side.
// In C++, we work exclusively with vectors.


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



// Function object implementing the ImplicitFunctionFO interface.
// Interfaces with R functions K1 and dfdu_solve for computations.
struct TvecHatImplicitFunctionFO : public ImplicitFunctionFO {
private:
  Rcpp::Function K1_r;         ///< R function for computing f(u, v)
  Rcpp::Function dfdu_solve_r; ///< R function for solving (df/du)^{-1} w
  
public:
  // Constructor that initializes R function callbacks.
  TvecHatImplicitFunctionFO(Rcpp::Function K1_func, Rcpp::Function dfdu_solve_func)
    : K1_r(K1_func), dfdu_solve_r(dfdu_solve_func) {}
  
  // f(u,v) calling K1
  // Overrides the f(u, v) method by invoking the R K1 function.
  a_vector f(const a_vector& tvec, const a_vector& theta) const override {
    ADrep res_r(K1_r(send_a_vector_to_advector(tvec), 
                     send_a_vector_to_advector(theta)));
    RADvector advector_result = Rcpp::as<RADvector>(res_r);
    return get_a_vector_from_advector(advector_result);
  }
  
  // dfdu_solve implemented via R function
  // Overrides the dfdu_solve method for double types by invoking the R dfdu_solve function.
  matrix<double> dfdu_solve(const ConstDoubleVecMap& tvec,
                            const ConstDoubleVecMap& theta,
                            const ConstDoubleVecMap& w) const override {
                              // Convert inputs to R-compatible formats
                              Rcpp::NumericVector tvec_r = Rcpp::wrap(tvec);
                              Rcpp::NumericVector theta_r = Rcpp::wrap(theta);
                              Rcpp::NumericVector w_r = Rcpp::wrap(w);
                              
                              
                              // Call the dfdu_solve R function
                              Rcpp::NumericVector solution_r = dfdu_solve_r(tvec_r, theta_r, w_r);
                              
                              
                              // Convert the solution back to matrix<double>
                              Eigen::VectorXd solution_eigen = Rcpp::as<Eigen::VectorXd>(solution_r);
                              matrix<double> res = solution_eigen;
                              return res;
                            }
  
  matrix<ad> dfdu_solve(const ConstADVecMap& tvec,
                        const ConstADVecMap& theta,
                        const ConstADVecMap& w) const override {
                          
                          
                          Eigen::Matrix<ad, Eigen::Dynamic, 1> tvec_converted = tvec;
                          Eigen::Matrix<ad, Eigen::Dynamic, 1> theta_converted = theta;
                          Eigen::Matrix<ad, Eigen::Dynamic, 1> w_converted = w;
                          
                          // Convert inputs to R-compatible formats
                          RADvector tvec_r = send_a_vector_to_advector(tvec_converted);
                          RADvector theta_r = send_a_vector_to_advector(theta_converted);
                          RADvector w_r = send_a_vector_to_advector(w_converted);
                          
                          
                          // Call the dfdu_solve R function
                          RADvector solution_r = dfdu_solve_r(tvec_r, theta_r, w_r);
                          
                          // Convert the solution back to matrix<ad>
                          a_vector solution = get_a_vector_from_advector(solution_r);
                          matrix<ad> res = solution;
                          return res;
                        }
};


template<class Type>
vector<Type> tvec_hat(const vector<Type>& tvec, 
                      const vector<Type>& theta, 
                      const vector<Type>& observations, 
                      Rcpp::Function K1,         
                      Rcpp::Function dfdu_solve        
){
  // Create an instance of the function object
  TvecHatImplicitFunctionFO* fobj = new TvecHatImplicitFunctionFO(K1, dfdu_solve);
  
  // std::unique_ptr<TvecHatImplicitFunctionFO> fobj(new TvecHatImplicitFunctionFO(K1, dfdu_solve));
  
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



//**** At the time of writing this, ADrep class is still not available in the master branch of RTMB.
// //**** Remember to revert to Rcpp::ComplexVector if ADrep remains unavailable. 
// [[Rcpp::export]]
ADrep tvec_hat_for_ad(ADrep theta, 
                      vec tvec,
                      vec observations,
                      SEXP K1_fn,    
                      SEXP K2_solve_fn 
){
  
  Rcpp::Function K1_r(K1_fn);
  Rcpp::Function K2_solve_r(K2_solve_fn);
  
  // Extract raw pointers to the underlying array of 'ad' elements.
  // 'adptr(x)' returns a pointer (ad*) to the first element of x.
  // We can then treat these pointers as arrays and index them.
  const a_scalar* theta_ptr = adptr(theta);
              // const a_scalar* tvec_ptr = adptr(tvec);
              // const a_scalar* observations_ptr = adptr(observations);
  
  
  
  // Convert ADrep (R objects) to vector<a_scalar> (C++ objects).
  // We copy the elements from the pointers into our C++ vectors.
  vector<a_scalar> tvec_ad(tvec.size()), theta_ad(theta.size()), observations_ad(observations.size());
  // for (size_t i=0; i<tvec.size(); i++) tvec_ad(i)=tvec_ptr[i];
  // for (size_t i=0; i<observations.size(); i++) observations_ad(i)=observations_ptr[i];
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(theta.size()); i++) theta_ad(i)=theta_ptr[i];
  for (Eigen::Index i=0; i<tvec.size(); i++) tvec_ad(i)=tvec(i);
  for (Eigen::Index i=0; i<observations.size(); i++) observations_ad(i)=observations(i);
  

  
  // Call the main tvec_hat function that works with vector<a_scalar>.
  vector<a_scalar> result = tvec_hat(tvec_ad, theta_ad, observations_ad, K1_r, K2_solve_r);
  ADrep out_res(result.size());
  a_scalar* out_ptr = adptr(out_res);
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(result.size()); i++) out_ptr[i] = result[i];
  return out_res;
}


