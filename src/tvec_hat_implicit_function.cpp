
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




// [[Rcpp::export]]
ADrep tvec_hat_from_tvec(ADrep theta, 
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
















































// We define a new function-object
struct SaddlepointFO  {
private:
  Rcpp::Function K1_r;           // R function for f(tvec, theta) = K1(...)
  Rcpp::Function dfdu_solve_r;   // R function for (df/dtvec)^{-1} * w
  Rcpp::Function saddlepoint_solve_r;  
  Rcpp::RObject cgf_obj;        
  Rcpp::NumericVector obs_;      // or whatever your data is
  
public:
  // Constructor
  SaddlepointFO(Rcpp::Function K1_func,
                Rcpp::Function dfdu_solve_func,
                Rcpp::Function saddlepoint_solve_func,
                Rcpp::RObject cgf_object,
                Rcpp::NumericVector obs)
    : K1_r(K1_func), 
      dfdu_solve_r(dfdu_solve_func), 
      saddlepoint_solve_r(saddlepoint_solve_func),
      cgf_obj(cgf_object), 
      obs_(obs)
  {}
  
  const Rcpp::NumericVector& getObs() const {return obs_;}
  const Rcpp::RObject& getCgfObj() const {return cgf_obj;}
  const Rcpp::Function& getSaddlepointSolveFn() const {return saddlepoint_solve_r;}
  
  
  //*************************
  //   f(tvec, theta) as AD
  //*************************
  a_vector f(const a_vector& tvec, const a_vector& theta) const {
    ADrep res_r(K1_r(send_a_vector_to_advector(tvec), 
                     send_a_vector_to_advector(theta)));
    RADvector advector_result = Rcpp::as<RADvector>(res_r);
    return get_a_vector_from_advector(advector_result);
  }
  
  //****************************
  //   df/du_solve(tvec,theta,w) 
  //****************************
  matrix<double> dfdu_solve(const ConstDoubleVecMap& tvec,
                            const ConstDoubleVecMap& theta,
                            const ConstDoubleVecMap& w) const  {
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
                        const ConstADVecMap& w) const {
    
    
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
  template <class Type>
  TMBad::ADFun<> get_f_adfun(const size_t d, const size_t p,
                             const std::vector<Type>& uv_init) const 
  {
    auto f_eval = [d, p, this](const std::vector<a_scalar>& X) {
      a_vector u_ad = Eigen::Map<const a_vector>(X.data(), d);
      a_vector v_ad = Eigen::Map<const a_vector>(X.data() + d, p);
      std::vector<a_scalar> res(d);
      Eigen::Map<a_vector>(res.data(), d) = this->f(u_ad, v_ad);
      return res;
    };
    
    return TMBad::ADFun<>(f_eval, uv_init);
  }
};






template<class dummy=void>
void saddlepoint_solve_internal(const CppAD::vector<TMBad::ad_aug>& x, 
                                CppAD::vector<TMBad::ad_aug>& y, 
                                const SaddlepointFO* fobj_ptr);


template <class FunctionObjectType>
struct SaddlepointSolveOp : TMBad::global::DynamicInputOutputOperator {
  using Base = TMBad::global::DynamicInputOutputOperator;
  
  const FunctionObjectType* fobj_;
  
  Rcpp::Function saddlepoint_solve_r;
  
  Rcpp::NumericVector obs_;
  
  // ---- Constructor ----
  SaddlepointSolveOp(TMBad::Index p, 
                     TMBad::Index d, 
                     const FunctionObjectType* fobj)
    : Base(p, d), 
      fobj_(fobj),
      saddlepoint_solve_r(fobj_->getSaddlepointSolveFn()), 
      obs_(fobj_->getObs())    {}
  
  const char* op_name() {
    return "SaddlepointSolve";
  }
  
  static const bool add_static_identifier = true;
  
  // Forward pass (0th order) with real doubles
  // ------------------------------------------
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {
    
    std::vector<double> x_vals(this->input_size());
    std::vector<double> param_dbl(this->output_size());
    for (size_t i = 0; i < x_vals.size(); ++i) x_vals[i] = _args_.x(i);
    for (size_t i = 0; i < param_dbl.size(); ++i) param_dbl[i] = x_vals[i];
    
    // Now call the R function saddlepoint.solve(...) to get tvec.
    Rcpp::NumericVector param_r = Rcpp::wrap(param_dbl);
    //  saddlepoint.solve(theta, obs, cgf, ...) -> returns numeric tvec
    Rcpp::NumericVector tvec_sol = saddlepoint_solve_r(param_r, obs_, fobj_->getCgfObj());
    
    // Copy result to _args_.y
    for (size_t i = 0; i < tvec_sol.size(); i++) _args_.y(i) = tvec_sol[i];
  }
  
  
  // -----------------------------------------------------
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {
    CppAD::vector<TMBad::Replay> x_vals(this->input_size());
    CppAD::vector<TMBad::Replay> y_vals(this->output_size());
    for (size_t i = 0; i < x_vals.size(); ++i) x_vals[i] = _args_.x(i);
    saddlepoint_solve_internal(x_vals, y_vals, fobj_);
    for (size_t i = 0; i < y_vals.size(); ++i) _args_.y(i) = y_vals[i];
  }
  
  
  // ------------------------------------
  template <class Type>
  void reverse(TMBad::ReverseArgs<Type> _args_) {
    size_t d = this->output_size();     
    size_t p = this->input_size() - d; 
    
    std::vector<Type> v(p), fvals(d), u_hat(d);
    for (size_t i = 0; i < p; ++i) v[i] = _args_.x(i);
    for (size_t i = 0; i < d; ++i) fvals[i] = _args_.x(p+i);
    for (size_t i = 0; i < d; ++i) u_hat[i] = _args_.y(i);
    
    std::vector<Type> py(d);
    for (size_t i = 0; i < d; ++i) py[i] = _args_.dy(i);
    
    std::vector<Type> w(d);
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> EigenVector_type;
    Eigen::Map<const EigenVector_type> u_hat_eigen(u_hat.data(), d);
    Eigen::Map<const EigenVector_type> v_eigen(v.data(), p);
    Eigen::Map<const EigenVector_type> py_eigen(py.data(), d);
    Eigen::Map<EigenVector_type> w_eigen(w.data(), d);
    
    w_eigen = - fobj_->dfdu_solve(u_hat_eigen, v_eigen, py_eigen);
    
    std::vector<Type> u_hat_and_v(d+p);
    for (size_t i = 0; i < d; ++i) u_hat_and_v[i] = u_hat[i];
    for (size_t i = 0; i < p; ++i) u_hat_and_v[d+i] = v[i];
    TMBad::ADFun<> f_adfun = fobj_->get_f_adfun(d, p, u_hat_and_v);
    std::vector<Type> pw = f_adfun.Jacobian(u_hat_and_v, w);
    
    
    std::vector<Type> px(this->input_size());
    for (size_t i = 0; i < p; ++i) px[i] = pw[i];
    for (size_t i = 0; i < d; ++i) px[p+i] = w[i];
    
    
    for (size_t i=0; i<px.size(); ++i) _args_.dx(i) += px[i];
  }
  
  void forward(TMBad::ForwardArgs<TMBad::Writer>& args) {
    if (!(false)) {
      Rcerr << "TMBad assertion failed.\n";
      Rcerr << "The following condition was not met: " << "false" << "\n";
      Rcerr << "Possible reason: " << "Unknown" << "\n";
      Rcerr << "For more info run your program through a debugger.\n";
      Rcpp::stop("TMB unexpected");
    }
  }
  void reverse(TMBad::ReverseArgs<TMBad::Writer>& args) {
    if (!(false)) {
      Rcerr << "TMBad assertion failed.\n";
      Rcerr << "The following condition was not met: " << "false" << "\n";
      Rcerr << "Possible reason: " << "Unknown" << "\n";
      Rcerr << "For more info run your program through a debugger.\n";
      Rcpp::stop("TMB unexpected");
    }
  }
};






template<class dummy>
void saddlepoint_solve_internal(const CppAD::vector<TMBad::ad_aug>& x, 
                                CppAD::vector<TMBad::ad_aug>& y, 
                                const SaddlepointFO* fobj_ptr
) {
  TMBad::Index p = x.size();     
  TMBad::Index d = fobj_ptr->getObs().size();
  
  typedef SaddlepointSolveOp<SaddlepointFO> OP;
  
  TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(p, d, fobj_ptr);
  std::vector<TMBad::ad_plain> X(&x[0], &x[0] + x.size());
  std::vector<TMBad::ad_plain> Y = TMBad::get_glob()->add_to_stack<OP>(pOp, X);
  for (size_t i = 0; i < Y.size(); ++i) y[i] = Y[i]; 
}


template<class Type>
vector<Type> saddlepoint_solve_hat(const vector<Type>& theta,
                                   const vector<Type>& observations,
                                   Rcpp::Function K1,
                                   Rcpp::Function dfdu_solve,
                                   Rcpp::Function saddlepoint_solve,
                                   Rcpp::RObject cgf_obj,
                                   Rcpp::NumericVector obs
){
  SaddlepointFO* fobj = new SaddlepointFO(K1, dfdu_solve, saddlepoint_solve, cgf_obj, obs);
  
  CppAD::vector<Type> arg(theta.size()+observations.size());
  for(int i=0;i<theta.size();++i){arg[i]=theta(i);}
  for(int i=0;i<observations.size();++i){arg[i+theta.size()]=observations(i);}
  
  CppAD::vector<Type> res(observations.size());
  saddlepoint_solve_internal(arg, res, fobj);
  vector<Type> result(res.size());
  for(size_t i=0;i<res.size();++i){result[i]=res[i];}
  return result;
}







// [[Rcpp::export]]
ADrep tapedSaddlepointSolve(ADrep theta,
                            vec observations,
                            SEXP K2_solve_fn,
                            SEXP saddlepoint_solve_fn,
                            SEXP cgf_obj
){
  
  // Rcpp::Function K1_r(K1_fn);
  Rcpp::Function K2_solve_r(K2_solve_fn);
  Rcpp::Function saddlepoint_solve_r(saddlepoint_solve_fn);

  Rcpp::RObject cgf_obj_r(cgf_obj);
  
  Rcpp::Environment cgf_env(cgf_obj_r);
  Rcpp::Function K1_r = cgf_env["K1"];
  
  const a_scalar* theta_ptr = adptr(theta);
  // const a_scalar* tvec_ptr = adptr(tvec);
  // const a_scalar* observations_ptr = adptr(observations);
  
  // Convert ADrep (R objects) to vector<a_scalar> (C++ objects).
  // We copy the elements from the pointers into our C++ vectors.
  vector<a_scalar> theta_ad(theta.size()), observations_ad(observations.size());
  // for (size_t i=0; i<tvec.size(); i++) tvec_ad(i)=tvec_ptr[i];
  // for (size_t i=0; i<observations.size(); i++) observations_ad(i)=observations_ptr[i];
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(theta.size()); i++) theta_ad(i)=theta_ptr[i];
  for (Eigen::Index i=0; i<observations.size(); i++) observations_ad(i)=observations(i);
  
  
  // Call the internal function which uses R's saddlepoint.solve
  // to compute the t-vector while capturing the AD dependency.
  vector<a_scalar> result = saddlepoint_solve_hat(theta_ad, 
                                                  observations_ad, 
                                                  K1_r, 
                                                  K2_solve_r, 
                                                  saddlepoint_solve_r, 
                                                  cgf_obj_r, 
                                                  Rcpp::wrap(observations));
  ADrep out_res(result.size());
  a_scalar* out_ptr = adptr(out_res);
  for (Eigen::Index i=0; i<static_cast<Eigen::Index>(result.size()); i++) out_ptr[i] = result[i];
  return out_res;
}
