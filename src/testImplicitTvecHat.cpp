
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
  
  // The function object that provides 'f(u,v)' and 'dfdu_solve(u,v,w)'.
  // In our example, that is a pointer to 'SaddlepointFO'.
  // We store it so we can call it in the reverse pass, etc.
  const FunctionObjectType* fobj_;
  
  // We also store the R function 'saddlepoint.solve' (numeric solver),
  // if needed. Or store an environment that can call it. 
  // For demonstration, I'll store it as an Rcpp::Function:
  Rcpp::Function saddlepoint_solve_r;
  
  // Possibly also store observations, etc., if your numeric solver needs them.
  Rcpp::NumericVector obs_;
  
  // ---- Constructor ----
  SaddlepointSolveOp(TMBad::Index p, 
                     TMBad::Index d, 
                     const FunctionObjectType* fobj)
    : Base(p, d), 
      fobj_(fobj)
  {
    // We store the R function 'saddlepoint.solve' for the forward pass.
    // We can also store the observations, etc., if needed.
    saddlepoint_solve_r = fobj->getSaddlepointSolveFn();
    obs_ = fobj->getObs();
  }
  
  const char* op_name() {
    return "SaddlepointSolve";
  }
  
  static const bool add_static_identifier = true;
  
  // Forward pass (0th order) with real doubles
  // ------------------------------------------
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {
    
    // We have p_ = dimension of param
    // The input vector x has length p_
    // The output vector y has length d_ (the tvec dimension).
    std::vector<double> param_dbl(this->input_size());
    for (size_t i = 0; i < this->input_size(); i++) {
      param_dbl[i] = _args_.x(i);
    }
    
    // Now call the R function saddlepoint.solve(...) to get tvec.
    Rcpp::NumericVector param_r = Rcpp::wrap(param_dbl);
    //  saddlepoint.solve(theta, obs, cgf, ...) -> returns numeric tvec
    Rcpp::NumericVector tvec_sol = saddlepoint_solve_r(param_r, obs_, fobj_->getCgfObj());
    
    // Copy result to _args_.y
    for (size_t i = 0; i < tvec_sol.size(); i++) {
      _args_.y(i) = tvec_sol[i];
    }
  }
  
  
  // -----------------------------------------------------
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {
    CppAD::vector<TMBad::Replay> x_vals(this->input_size());
    CppAD::vector<TMBad::Replay> y_vals(this->output_size());
    for (size_t i = 0; i < x_vals.size(); ++i) x_vals[i] = _args_.x(i);
    saddlepoint_solve_internal(x_vals, y_vals, fobj_);
    for (size_t i = 0; i < y_vals.size(); ++i) _args_.y(i) = y_vals[i];
  }
  
  
  
  
  
  // Reverse mode (derivatives wrt param)
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
    for (size_t i = 0; i < d; ++i) px[i] = 0;
    for (size_t i = 0; i < p; ++i) px[d+i] = pw[d+i];
    for (size_t i = 0; i < d; ++i) px[d+p+i] = w[i];
    
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
                                ) 
{
  TMBad::Index p = x.size();     
  TMBad::Index d = fobj_ptr->getObs().size();
  
  typedef SaddlepointSolveOp<SaddlepointFO> OP;
  
  // Create the operator instance from global stack
  // Then convert input from ad_aug to ad_plain for stack addition
  TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(p, d, fobj_ptr);
  std::vector<TMBad::ad_plain> X(&x[0], &x[0] + x.size());
  std::vector<TMBad::ad_plain> Y = TMBad::get_glob()->add_to_stack<OP>(pOp, X);
  for (size_t i = 0; i < Y.size(); ++i) y[i] = Y[i]; 
}











