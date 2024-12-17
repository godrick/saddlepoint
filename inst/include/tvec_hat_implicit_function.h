#ifndef TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED
#define TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED

// #include "baseCGF.h"
# include "saddlepoint_types.h"


template <class CGF_type>
class K1ImplicitFunctionObject {
private:
  const CGF_type* modelCGFInstance;
  
public:
  K1ImplicitFunctionObject(const CGF_type* instance) : modelCGFInstance(instance) {}
  
  template <class tvec_type, class theta_type, class w_type>
  auto operator()(const tvec_type& tvec, const theta_type& theta, const w_type& w) const {
    typedef decltype(tvec.sum()) scalar_type;
    
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> tvec_converted = tvec;
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> theta_converted = theta;
    
    matrix<scalar_type> K2_val = modelCGFInstance->K2(tvec_converted, theta_converted);
    matrix<scalar_type> w_converted = w;
    return atomic::matmul(atomic::matinv(K2_val), w_converted);
    
  }
};




// Forward declaration
template<class dummy=void>
void tvec_hat_internal(const CppAD::vector<TMBad::ad_aug>& tx, CppAD::vector<TMBad::ad_aug>& ty, const CGF_with_AD* cgf_ptr);


template <class CGF_type>
struct tvec_hatOp : TMBad::global::DynamicInputOutputOperator {
  typedef TMBad::global::DynamicInputOutputOperator Base;
  
  const CGF_type* cgf_ptr_instance;
  K1ImplicitFunctionObject<CGF_type> dfdu_solve_func;
  
  tvec_hatOp(TMBad::Index n, TMBad::Index m, const CGF_type* cgf_instance)
    : Base(n, m), cgf_ptr_instance(cgf_instance), dfdu_solve_func(cgf_instance) {}
  
  const char* op_name() {
    return "tvec_hat";
  }
  
  static const bool add_static_identifier = true;
  
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {
    CppAD::vector<TMBad::Scalar> tx(this->input_size());
    CppAD::vector<TMBad::Scalar> ty(this->output_size());
    
    for (size_t i = 0; i < tx.size(); ++i) tx[i] = _args_.x(i);
    for (size_t i = 0; i < ty.size(); ++i) ty[i] = tx[i];
    for (size_t i = 0; i < ty.size(); ++i) _args_.y(i) = ty[i];
  }
  
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {
    CppAD::vector<TMBad::Replay> tx(this->input_size());
    CppAD::vector<TMBad::Replay> ty(this->output_size());
    for (size_t i = 0; i < tx.size(); ++i) tx[i] = _args_.x(i);
    tvec_hat_internal(tx, ty, cgf_ptr_instance);
    for (size_t i = 0; i < ty.size(); ++i) _args_.y(i) = ty[i];
  }
  
  template <class Type>
  void reverse(TMBad::ReverseArgs<Type> _args_) {
    //    if (isDouble<Type>::value && this->output_size() == 1 && _args_.dy(0) == Type(0)) {
    //      return;
    //    }
    // Skip this to avoid strange and unlikely special cases in debugging?
    
    // Note: x is packed as (tvec, theta, obs)
    
    // this->input_size() and this->output_size() should in principle be available from sizes of x and y
    // but not clear how this information is accessed
    // TO DO: determine from workings of TMBad::ReverseArgs<Type>
    size_t d = this->output_size();
    size_t p = this->input_size() - 2*d;
    
    // Copy already-computed values
    std::vector<Type> tvec(d);
    std::vector<Type> theta(p);
    std::vector<Type> obs(d);
    std::vector<Type> tvec_hat(d);
    
    // From x:
    for (size_t i = 0; i < d; ++i) tvec[i] = _args_.x(i);
    for (size_t i = 0; i < p; ++i) theta[i] = _args_.x(d+i);
    for (size_t i = 0; i < d; ++i) obs[i] = _args_.x(d+p+i);
    // From y:
    for (size_t i = 0; i < d; ++i) tvec_hat[i] = _args_.y(i);
    
    // Copy supplied partial_y values:
    std::vector<Type> py(this->output_size());
    for (size_t i = 0; i < py.size(); ++i) py[i] = _args_.dy(i);
    
    // Initialize partial_x to be computed
    std::vector<Type> px(this->input_size());
    // Notionally, px = partial_x consists of (ptvec, ptheta, pobs)
    // We will set these values at the end of reverse()
    
    // Next compute py * (dK1/dt)^{-1}, to be stored as std::vector for later use
    // (where we understand all p... vectors as row vectors)
    // Use tvec_hat values, not tvec
    std::vector<Type> w(d);
    // Use Eigen mappings
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> EigenVector_type;
    Eigen::Map<EigenVector_type> mapped_w(w.data(), d);
    Eigen::Map<const EigenVector_type> mapped_tvec_hat(tvec_hat.data(), tvec.size());
    Eigen::Map<const EigenVector_type> mapped_theta(theta.data(), theta.size());
    Eigen::Map<const EigenVector_type> mapped_py(py.data(), py.size());
    // Assign to mapped w
    mapped_w = - dfdu_solve_func(mapped_tvec_hat, mapped_theta, mapped_py);
    
    // Next create ADFun object for K1
    // Prepare combined argument vector
    std::vector<Type> tvec_hat_and_theta(d+p);
    // Store tvec_hat into combined vector **note: not tvec**
    for (size_t i = 0; i < d; ++i) tvec_hat_and_theta[i] = tvec_hat[i];
    // followed by theta
    for (size_t i = 0; i < p; ++i) tvec_hat_and_theta[d+i] = theta[i];
    
    auto K1func = [d, p, this](const std::vector<a_scalar>& x) {
      a_vector tvec_ad = Eigen::Map<const a_vector>(x.data(), d);
      a_vector theta_ad = Eigen::Map<const a_vector>(x.data() + d, p);
      std::vector<a_scalar> res(d);
      Eigen::Map<a_vector>(res.data(), d) = cgf_ptr_instance->K1(tvec_ad, theta_ad);
      return res;
    };
    TMBad::ADFun<> K1_adf(K1func, tvec_hat_and_theta);
    std::vector<Type> pt_and_ptheta = K1_adf.Jacobian(tvec_hat_and_theta, w);
    // We understand pt_and_ptheta as the row vector (w * dK1/dt, w * dK1/dtheta) in block form
    // ptheta will be set using the second block
    // The first block will be discarded and ptvec will be set to 0
    
    // Set values in px
    // ptvec == 0
    for (size_t i = 0; i < 1+d; ++i) px[i] = 0;
    // ptheta == second block of pt_and_ptheta 
    // (entries d,d+1,...,d+p-1 of pt_and_ptheta)
    for (size_t i = d; i < d+p; ++i) px[i] = pt_and_ptheta[i];
    // pobs == -w
    for (size_t i = 0; i < d; ++i) px[d+p+i] = w[i];
    
    // Results copied into _args_.dx
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
void tvec_hat_internal(const CppAD::vector<TMBad::ad_aug>& tx, CppAD::vector<TMBad::ad_aug>& ty, const CGF_with_AD* cgf_ptr) {
  //ty = tvec_hat_internal(tx, cgf_ptr);
  // Don't do this; x-only form of tvec_hat_internal is more complicated
  // Instead: mostly copying body, with ty.size() instead of tx[0]
  
  TMBad::Index n = tx.size();
  TMBad::Index m = ty.size();
  
  // For simplicity skip to not-all-constant branch
  typedef tvec_hatOp<CGF_with_AD> OP;
  TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m, cgf_ptr);
  std::vector<TMBad::ad_plain> x(&tx[0], &tx[0] + tx.size());
  std::vector<TMBad::ad_plain> y = TMBad::get_glob()->add_to_stack<OP>(pOp, x);
  for (size_t i = 0; i < y.size(); ++i) ty[i] = y[i];
}


a_vector tvec_hat(a_vector tvec, a_vector theta, a_vector observations, const CGF_with_AD* cgf_ptr){
  CppAD::vector<a_scalar> arg(tvec.size()+theta.size()+observations.size());
  // arg[0] = tvec.size();
  for(int i=0;i<tvec.size();++i){arg[i]=tvec(i);}
  for(int i=0;i<theta.size();++i){arg[i+tvec.size()]=theta(i);}
  for(int i=0;i<observations.size();++i){arg[i+tvec.size()+theta.size()]=observations(i);}
  // for(int i=0;i<a_c.size();++i){arg[1+i+tvec.size()+theta.size()+observations.size()]=a_c(i);}
  
  CppAD::vector<a_scalar> res(tvec.size());
  tvec_hat_internal(arg, res, cgf_ptr);
  a_vector result(res.size());
  for(size_t i=0;i<res.size();++i){result(i)=res[i];}
  return result;
}


#endif // TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED
