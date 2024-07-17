#ifndef TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED
#define TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED

// #include "baseCGF.h"
# include "saddlepoint_types.h"

template <class T>
  std::string display_name();
template <>
  std::string display_name<double>() {return "double";}
template <>
  std::string display_name<TMBad::ad_aug>() {return "ad_aug";}


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


template<class dummy=void>
CppAD::vector<TMBad::ad_aug> tvec_hat_internal(const CppAD::vector<TMBad::ad_aug>& x, const CGF_with_AD* cgf_ptr);

template<class dummy=void>
CppAD::vector<double> tvec_hat_internal (const CppAD::vector<double> &tx, CppAD::vector<TMBad::ad_aug>& ty, const CGF_with_AD* cgf_ptr);

template<class dummy>
CppAD::vector<TMBad::ad_aug> tvec_hat_internal(const CppAD::vector<TMBad::ad_aug> &x, const CGF_with_AD* cgf_ptr);

template<class dummy>
CppAD::vector<double> tvec_hat_internal (const CppAD::vector<double> &tx, const CGF_with_AD* cgf_ptr) {
  CppAD::vector<double> ty(CppAD::Integer(tx[0]));
  for (size_t i = 0; i < ty.size(); i++) ty[i] = tx[i+1];
  return ty;
}

a_vector tvec_hat(a_vector tvec, a_vector theta, a_vector observations, const CGF_with_AD* cgf_ptr);


template <class CGF_type>
struct tvec_hatOp : TMBad::global::DynamicInputOutputOperator {
  typedef TMBad::global::DynamicInputOutputOperator Base;
  
  const CGF_type* cgf_ptr_instance;
  K1ImplicitFunctionObject<CGF_with_AD> dfdu_solve_func;
  
  tvec_hatOp(TMBad::Index n, TMBad::Index m, const CGF_type* cgf_instance)
  : Base(n, m), cgf_ptr_instance(cgf_instance), dfdu_solve_func(cgf_instance) {}
  
  const char* op_name() {
    return "tvec_hat";
  }
  
  static const bool add_static_identifier = true;
  
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {
    CppAD::vector<TMBad::Scalar> tx(this->input_size());
    CppAD::vector<TMBad::Scalar> ty(this->output_size());
    
    // Rcout << "In forward: " << this->input_size() << " " << this->output_size() << std::endl;
    for (size_t i = 0; i < ty.size(); i++) ty[i] = _args_.x(i+1);
    // for (size_t i = 0; i < tx.size(); i++) tx[i] = _args_.x(i);
    // tvec_hat_internal(tx, ty, cgf_ptr_instance);
    for (size_t i = 0; i < ty.size(); i++) _args_.y(i) = ty[i];
  }
  
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {
    // Rcout << "In forward replay: " << this->input_size() << " " << this->output_size() << std::endl;
    
    CppAD::vector<TMBad::Replay> tx(this->input_size());
    for (size_t i = 0; i < tx.size(); i++) tx[i] = _args_.x(i);
    
    size_t d = this->output_size();
    size_t p = tx.size() - 3*d - 1;
    a_vector tvec(d);
    a_vector theta(p);
    a_vector observations(d);
    for (size_t i = 0; i < d; i++) tvec[i] = _args_.x(i+1);
    for (size_t i = 0; i < p; i++) theta[i] = _args_.x(i+1+d);
    for (size_t i = 0; i < d; i++) observations[i] = _args_.x(i+1+d+p);
    a_vector ty = tvec_hat(tvec, theta, observations, cgf_ptr_instance);
    for (int i = 0; i < ty.size(); i++) _args_.y(i) = ty(i);
    
    //CppAD::vector<TMBad::Replay> ty = tvec_hat_internal(tx, cgf_ptr_instance);
    
    //for (int i = 0; i < ty.size(); i++) _args_.y(i) = ty[i];
  }
  
  template <class Type>
  void reverse(TMBad::ReverseArgs<Type> _args_) {
    // Rcout << "reverse, Type = " << display_name<Type>() << std::endl;
    if (isDouble<Type>::value && this->output_size() == 1 && _args_.dy(0) == Type(0)) {
      return;
    }
    
    CppAD::vector<Type> tx(this->input_size());
    CppAD::vector<Type> ty(this->output_size());
    CppAD::vector<Type> px(this->input_size());
    CppAD::vector<Type> py(this->output_size());
    
    for (size_t i = 0; i < tx.size(); i++) tx[i] = _args_.x(i);
    for (size_t i = 0; i < ty.size(); i++) ty[i] = _args_.y(i);
    for (size_t i = 0; i < py.size(); i++) py[i] = _args_.dy(i);
    
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> EigenVector_type;
    
    size_t n = ty.size();
    size_t m = tx.size() - 3*n - 1;
    // Rcout << "In reverse: " << n << " " << m << " tx.size() = " << tx.size() << ", py.size() = " << py.size() << std::endl;
    Eigen::Map<EigenVector_type> partial_u(px.data()+1, n);
    partial_u = Eigen::Map<const EigenVector_type>(py.data(), n);
    Eigen::Map<EigenVector_type> partial_v(px.data()+n+1, m+n); // theta and x together
    partial_v.setConstant(0);
    Eigen::Map<EigenVector_type> partial_c(px.data()+2*n+m+1, n);
    
    //Eigen::Map<const EigenVector_type> mapped_u(ty.data(), n);
    Eigen::Map<const EigenVector_type> mapped_u(tx.data()+1, n); 
    // To check: take t from x or y? Apparently from x
    // TO DO: understand why this is correct
    
    Eigen::Map<const EigenVector_type> mapped_v(tx.data()+n+1, m); // just theta 
    partial_c = - dfdu_solve_func(mapped_u, mapped_v, Eigen::Map<const EigenVector_type>(py.data(), n));
    px[0] = 0;
    
    // Rcout << "reverse before for loop" << std::endl;
    for (size_t i=0; i<px.size(); i++) _args_.dx(i) += px[i];
    // Rcout << "reverse finished" << std::endl;
    // Rcout << "In reverse, values in px: ";
    // for (size_t i = 0; i < px.size(); i++) Rcout << px[i] << " ";
    // Rcout << std::endl;
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
CppAD::vector<TMBad::ad_aug> tvec_hat_internal(const CppAD::vector<TMBad::ad_aug> &tx, const CGF_with_AD* cgf_ptr) {
  // Rcout << "In tvec_hat_internal(tx, CGF) ad_aug version: " << tx.size() << std::endl;
  
  TMBad::Index n = tx.size();
  TMBad::Index m = CppAD::Integer(tx[0]) ;
  // Rcout << "In tvec_hat_internal(tx, CGF) ad_aug version: m = " << m << std::endl;
  
  typedef tvec_hatOp<CGF_with_AD> OP;
  bool all_constant = true;
  for (size_t i = 0; i < tx.size(); i++) all_constant &= tx[i].constant();
  
  CppAD::vector<TMBad::ad_aug> ty(m);
  
  if (all_constant) {
    for (size_t i = 0; i < m; i++) ty[i] = tx[i+1].Value();
  } else {
    TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m, cgf_ptr);
    std::vector<TMBad::ad_plain> x(&tx[0], &tx[0] + tx.size());
    std::vector<TMBad::ad_plain> y = TMBad::get_glob()->add_to_stack<OP>(pOp, x);
    for (size_t i = 0; i < y.size(); i++) ty[i] = y[i];
  }
  
  return ty;
}

template<class dummy=void>
void tvec_hat_internal(const CppAD::vector<TMBad::ad_aug>& tx, CppAD::vector<TMBad::ad_aug>& ty, const CGF_with_AD* cgf_ptr) {
  // Rcout << "In tvec_hat_internal(tx, ty, CGF) ad_aug version: " << tx.size() << " " << ty.size() << std::endl;
  // //ty = tvec_hat_internal(tx, cgf_ptr);
  
  TMBad::Index n = tx.size();
  TMBad::Index m = ty.size();
  
  typedef tvec_hatOp<CGF_with_AD> OP;
  bool all_constant = true;
  for (size_t i = 0; i < tx.size(); i++) all_constant &= tx[i].constant();
  
  if (all_constant) {
    for (size_t i = 0; i < m; i++) ty[i] = tx[i+1].Value();
  } else {
    TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m, cgf_ptr);
    std::vector<TMBad::ad_plain> x(&tx[0], &tx[0] + tx.size());
    std::vector<TMBad::ad_plain> y = TMBad::get_glob()->add_to_stack<OP>(pOp, x);
    for (size_t i = 0; i < y.size(); i++) ty[i] = y[i];
  }
}
template<class dummy=void>
void tvec_hat_internal(const CppAD::vector<double>& tx, CppAD::vector<double>& ty, const CGF_with_AD* cgf_ptr) {
  // // Rcout << "In tvec_hat_internal(tx, ty, CGF) double version: " << tx.size() << " " << ty.size() << std::endl;
  //ty = tvec_hat_internal(tx, cgf_ptr);
  
  // // TMBad::Index n = tx.size();
  // TMBad::Index m = ty.size();
  
  for (TMBad::Index i = 0; i < ty.size(); ++i) ty[i] = tx[i+1];
}

a_vector tvec_hat(a_vector tvec, a_vector theta, a_vector observations, const CGF_with_AD* cgf_ptr){
  a_vector a_c = cgf_ptr->K1(tvec, theta) - observations;
  
  CppAD::vector<a_scalar> arg(1+tvec.size()+theta.size()+observations.size()+a_c.size());
  arg[0] = tvec.size();
  for(int i=0;i<tvec.size();i++){arg[1+i]=tvec(i);}
  for(int i=0;i<theta.size();i++){arg[1+i+tvec.size()]=theta(i);}
  for(int i=0;i<observations.size();i++){arg[1+i+tvec.size()+theta.size()]=observations(i);}
  for(int i=0;i<a_c.size();i++){arg[1+i+tvec.size()+theta.size()+observations.size()]=a_c(i);}
  
  // Rcout << "In tvec_hat(tvec,theta,obs,CGF): " << tvec.size() << " " << theta.size() << " " << observations.size() << " " << a_c.size() << " " << arg.size() << " " << std::endl;
  CppAD::vector<a_scalar> res(tvec.size());
  tvec_hat_internal(arg, res, cgf_ptr);
  a_vector result(res.size());
  for(size_t i=0;i<res.size();i++){result(i)=res[i];}
  return result;
}


#endif // TVEC_ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED