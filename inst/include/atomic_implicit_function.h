#ifndef ATOMIC_IMPLICIT_FUNCTION_H
#define ATOMIC_IMPLICIT_FUNCTION_H


# include "implicit_function_f0.h"


// Forward declaration
template<class dummy=void>
void implicit_function_internal(const CppAD::vector<TMBad::ad_aug>& x, 
                                CppAD::vector<TMBad::ad_aug>& y, 
                                const ImplicitFunctionFO* fobj_ptr);


template <class FunctionObjectType>
struct implicit_functionOp : TMBad::global::DynamicInputOutputOperator {
  typedef TMBad::global::DynamicInputOutputOperator Base;
  
  // The function object that provides:
  //   f(u,v)        : compute function f
  //   dfdu_solve(u,v,w) : solve (df/du)^{-1} w
  const FunctionObjectType* function_object;
  
  // Constructor: 
  // Here 'n' is the number of outputs (dimension of u)
  // and 'm' is the number of inputs.
  implicit_functionOp(TMBad::Index n, TMBad::Index m, const FunctionObjectType* fobj)
  : Base(n, m), function_object(fobj) {}
  
  const char* op_name() {
    return "implicit_function";
  }
  
  static const bool add_static_identifier = true;
  
  // Forward mode (no replay) - zeroth order
  // Simply copies inputs to outputs
  // We start from the assumption that output y = u (initial guess).
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {
    CppAD::vector<TMBad::Scalar> x_vals(this->input_size());
    CppAD::vector<TMBad::Scalar> y_vals(this->output_size());
    for (size_t i = 0; i < x_vals.size(); ++i) x_vals[i] = _args_.x(i);
    for (size_t i = 0; i < y_vals.size(); ++i) {
      // At zero-order: output y = input portion that represents u
      y_vals[i] = x_vals[i];
    }
    for (size_t i = 0; i < y_vals.size(); ++i) _args_.y(i) = y_vals[i];
  }
  
  // Forward mode with replay (for recording the operation)
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {
    CppAD::vector<TMBad::Replay> x_vals(this->input_size());
    CppAD::vector<TMBad::Replay> y_vals(this->output_size());
    for (size_t i = 0; i < x_vals.size(); ++i) x_vals[i] = _args_.x(i);
    implicit_function_internal(x_vals, y_vals, function_object);
    for (size_t i = 0; i < y_vals.size(); ++i) _args_.y(i) = y_vals[i];
  }
  
  // Reverse mode differentiation:
  // Compute partial derivatives with respect to the inputs, given partial derivatives with respect to the outputs.
  template <class Type>
  void reverse(TMBad::ReverseArgs<Type> _args_) {
    // The input vector is arranged as (u, v, fvals),
    // where u and fvals each have dimension d (the output size),
    // and v has dimension p such that input_size() = 2d + p.
    
    size_t d = this->output_size();       // dimension of u (and fvals)
    size_t p = this->input_size() - 2*d;  // dimension of v
    
    // Extract values from _args_:
    // u: first d entries of input
    // v: next p entries
    // fvals: last d entries (f(u,v))
    std::vector<Type> u(d), v(p), fvals(d), u_hat(d);
    for (size_t i = 0; i < d; ++i) u[i] = _args_.x(i);
    for (size_t i = 0; i < p; ++i) v[i] = _args_.x(d+i);
    for (size_t i = 0; i < d; ++i) fvals[i] = _args_.x(d+p+i);
    for (size_t i = 0; i < d; ++i) u_hat[i] = _args_.y(i);
    
    // Partial derivatives of final scalar function wrt output y:
    std::vector<Type> py(d);
    for (size_t i = 0; i < d; ++i) py[i] = _args_.dy(i);
    
    // We now compute w = - (df/du)^{-1} py at (u_hat, v)
    std::vector<Type> w(d);
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> EigenVector_type;
    Eigen::Map<const EigenVector_type> u_hat_eigen(u_hat.data(), d);
    Eigen::Map<const EigenVector_type> v_eigen(v.data(), p);
    Eigen::Map<const EigenVector_type> py_eigen(py.data(), d);
    Eigen::Map<EigenVector_type> w_eigen(w.data(), d);
    
    w_eigen = - function_object->dfdu_solve(u_hat_eigen, v_eigen, py_eigen);
    
    // Next, we need to compute the Jacobian of f(u_hat, v) wrt (u_hat, v),
    // and multiply by w to get partials wrt v and fvals.
    // First we obtain an ADFun for f(u,v)
    // Combine u_hat and v into one vector for ADFun:
    std::vector<Type> u_hat_and_v(d+p);
    for (size_t i = 0; i < d; ++i) u_hat_and_v[i] = u_hat[i];
    for (size_t i = 0; i < p; ++i) u_hat_and_v[d+i] = v[i];
    TMBad::ADFun<> f_adfun = function_object->get_f_adfun(d, p, u_hat_and_v);
    // Compute Jacobian * w:
    std::vector<Type> pw = f_adfun.Jacobian(u_hat_and_v, w);
    
    // pw represents (w * df/du_hat, w * df/dv).
    // For the implicit function, we set partial wrt u to 0,
    // partial wrt v is pw block after the first d entries,
    // and partial wrt fvals is -w.
    
    std::vector<Type> px(this->input_size());
    // px = (pu, pv, pfvals)
    // pu = 0 (no dependence on original u)
    for (size_t i = 0; i < d; ++i) px[i] = 0;
    // pv = the second block of pw (from d to d+p-1)
    for (size_t i = 0; i < p; ++i) px[d+i] = pw[d+i];
    // pfvals = -w
    for (size_t i = 0; i < d; ++i) px[d+p+i] = w[i];
    
    // Add results back into _args_:
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
void implicit_function_internal(const CppAD::vector<TMBad::ad_aug>& x, 
                                CppAD::vector<TMBad::ad_aug>& y, 
                                const ImplicitFunctionFO* fobj_ptr) 
{
  // x are the inputs: (u, v, f(u,v))
  // y are the outputs: u_hat
  
  TMBad::Index n = x.size();     // Total number of inputs
  TMBad::Index m = y.size();     // Number of outputs (dimension of u)
  
  // Instantiate operator:
  // Here, ImplicitFunctionFO is the type of your function object that provides:
  typedef implicit_functionOp<ImplicitFunctionFO> OP;
  
  // Create the operator instance from global stack
  // Then convert input from ad_aug to ad_plain for stack addition
  TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m, fobj_ptr);
  std::vector<TMBad::ad_plain> X(&x[0], &x[0] + x.size());
  std::vector<TMBad::ad_plain> Y = TMBad::get_glob()->add_to_stack<OP>(pOp, X);
  for (size_t i = 0; i < Y.size(); ++i) y[i] = Y[i]; 
}








#endif // ATOMIC_IMPLICIT_FUNCTION_H
