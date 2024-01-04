#ifndef ADAPTORSWITHR_H_INCLUDED
#define ADAPTORSWITHR_H_INCLUDED

#include "parametric_submodelCGF.h"
//[[Rcpp::depends(RcppEigen)]]
# include <Rcpp.h>
# include <RcppEigen.h>

namespace saddlepoint {
namespace CGFs_with_Rcpp {

class AdaptorFromAtomicFunction
  : public CGFs_with_AD::Adaptor {

private:
  Rcpp::Function r_funcH;
  CppAD::atomic_three<double>* f_atom;

public:
  AdaptorFromAtomicFunction(Rcpp::Function r_H, CppAD::atomic_three<double>* f_atomic)
    : r_funcH(r_H), f_atom(f_atomic) {}
  // The user supplies a function that transforms the parameter vector
  // when this vector is supplied as a const vec&

  // The user also supplies an atomic function that transforms the parameter vector
  // when the parameter vector is supplied as const a_vector&
  // This atomic function should perform the same operations as the function
  // and must return a result of the same dimension as the function

  vec operator()(const vec& model_parameter_vector) const override {
    Eigen::Map<Eigen::VectorXd> func_Hvalues = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_funcH(model_parameter_vector));
    return func_Hvalues;
  }
  a_vector operator()(const a_vector& model_parameter_vector) const override {
    // Note: atomic functions do not create their output vector but rather write their results
    // into a vector that must already be allocated.
    // For this reason, the input vector is first converted to a non-AD equivalent.
    // Then r_funcH is invoked on this equivalent vector.
    // The result of this computation is needed only for its size, in order to allocate
    // the output a_vector of the correct size.
    // Then the atomic function result can be written into this vector.

    // Note: this is inefficient.
    // However, code operating on AD types should only execute during the recording of an
    // ADFun object, which should only be recorded once and then used repeatedly.
    // So this inefficiency should not make too much difference.

    size_t input_dim = model_parameter_vector.size();
    // Convert model_parameter_vector to non-AD equivalent
    vec mpv(input_dim);
    for (size_t i = 0; i < input_dim; ++i) {
      mpv[i] = CppAD::Value(CppAD::Var2Par(model_parameter_vector[i]));
    }
    Rcpp::NumericVector val0 = (r_funcH(mpv));
    size_t output_dim = val0.size();
    a_vector result(output_dim);
    // Write output of atomic function into result
    (*f_atom)(model_parameter_vector, result);
    return(result);
  }
};

class AdaptorFromAtomicFunctionWithRFunctions
  : public CGFs_with_AD::Adaptor {

private:
  CppAD::atomic_three<double>* f_atom;
  Rcpp::Function r_funcH;

public:
  AdaptorFromAtomicFunctionWithRFunctions(Rcpp::Function r_H, CppAD::atomic_three<double>* f_atomic)
    : r_funcH(r_H), f_atom(f_atomic) {}
  // The user supplies a function that transforms the parameter vector
  // when this vector is supplied as a const vec&

  // The user also supplies an atomic function that transforms the parameter vector
  // when the parameter vector is supplied as const a_vector&
  // This atomic function should perform the same operations as the function
  // and must return a result of the same dimension as the function

  vec operator()(const vec& model_parameter_vector) const override {
    Eigen::Map<Eigen::VectorXd> func_Hvalues = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_funcH(model_parameter_vector));
    return func_Hvalues;
  }
  a_vector operator()(const a_vector& model_parameter_vector) const override {
    // Note: atomic functions do not create their output vector but rather write their results
    // into a vector that must already be allocated.
    // For this reason, the input vector is first converted to a non-AD equivalent.
    // Then r_funcH is invoked on this equivalent vector.
    // The result of this computation is needed only for its size, in order to allocate
    // the output a_vector of the correct size.
    // Then the atomic function result can be written into this vector.

    // Note: this is inefficient.
    // However, code operating on AD types should only execute during the recording of an
    // ADFun object, which should only be recorded once and then used repeatedly.
    // So this inefficiency should not make too much difference.

    size_t input_dim = model_parameter_vector.size();
    // Convert model_parameter_vector to non-AD equivalent
    vec mpv(input_dim);
    for (size_t i = 0; i < input_dim; ++i) {
      mpv[i] = CppAD::Value(CppAD::Var2Par(model_parameter_vector[i]));
    }
    Rcpp::NumericVector val0 = (r_funcH(mpv));
    size_t output_dim = val0.size();
    a_vector result(output_dim);
    // Write output of atomic function into result
    (*f_atom)(model_parameter_vector, result);
    return(result);
  }
};



class AdaptorWithStoredData
  : public CGFs_with_AD::Adaptor {
public:
  // Convenience typedefs for pointers to functions
  typedef vec (*d_func_pointer)(const vec&, Rcpp::List);
  typedef a_vector (*ad_func_pointer)(const a_vector&, Rcpp::List);

private:
  Rcpp::List mylist;
  d_func_pointer f_d;
  ad_func_pointer f_ad;

public:
  AdaptorWithStoredData(d_func_pointer f_double,
                        ad_func_pointer f_ad_double,
                        Rcpp::List mylist0){
    f_d = f_double;
    f_ad = f_ad_double;
    mylist = mylist0;
  }
  // The user supplies two functions that transform the parameter vector
  // One function uses double as its scalar type, the other AD<double>

  vec operator()(const vec& model_parameter_vector) const override {
    return (*f_d)(model_parameter_vector, mylist);
  }
  a_vector operator()(const a_vector& model_parameter_vector) const override {
    return (*f_ad)(model_parameter_vector, mylist);
  }
};

} // namespace CGFs_with_Rcpp
} // namespace saddlepoint

#endif // ADAPTORSWITHR_H_INCLUDED
