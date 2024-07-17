#ifndef ADAPTORSWITHR_H_INCLUDED
#define ADAPTORSWITHR_H_INCLUDED

#include "parametric_submodelCGF.h"
// //[[Rcpp::depends(RcppEigen)]]
// # include <Rcpp.h>
// # include <RcppEigen.h>

typedef Rcpp::ComplexVector RADvector;

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


namespace saddlepoint {
namespace CGFs_with_Rcpp {



class AdaptorUsingRFunctions: public CGFs_with_AD::Adaptor {
  
private:
  Rcpp::Function r_function;
  
public:
  AdaptorUsingRFunctions(Rcpp::Function r_fun)
    : r_function(r_fun) {}
  // The user supplies a function that transforms the parameter vector
  // when this vector is supplied as a const vec&
  // The supplied function must be compatible with RTMB methods (that keep the AD class attribute)
  
  // double version
  vec operator()(const vec& model_parameter_vector) const override {
    Eigen::Map<Eigen::VectorXd> res = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_function(model_parameter_vector));
    return res;
  }
  // a_vector version
  a_vector operator()(const a_vector& model_parameter_vector) const override {
    RADvector f_advector = r_function(send_a_vector_to_advector(model_parameter_vector));
    return get_a_vector_from_advector(f_advector);
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
