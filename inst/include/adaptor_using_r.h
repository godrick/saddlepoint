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
RADvector numeric2RAD(const Rcpp::NumericVector &x) {
  RADvector ans(x.size());
  for (int i=0; i<x.size(); i++) ans[i] = ad2cplx(ad(x[i]));
  return as_advector(ans);
}

RADvector SEXP2RADvector(const SEXP &x) {
  // TYPEOF(x) - REALSXP/14, INTSXP/13, CPLXSXP/15(RADvector), STRSXP/16, VECSXP/19(list), LGLSXP/10
  // Rcpp::Rcout << "TYPEOF(x) = " << TYPEOF(x) << " (" << Rf_type2char(TYPEOF(x)) << ")" << std::endl;
  if (TYPEOF(x) == REALSXP) return numeric2RAD(Rcpp::as<Rcpp::NumericVector>(x));
  if (TYPEOF(x) == INTSXP){
    Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(x);
    Rcpp::NumericVector nv(iv.begin(), iv.end());
    return numeric2RAD(nv);
  } 
  return Rcpp::as<RADvector>(x);
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
    SEXP r_function_result = r_function(model_parameter_vector);
    Rcpp::NumericVector numeric_result;
    if (TYPEOF(r_function_result) == INTSXP) { // from IntegerVector to NumericVector
      Rcpp::IntegerVector int_result = Rcpp::as<Rcpp::IntegerVector>(r_function_result);
      numeric_result = Rcpp::NumericVector(int_result.begin(), int_result.end());
    } 
    else if (TYPEOF(r_function_result) == REALSXP) {
      numeric_result = Rcpp::as<Rcpp::NumericVector>(r_function_result);
    } 
    else {
      Rcpp::stop("The result of the adaptor function is not of a supported type");
    }
    Eigen::Map<Eigen::VectorXd> res(numeric_result.begin(), numeric_result.size());
    return res;
  }
  // a_vector version for automatic differentiation
  a_vector operator()(const a_vector& model_parameter_vector) const override {

    SEXP r_function_result = r_function(send_a_vector_to_advector(model_parameter_vector));

    RADvector advector_result = SEXP2RADvector(r_function_result);

    return get_a_vector_from_advector(advector_result);
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
