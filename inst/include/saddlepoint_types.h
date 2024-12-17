#ifndef SADDLEPOINT_TYPES_H_INCLUDED
#define SADDLEPOINT_TYPES_H_INCLUDED

// [[Rcpp::depends(TMB)]]
#include <RTMB.h>

// [[Rcpp::depends(RcppEigen)]]
# include <RcppEigen.h>

using a_scalar = ad;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, 1> a_vector;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, Eigen::Dynamic> a_matrix;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vec;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat;

namespace saddlepoint{
namespace CGFs_with_AD{
class Adaptor;
class CGF_with_AD;
class parametric_submodelCGF;
}

namespace CGFs_via_virtual_functions{
}
namespace CGFs_via_templates{}
}

using namespace saddlepoint::CGFs_with_AD;
// using namespace saddlepoint::CGFs_with_Rcpp;
// using namespace saddlepoint::atomic_funcs;
// //using saddlepoint::CGFs_via_virtual_functions::CGF_base;

template<class T, class... ArgTypes>
void attach_attributes(Rcpp::XPtr<T>& ptr, const ArgTypes&... args) {
  // This method is used in Rcpp-exported functions to help each function keep track of all the pointers it needs internally
  Rcpp::List L = Rcpp::List::create(args...);
  ptr.attr("internal_resources") = L;
}





#endif // SADDLEPOINT_TYPES_H_INCLUDED
