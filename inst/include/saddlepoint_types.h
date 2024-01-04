#ifndef SADDLEPOINT_TYPES_H_INCLUDED
#define SADDLEPOINT_TYPES_H_INCLUDED

// [[Rcpp::plugins("cpp17")]]
// [[Rcpp::depends(RcppEigen)]]
# include <Rcpp.h>
# include <RcppEigen.h>

#include <cppad/cppad.hpp>


typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

typedef CppAD::AD<double> a_scalar;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, 1> a_vector;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, Eigen::Dynamic> a_matrix;

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
