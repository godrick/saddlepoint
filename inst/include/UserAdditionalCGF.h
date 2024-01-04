#ifndef USERADDITIONALCGF_H_INCLUDED
#define USERADDITIONALCGF_H_INCLUDED

#include "baseCGF.h"
#include "CGF_Defaults.h"

#include "parametric_submodelCGF.h"

template<typename T>
using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

typedef CppAD::AD<double> a_scalar;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, 1> a_vector;
typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, Eigen::Dynamic> a_matrix;

using namespace saddlepoint::CGFs_with_AD;
using namespace saddlepoint::CGFs_via_templates;



namespace saddlepoint {
namespace CGFs_with_AD {
// typedef for the user

template <class CGF_Class>
using UserClass = CGF_with_AD_from_template<CGF_Class>;


//template <class CGF_Class>
//using CGFWithAdditionalMatrix = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::AdaptedParametersCallbackCGF<CGF_Class, ProvideConstantMat>>;


}

}// namespace saddlepoint

// The makeUserClass function is an exported function that when called,
// creates a new CGF object from a template with any class as input.
// It returns Rcpp::Xptr
template<typename T>
Rcpp::XPtr<CGF_with_AD> makeUserClass() {
   CGF_with_AD* CGF_base_ptr = new UserClass<T>();
   Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
   return ptr;
}



#endif // USERADDITIONALCGF_H_INCLUDED
