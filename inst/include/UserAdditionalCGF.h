#ifndef USERADDITIONALCGF_H_INCLUDED
#define USERADDITIONALCGF_H_INCLUDED

#include "baseCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"

#include "CGF_Defaults.h"

using namespace saddlepoint::CGFs_with_AD;
using namespace saddlepoint::CGFs_via_templates;


template<typename T>
using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// typedef Eigen::VectorXd vec;
// typedef Eigen::MatrixXd mat;
// 
// typedef CppAD::AD<double> a_scalar;
// typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, 1> a_vector;
// typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, Eigen::Dynamic> a_matrix;


namespace saddlepoint {
namespace CGFs_with_AD {

template <class CGF_Class>
using UserClass = CGF_with_AD_from_template<CGF_Class>;

} // namespace CGFs_with_AD
} // namespace saddlepoint



// The makeUserClass function is an exported function that when called,
// creates a new CGF object from a template with any class as input.
// It returns Rcpp::Xptr
template<typename T>
Rcpp::XPtr<CGF_with_AD> makeUserClass() {
  CGF_with_AD* CGF_base_ptr = new UserClass<T>();
  Rcpp::XPtr<CGF_with_AD> ptr(CGF_base_ptr);
  return ptr;
}






// BaseVectorisedScalarCGF is a template class designed to provide a streamlined interface for users 
// who have their CGF classes in an iid vectorised form.
// The Derived class is expected to be a templated version of ScalarCGF with specific member functions like vectorisedK, vectorisedK1, etc.
template <typename Derived>
class BaseVectorisedScalarCGF : public ScalarVectorisedIID_Defaults<Derived> {
private:
  // This private method getInstance() manages a static instance of the Derived class.
  // It ensures that only one instance of Derived is created and used throughout the lifetime of the program.
  // singleton design pattern.
  static Derived& getInstance() {
    static Derived instance;
    return instance;
  }
  
public:
  // The following methods are static and provide an interface to the vectorised versions of the K functions.
  // K_vectorised_iid() calls the vectorisedK() method of the Derived class for a given vector and parameters.
  template<class t_vector_type, class param_type>
  static auto K_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().vectorisedK(tvec, params);
  }
  
  
  template<class t_vector_type, class param_type>
  static auto K1_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().vectorisedK1(tvec, params);
  }
  
  
  template<class t_vector_type, class param_type>
  static auto K2_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().vectorisedK2(tvec, params);
  }
  
  
  template<class t_vector_type, class param_type>
  static auto K3_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().vectorisedK3(tvec, params);
  }
  
  
  template<class t_vector_type, class param_type>
  static auto K4_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().vectorisedK4(tvec, params);
  }
  
  template<class t_vector_type, class param_type>
  static auto ineq_constraint_vectorised_iid(const t_vector_type &tvec, const param_type &params) {
    return getInstance().VectorisedIneqConstraint(tvec, params);
  }
  
};

#endif // USERADDITIONALCGF_H_INCLUDED
