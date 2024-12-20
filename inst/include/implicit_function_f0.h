#ifndef IMPLICIT_FUNCTION_FO_H
#define IMPLICIT_FUNCTION_FO_H

#include "saddlepoint_types.h"

using ConstDoubleVecMap = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>>;
using ConstADVecMap = Eigen::Map<const Eigen::Matrix<a_scalar, Eigen::Dynamic, 1>>;

          // using a_scalar = ad; // TMBad::ad_aug //
          // typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, 1> a_vector;
          // typedef Eigen::Matrix<a_scalar, Eigen::Dynamic, Eigen::Dynamic> a_matrix;


// This class defines an interface for an implicit function object
// that can be used with an implicit function atomic operator.
// It provides a generic form that we can inherit from and override.
//
// Requirements for the derived class:
// 1. Implement f(u,v): returns f(u,v) as an as a Eigen-based (ad) vector, given u and v.
// 2. Implement dfdu_solve(u,v,w): returns (df/du)^{-1} w, given u, v, and w. (we need both Types for this - ad/double)
//
// This base class also provides a helper method get_f_adfun(...) that constructs
// a TMBad::ADFun<> representing f(u,v) for a given size and initial points.
struct ImplicitFunctionFO {
  virtual ~ImplicitFunctionFO() {}

  // for AD types
  virtual a_vector f(const a_vector &u, const a_vector &v) const = 0;
  virtual matrix<ad> dfdu_solve(const ConstADVecMap &u, const ConstADVecMap &v, const ConstADVecMap &w) const = 0;
  
  // for double types
  virtual matrix<double> dfdu_solve(const ConstDoubleVecMap &u, const ConstDoubleVecMap &v, const ConstDoubleVecMap &w) const = 0;
  
  // get_f_adfun:
  // Given dimensions d (size of u), p (size of v),
  // and initial points u_init and v_init (as vectors of a_scalars),
  // this method builds and returns a TMBad::ADFun<> that computes f(u,v).
  //
  // This is useful when you need to compute the Jacobian of f(u,v) at a particular point,
  // as you can call ADFun.Jacobian(...) afterward.
  template <class Type>
  TMBad::ADFun<> get_f_adfun(const size_t d, const size_t p,
                             const std::vector<Type>& uv_init) const 
  {
    auto f_eval = [d, p, this](const std::vector<a_scalar>& X) {
      a_vector u_ad = Eigen::Map<const a_vector>(X.data(), d);
      a_vector v_ad = Eigen::Map<const a_vector>(X.data() + d, p);
      std::vector<a_scalar> res(d);
      Eigen::Map<a_vector>(res.data(), d) = this->f(u_ad, v_ad);
      return res;
    };
    
    return TMBad::ADFun<>(f_eval, uv_init);
  }
};

#endif // IMPLICIT_FUNCTION_FO_H