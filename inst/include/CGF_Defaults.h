#ifndef CGF_DEFAULTS_H_INCLUDED
#define CGF_DEFAULTS_H_INCLUDED

#define _USE_MATH_DEFINES
#include <utility>
#include <cmath>
#include <Eigen/Eigen>

// #include "saddlepoint_types.h"

namespace saddlepoint {
namespace CGFs_via_templates {

#include "needs_full_list_of_CGF_methods.h"
template <class BasicCGF>
class CGF_Defaults {
public:
  // Provides default behaviour for certain CGF methods, when these methods can be implemented in terms of methods from BasicCGF
  // BasicCGF must provide methods K, K1, K2, K3operator, K4operator

  // Usage: compile-time polymorphism/curiously recurring template pattern
  // template <class T> class newCGF : public CGF_defaults<newCGF<T>> { /* class definition */ };
  // Any method not implemented in newCGF directly will call the version in CGF_Defaults
  // Do not use CGF_Defaults in any other way, as this will invalidate the static_cast calls used by CGF_Defaults

  template <class t_vector_type, class x_type, class y_type, class... ParameterTypes>
  auto K2operator(const t_vector_type& tvec, const x_type& x, const y_type& y, ParameterTypes&&... params) const
    //        -> decltype(x.transpose() * static_cast<const BasicCGF*>(this)->K2(tvec, std::forward<ParameterTypes>(params)...) * y)
  {
    return (x.transpose() * static_cast<const BasicCGF*>(this)->K2(tvec, std::forward<ParameterTypes>(params)...) * y).value();
  }

  template <class t_vector_type, class... ParameterTypes, class matrix_type>
  auto K2operatorAK2AT(const t_vector_type& tvec, const matrix_type& A, ParameterTypes&&... params) const
    //        -> decltype(A * static_cast<const BasicCGF*>(this)->K2(tvec, std::forward<ParameterTypes>(params)...) * A.transpose())
  {
    return (A * static_cast<const BasicCGF*>(this)->K2(tvec, std::forward<ParameterTypes>(params)...) * A.transpose()).eval();
  }

  template <class t_vector_type, class... ParameterTypes>
  auto tilting_exponent(const t_vector_type& tvec, const ParameterTypes&... params) const
    //        -> decltype(static_cast<const BasicCGF*>(this)->K(tvec, params...) - ( tvec.array() * static_cast<const BasicCGF*>(this)->K1(tvec, params...).array() ).sum())
  {
    return static_cast<const BasicCGF*>(this)->K(tvec, params...)
    - ( tvec.array() * static_cast<const BasicCGF*>(this)->K1(tvec, params...).array() ).sum();
    // Note: parameters passed as const reference, rather than being forwarded via std::forward, since they are used more than once
    // TO DO: Check whether this is appropriate
  }

  template <class t_vector_type, class... ParameterTypes>
  auto neg_ll(const t_vector_type& tvec, const ParameterTypes&... params) const  {
    typedef decltype(tvec.sum()) scalar_type;
    return 0.5 * atomic::logdet(matrix<scalar_type>(static_cast<const BasicCGF*>(this)->K2(tvec, params...)))
    + 0.5*log(2*M_PI)*tvec.size() - tilting_exponent(tvec, params...);
  }

  template <class t_vector_type, class... ParameterTypes>
  auto func_T(const t_vector_type& tvec, const ParameterTypes&... params) const
    // decltype for compatibility only
    // In effect it replicates the effect of the code, but is unwieldy to look at
    //        -> decltype(K4operatorAABB_factored(tvec, decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), params...)/8 - K3K3operatorAABBCC_factored(tvec, decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), params...)/8 - K3K3operatorABCABC_factored(tvec, decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), decltype( static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().reconstructedMatrix() )::Identity(tvec.size(), tvec.size()) * static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().transpositionsP().transpose(), static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt().vectorD(), params...)/12 )
  {

    // Returns the first higher-order correction term in the saddlepoint expansion

    // Specify Q = K2^{-1} by computing an LDLT decomposition for K2
    // and using it to compute the A and d arguments to the _factored form of the operator methods
    // To save repeated inversion, this is performed once rather than delegated to the _factored_inverse methods
    // Note that with K2 = P^T L D L^T P, we have Q = K2^{-1} = P^{-1} L^{-1}^T D^{-1} L^{-1} P^{-1}^T
    // and since P is a permutation matrix P^T = P^{-1} and this reduces to
    // Q = (L^{-1} P)^T D^{-1} L^{-1} P, i.e., A=(L^{-1}P)^T = P^T (L^T)^{-1}

    // auto k2ldlt = static_cast<const BasicCGF*>(this)->K2(tvec, params...).ldlt();
    //
    // auto dim = tvec.size();
    // typedef decltype( k2ldlt.reconstructedMatrix() ) matrix_type;
    // matrix_type A = matrix_type::Identity(dim, dim) * k2ldlt.transpositionsP().transpose();
    // k2ldlt.matrixU().template solveInPlace<Eigen::OnTheRight>(A); // matrixU() provides U=L^T and solveInPlace changes A into A U^{-1}


    auto inverse_k2ldlt = static_cast<const BasicCGF*>(this)->K2(tvec, params...).inverse().ldlt();
    typedef decltype( inverse_k2ldlt.reconstructedMatrix() ) matrix_type;
    matrix_type inverse_k2ldlt_L = inverse_k2ldlt.matrixL();
    matrix_type A = inverse_k2ldlt.transpositionsP().transpose() * inverse_k2ldlt_L;

    auto d = inverse_k2ldlt.vectorD();

    auto K4operatorAABB_val = K4operatorAABB_factored(tvec, A, d, A, d, params...);
    auto K3K3operatorAABBCC_val = K3K3operatorAABBCC_factored(tvec, A, d, A, d, A, d, params...);
    auto K3K3operatorABCABC_val = K3K3operatorABCABC_factored(tvec, A, d, A, d, A, d, params...);

    // auto inverse_k2 = static_cast<const BasicCGF*>(this)->K2(tvec, params...).inverse();
    //
    // auto K4operatorAABB_val = K4operatorAABB(tvec, inverse_k2, inverse_k2, params...);
    // auto K3K3operatorAABBCC_val = K3K3operatorAABBCC(tvec, inverse_k2, inverse_k2,inverse_k2, params...);
    // auto K3K3operatorABCABC_val = K3K3operatorABCABC(tvec,inverse_k2,inverse_k2,inverse_k2, params...);

    return K4operatorAABB_val/8 - K3K3operatorAABBCC_val/8 - K3K3operatorABCABC_val/12;
  }

  template <class t_vector_type, class... ParameterTypes>
  auto ineq_constraint(const t_vector_type& tvec, const ParameterTypes&... params) const
  {
    typedef typename std::remove_cv_t<typename std::remove_reference_t< decltype(tvec.eval()) > > vector_type;

    vector_type res;
    return res;
  }


  template <class t_vector_type, class Q_matrix_type, class... ParameterTypes>
  auto K4operatorAABB(const t_vector_type& tvec, const Q_matrix_type& Q1, const Q_matrix_type& Q2, const ParameterTypes&... params) const
    //        -> decltype( K4operatorAABB_factored(tvec, Q1.ldlt().matrixL().eval(), Q1.ldlt().vectorD(), Q1.ldlt().matrixL().eval(), Q2.ldlt().vectorD(), params...) )
  {
    auto ldlt1 = Q1.ldlt();
    auto ldlt2 = Q2.ldlt();

    typedef decltype(ldlt1.reconstructedMatrix()) matrix_type;
    auto extract_A = [](const decltype(ldlt1)& ldlt) {
      matrix_type res = ldlt.matrixL();
      res = ldlt.transpositionsP().transpose() * res;
      return res;
    };

    return K4operatorAABB_factored(tvec, extract_A(ldlt1), ldlt1.vectorD(), extract_A(ldlt2), ldlt2.vectorD(), params...);
  }

  template <class t_vector_type, class Q_matrix_type, class... ParameterTypes>
  auto K3K3operatorAABBCC(const t_vector_type& tvec, const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                          const ParameterTypes&... params) const
    //        -> decltype( K3K3operatorAABBCC_factored(tvec, Q1.ldlt().matrixL().eval(), Q1.ldlt().vectorD(), Q1.ldlt().matrixL().eval(), Q2.ldlt().vectorD(), Q1.ldlt().matrixL().eval(), Q3.ldlt().vectorD(), params...) )
  {
    auto ldlt1 = Q1.ldlt();
    auto ldlt2 = Q2.ldlt();
    auto ldlt3 = Q3.ldlt();

    typedef decltype(ldlt1.reconstructedMatrix()) matrix_type;
    auto extract_A = [](const decltype(ldlt1)& ldlt) {
      matrix_type res = ldlt.matrixL();
      res = ldlt.transpositionsP().transpose() * res;
      return res;
    };
    return K3K3operatorAABBCC_factored(tvec,
                                       extract_A(ldlt1), ldlt1.vectorD(),
                                       extract_A(ldlt2), ldlt2.vectorD(),
                                       extract_A(ldlt3), ldlt3.vectorD(),
                                       params...);
  }

  template <class t_vector_type, class Q_matrix_type, class... ParameterTypes>
  auto K3K3operatorABCABC(const t_vector_type& tvec, const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                          const ParameterTypes&... params) const
    //        -> decltype( K3K3operatorABCABC_factored(tvec, Q1.ldlt().matrixL().eval(), Q1.ldlt().vectorD(), Q1.ldlt().matrixL().eval(), Q2.ldlt().vectorD(), Q1.ldlt().matrixL().eval(), Q3.ldlt().vectorD(), params...) )
  {
    auto ldlt1 = Q1.ldlt();
    auto ldlt2 = Q2.ldlt();
    auto ldlt3 = Q3.ldlt();

    typedef decltype(ldlt1.reconstructedMatrix()) matrix_type;
    auto extract_A = [](const decltype(ldlt1)& ldlt) {
      matrix_type res = ldlt.matrixL();
      res = ldlt.transpositionsP().transpose() * res;
      return res;
    };
    return K3K3operatorABCABC_factored(tvec,
                                       extract_A(ldlt1), ldlt1.vectorD(),
                                       extract_A(ldlt2), ldlt2.vectorD(),
                                       extract_A(ldlt3), ldlt3.vectorD(),
                                       params...);
  }

  template <class t_vector_type, class matrix_type, class d_vector_type, class... ParameterTypes>
  auto K4operatorAABB_factored(const t_vector_type& tvec,
                               const matrix_type& A1, const d_vector_type& d1,
                               const matrix_type& A2, const d_vector_type& d2,
                               const ParameterTypes&... params) const
    //        -> decltype( d1(0)*d2(0)*static_cast<const BasicCGF*>(this)->K4operator(tvec, A1.col(0), A1.col(0), A2.col(0), A2.col(0), params...) )
  {
    // Key identity:
    // sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
    // = sum_{m1=1,...,r1} sum_{m2=1,...,r2} d1(m1)*d2(m2)*( sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)*A2(i4,m2) )
    // Complexity: r1*r2 calls to K4operator
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();

    typedef decltype( d1(0)*d2(0)*static_cast<const BasicCGF*>(this)->K4operator(tvec, A1.col(0), A1.col(0), A2.col(0), A2.col(0), params...) )
      scalar_type;
    scalar_type res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for (index_type m2 = 0; m2 < r2; ++m2) {
        res += d1(m1) * d2(m2) * static_cast<const BasicCGF*>(this)->K4operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), A2.col(m2), params...);
      }
    }
    return res;
  }

  template <class t_vector_type, class matrix_type, class d_vector_type, class... ParameterTypes>
  auto K3K3operatorAABBCC_factored(const t_vector_type& tvec,
                                   const matrix_type& A1, const d_vector_type& d1,
                                   const matrix_type& A2, const d_vector_type& d2,
                                   const matrix_type& A3, const d_vector_type& d3,
                                   const ParameterTypes&... params) const
    //        -> decltype( static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(0), A1.col(0), A2.col(0), params...) )
  {
    // Key identity:
    // sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
    // = sum_{m2=1,...,r2} d2(m2)*( sum_{m1=1,...,r1} d1(m1)*(sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)) )
    //                           *( sum_{m3=1,...,r3} d3(m3)*(sum_{i4,...i6=1,...,d} K3(i4,i5,i6)*A2(i4,m2)*A3(i5,m3)*A3(i6,m3)) )
    // Complexity: r2*(r1+r3) calls to K3operator
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();

    typedef decltype( static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(0), A1.col(0), A2.col(0), params...) )
      scalar_type;
    scalar_type res = 0;
    for(index_type m2 = 0; m2 < r2; ++m2) {
      scalar_type factor1 = 0;
      for(index_type m1 = 0; m1 < r1; ++m1) {
        factor1 += d1(m1) * static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), params...);
      }
      scalar_type factor2 = 0;
      for(index_type m3 = 0; m3 < r3; ++m3) {
        factor2 += d3(m3) * static_cast<const BasicCGF*>(this)->K3operator(tvec, A2.col(m2), A3.col(m3), A3.col(m3), params...);
      }
      res += d2(m2)*factor1*factor2;
    }
    return res;
  }

  template <class t_vector_type, class matrix_type, class d_vector_type, class... ParameterTypes>
  auto K3K3operatorABCABC_factored(const t_vector_type& tvec,
                                   const matrix_type& A1, const d_vector_type& d1,
                                   const matrix_type& A2, const d_vector_type& d2,
                                   const matrix_type& A3, const d_vector_type& d3,
                                   const ParameterTypes&... params) const
    //        -> decltype( static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(0), A2.col(0), A3.col(0), params...) )
  {
    // Key identity:
    // sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
    // = sum_{m1=1,...,r1} sum_{m2=1,...,r2} sum_{m3=1,...,r3} d1(m1)*d2(m2)*d3(m3)*
    //          ( sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A2(i2,m2)*A3(i3,m3) )^2
    // since the sum involving i4,i5,i6 is identical to the sum involving i1,i2,i3 except for the labelling of the variables of summation
    // Complexity: r1*r2*r3 calls to K3operator
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();

    typedef decltype( static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(0), A2.col(0), A3.col(0), params...) )
      scalar_type;
    scalar_type res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for(index_type m2 = 0; m2 < r2; ++m2) {
        for(index_type m3 = 0; m3 < r3; ++m3) {
          auto square = [](scalar_type x) {return x*x;};
          res += d1(m1)*d2(m2)*d3(m3)*square(static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(m1), A2.col(m2), A3.col(m3), params...));
        }
      }
    }
    return res;
  }
//   template <class t_vector_type, class matrix_type, class d_vector_type, class... ParameterTypes>
//   auto K3K3operatorABCABC_factored(const t_vector_type& tvec,
//                                    const matrix_type& A1, const d_vector_type& d1,
//                                    const matrix_type& A2, const d_vector_type& d2,
//                                    const matrix_type& A3, const d_vector_type& d3,
//                                    const ParameterTypes&... params) const
//   {
//     typedef decltype(d1.size()) index_type;
//     index_type r1 = d1.size();
//     index_type r2 = d2.size();
//     index_type r3 = d3.size();
//     
//     typedef decltype(static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(0), A2.col(0), A3.col(0), params...))
//       scalar_type;
//     
//     // Initialize the global result
//     scalar_type global_res = 0;
//     
// #pragma omp parallel
// {
//   // Create a local copy of res for each thread
//   scalar_type local_res = 0;
//   
// #pragma omp for nowait // Distribute loop iterations without waiting at the end
//   for(index_type m1 = 0; m1 < r1; ++m1) {
//     for(index_type m2 = 0; m2 < r2; ++m2) {
//       for(index_type m3 = 0; m3 < r3; ++m3) {
//         auto square = [](scalar_type x) { return x*x; };
//         local_res += d1(m1)*d2(m2)*d3(m3)*square(static_cast<const BasicCGF*>(this)->K3operator(tvec, A1.col(m1), A2.col(m2), A3.col(m3), params...));
//       }
//     }
//   }
//   
//   // Safely add the result from this thread to the global result
// #pragma omp critical
//   global_res += local_res;
// }
// 
// return global_res;
//   }
  

}; // class CGF_Defaults

} // namespace CGFs_via_templates
} // namespace saddlepoint

#endif // CGF_DEFAULTS_H_INCLUDED
