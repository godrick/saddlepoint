#include "baseCGF.h"

// using namespace saddlepoint::CGFs_with_AD;
// using namespace saddlepoint::CGFs_via_templates;



namespace saddlepoint {
namespace CGFs_with_AD {
/**
 * 
 * 
 * This class `CGF_with_AD_from_r_functions` extends `CGF_with_AD` and allows for the utilization of
 * vectorized scalar CGFs computations directly from R functions. 
 * It's meant for users who want to incorporate their own CGF functions into the saddlepoint package framework.
 * The advantage with this is that it allows for the use of R functions directly as opposed to C++ functions.
 * 
 * 
 * Current limitations and future work:
 * - This class assumes that all CGF-related functions are scalar and vectorized (iid vectorized). 
 *   CGF-related functions for multivariate rv.s are not supported; that can still be done using C++ functions.
 * - For instance, this class currently expects that the R function for K2 returns a vector which will be used to form the diagonal of a matrix.
 *   This implies that off-diagonal elements are zero, which is efficient for diagonal matrices but restrictive
 *   for full matrices (in multivariate rv.s). The class can be extended to handle full matrices directly from R.
 * 
 * Important:
 * All R functions must be compatible with codes that interact with the RTMB package to ensure
 * AD functionality.
 */
class CGF_with_AD_from_r_functions : public CGF_with_AD {
private:
  Rcpp::Function Kvectorized_rfunction;
  Rcpp::Function K1vectorized_rfunction;
  Rcpp::Function K2vectorized_rfunction;
  Rcpp::Function K3vectorized_rfunction;  
  Rcpp::Function K4vectorized_rfunction;  
  Rcpp::Function ineq_constraint_vectorized_rfunction; 
  
public:
  CGF_with_AD_from_r_functions(Rcpp::Function K_rfunc, Rcpp::Function K1_rfunc, Rcpp::Function K2_rfunc, 
                               Rcpp::Function K3_rfunc, Rcpp::Function K4_rfunc, Rcpp::Function ineq_constraint_rfunc)
    : Kvectorized_rfunction(K_rfunc), K1vectorized_rfunction(K1_rfunc), K2vectorized_rfunction(K2_rfunc),
      K3vectorized_rfunction(K3_rfunc), K4vectorized_rfunction(K4_rfunc), ineq_constraint_vectorized_rfunction(ineq_constraint_rfunc) {}
  
  
  
  double K(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector result = Kvectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    return Rcpp::sum(result); 
  }
  a_scalar K(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector f_advector = Kvectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    // a_vector a_res = get_a_vector_from_advector(f_advector);
    return get_a_vector_from_advector(f_advector).sum();
  }
  // ---------------------------
  // ---------------------------
  vec K1(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector result = K1vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec res(result.size());
    for(int i = 0; i < result.size(); i++) { res[i] = result[i]; }
    return res;
  }
  a_vector K1(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector f_advector = K1vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    return get_a_vector_from_advector(f_advector);
  }
  // ---------------------------
  // ---------------------------
  
  // 
  mat K2(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector result = K2vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    // // int index_ = 0;
    // mat res(tvec.size(), tvec.size());
    // for(int j = 0; j < res.cols(); j++) { 
    //   for (int i = 0; i < res.rows(); i++){
    //     // res(i,j) = result[index_++]; 
    //     res(i,j) = result(i,j);
    //   }
    // }
    vec res(result.size());
    for(int i = 0; i < result.size(); i++) { res[i] = result[i]; }
    return res.matrix().asDiagonal();
  }
  a_matrix K2(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector f_advector = K2vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    return get_a_vector_from_advector(f_advector).matrix().asDiagonal();
  }
  
  double tilting_exponent(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector res_K = Kvectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    Rcpp::NumericVector res_K1 = K1vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K[i] = res_K[i]; }
    vec vecres_K1(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K1[i] = res_K1[i]; }
    return (vecres_K.array() - tvec.array()*vecres_K1.array()).sum();
  }
  a_scalar tilting_exponent(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector rfK_advector = Kvectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    RADvector rfK1_advector = K1vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K_vals = get_a_vector_from_advector(rfK_advector);
    a_vector K1_vals = get_a_vector_from_advector(rfK1_advector);
    return (K_vals.array() - tvec.array()*K1_vals.array()).sum();
  }
  
  
  double neg_ll(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector K_vals = Kvectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    Rcpp::NumericVector K1_vals = K1vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    Rcpp::NumericVector K2_vals = K2vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K[i] = K_vals[i]; }
    vec vecres_K1(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K1[i] = K1_vals[i]; }
    vec vecres_K2(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K1[i] = K2_vals[i]; }
  
    return (0.5*(2*M_PI*vecres_K2.array()).log() - (vecres_K.array() - tvec.array()*vecres_K1.array())).sum();
  }
  a_scalar neg_ll(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector rfK_advector = Kvectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    RADvector rfK1_advector = K1vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    RADvector rfK2_advector = K2vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K_vals = get_a_vector_from_advector(rfK_advector);
    a_vector K1_vals = get_a_vector_from_advector(rfK1_advector);
    a_vector K2_vals = get_a_vector_from_advector(rfK2_advector);
    return (0.5*(2*M_PI*K2_vals.array()).log() - (K_vals.array() - tvec.array()*K1_vals.array())).sum();
  }
  
  
  double func_T(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector K2_vals = K2vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    Rcpp::NumericVector K3_vals = K3vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    Rcpp::NumericVector K4_vals = K4vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K2(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K2[i] = K2_vals[i]; }
    vec vecres_K3(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K3[i] = K3_vals[i]; }
    vec vecres_K4(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K4[i] = K4_vals[i]; }
    vec k2sq_val = vecres_K2.array() * vecres_K2.array();
    return (vecres_K4.array()/(8*k2sq_val.array()) - 5*vecres_K3.array()*vecres_K3.array()/(24*k2sq_val.array()*vecres_K2.array())).sum();
  }
  a_scalar func_T(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector rfK2_advector = K2vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    RADvector rfK3_advector = K3vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    RADvector rfK4_advector = K4vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K2_vals = get_a_vector_from_advector(rfK2_advector);
    a_vector K3_vals = get_a_vector_from_advector(rfK3_advector);
    a_vector K4_vals = get_a_vector_from_advector(rfK4_advector);
    a_vector k2sq_val = K2_vals.array() * K2_vals.array();
    return (K4_vals.array()/(8*k2sq_val.array()) - 5*K3_vals.array()*K3_vals.array()/(24*k2sq_val.array()*K2_vals.array())).sum();
  }
  
  double K2operator(const vec& tvec, const vec& x, const vec& y, const vec& parameter_vector) const override {
    return (x.transpose() * K2(tvec, parameter_vector) * y).value();
  }
  a_scalar K2operator(const a_vector& tvec, const a_vector& x, const a_vector& y, const a_vector& parameter_vector) const override {
    return (x.transpose() * K2(tvec, parameter_vector) * y).value();
  }
  
  mat K2operatorAK2AT(const vec& tvec, const mat& A, const vec& parameter_vector) const override {
    return A * K2(tvec, parameter_vector) * A.transpose();
  }
  a_matrix K2operatorAK2AT(const a_vector& tvec, const a_matrix& A, const a_vector& parameter_vector) const override {
    return A * K2(tvec, parameter_vector) * A.transpose();
  }
  
  double K3operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& parameter_vector) const override {
    Rcpp::NumericVector K3_vals = K3vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K3(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K3[i] = K3_vals[i]; }
    return ( vecres_K3.array() * v1.array() * v2.array() * v3.array() ).sum();
  }
  a_scalar K3operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& parameter_vector) const override {
    RADvector rfK3_advector = K3vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K3_vals = get_a_vector_from_advector(rfK3_advector);
    return ( K3_vals.array() * v1.array() * v2.array() * v3.array() ).sum();
  }
  
  double K4operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& v4, const vec& parameter_vector) const override {
    Rcpp::NumericVector K4_vals = K4vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K4(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K4[i] = K4_vals[i]; }
    return ( vecres_K4.array() * v1.array() * v2.array() * v3.array() * v4.array() ).sum();
  }
  a_scalar K4operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& v4, const a_vector& parameter_vector) const override {
    RADvector rfK4_advector = K4vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K4_vals = get_a_vector_from_advector(rfK4_advector);
    return ( K4_vals.array() * v1.array() * v2.array() * v3.array() * v4.array() ).sum();
  }
  
  
  double K4operatorAABB(const vec& tvec, const mat& Q1, const mat& Q2, const vec& parameter_vector) const override {
    // K4 is a diagonal tensor so the result only uses the diagonals of Q1,Q2
    Rcpp::NumericVector K4_vals = K4vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K4(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K4[i] = K4_vals[i]; }
    return ( vecres_K4.array() * Q1.diagonal().array() * Q2.diagonal().array() ).sum();
  }
  a_scalar K4operatorAABB(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_vector& parameter_vector) const override {
    RADvector rfK4_advector = K4vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K4_vals = get_a_vector_from_advector(rfK4_advector);
    return ( K4_vals.array() * Q1.diagonal().array() * Q2.diagonal().array() ).sum();
  }
  
  
  double K3K3operatorAABBCC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {
    Rcpp::NumericVector K3_vals = K3vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K3(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K3[i] = K3_vals[i]; }
    return ( (Q1.diagonal().array()*vecres_K3.array() ).matrix().asDiagonal() * Q2 * (Q3.diagonal().array()*vecres_K3.array() ).matrix().asDiagonal() ).sum() ;
  }
  a_scalar K3K3operatorAABBCC(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {
    RADvector rfK3_advector = K3vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K3_vals = get_a_vector_from_advector(rfK3_advector);
    return ( (Q1.diagonal().array()*K3_vals.array() ).matrix().asDiagonal() * Q2 * (Q3.diagonal().array()*K3_vals.array() ).matrix().asDiagonal() ).sum() ;
  }
  
  
  double K3K3operatorABCABC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {
    Rcpp::NumericVector K3_vals = K3vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec vecres_K3(tvec.size()); for(int i = 0; i < tvec.size(); i++) { vecres_K3[i] = K3_vals[i]; }
    return ( vecres_K3.matrix().asDiagonal() * (Q1.diagonal().array()*Q2.diagonal().array()*Q3.diagonal().array()).matrix() * vecres_K3.matrix().asDiagonal() ).sum() ;
  }
  a_scalar K3K3operatorABCABC(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {
    RADvector rfK3_advector = K3vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    a_vector K3_vals = get_a_vector_from_advector(rfK3_advector);
    return (K3_vals.matrix().asDiagonal() * (Q1.diagonal().array()*Q2.diagonal().array()*Q3.diagonal().array()).matrix() * K3_vals.matrix().asDiagonal() ).sum() ;
  }
  
  double K4operatorAABB_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const vec& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    
    double res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for (index_type m2 = 0; m2 < r2; ++m2) {
        res += d1(m1) * d2(m2) * K4operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), A2.col(m2), parameter_vector);
      }
    }
    return res;
  }
  a_scalar K4operatorAABB_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_vector& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    
    a_scalar res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for (index_type m2 = 0; m2 < r2; ++m2) {
        res += d1(m1) * d2(m2) * K4operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), A2.col(m2), parameter_vector);
      }
    }
    return res;
  }
  double K3K3operatorAABBCC_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const mat& A3, const vec& d3, const vec& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();
    
    double res = 0;
    for(index_type m2 = 0; m2 < r2; ++m2) {
      double factor1 = 0;
      for(index_type m1 = 0; m1 < r1; ++m1) {
        factor1 += d1(m1) * K3operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), parameter_vector);
      }
      double factor2 = 0;
      for(index_type m3 = 0; m3 < r3; ++m3) {
        factor2 += d3(m3) * K3operator(tvec, A2.col(m2), A3.col(m3), A3.col(m3), parameter_vector);
      }
      res += d2(m2)*factor1*factor2;
    }
    return res;
  }
  a_scalar K3K3operatorAABBCC_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_matrix& A3, const a_vector& d3, const a_vector& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();
    
    a_scalar res = 0;
    for(index_type m2 = 0; m2 < r2; ++m2) {
      a_scalar factor1 = 0;
      for(index_type m1 = 0; m1 < r1; ++m1) {
        factor1 += d1(m1) * K3operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), parameter_vector);
      }
      a_scalar factor2 = 0;
      for(index_type m3 = 0; m3 < r3; ++m3) {
        factor2 += d3(m3) * K3operator(tvec, A2.col(m2), A3.col(m3), A3.col(m3), parameter_vector);
      }
      res += d2(m2)*factor1*factor2;
    }
    return res;
  }
  
  
  
  vec ineq_constraint(const vec& tvec, const vec& parameter_vector) const override {
    Rcpp::NumericVector result = ineq_constraint_vectorized_rfunction(Rcpp::wrap(tvec), Rcpp::wrap(parameter_vector));
    vec res(result.size());
    for(int i = 0; i < result.size(); i++) { res[i] = result[i]; }
    return res;
  }
  a_vector ineq_constraint(const a_vector& tvec, const a_vector& parameter_vector) const override {
    RADvector f_advector = ineq_constraint_vectorized_rfunction(send_a_vector_to_advector(tvec), send_a_vector_to_advector(parameter_vector));
    return get_a_vector_from_advector(f_advector);
  }
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  double K3K3operatorABCABC_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const mat& A3, const vec& d3, const vec& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();
    
    double res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for(index_type m2 = 0; m2 < r2; ++m2) {
        for(index_type m3 = 0; m3 < r3; ++m3) {
          auto square = [](double x) {return x*x;};
          res += d1(m1)*d2(m2)*d3(m3)*square(K3operator(tvec, A1.col(m1), A2.col(m2), A3.col(m3), parameter_vector));
        }
      }
    }
    return res;
  }
  a_scalar K3K3operatorABCABC_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_matrix& A3, const a_vector& d3, const a_vector& parameter_vector) const override {
    typedef decltype(d1.size()) index_type;
    index_type r1 = d1.size();
    index_type r2 = d2.size();
    index_type r3 = d3.size();
    
    a_scalar res = 0;
    for(index_type m1 = 0; m1 < r1; ++m1) {
      for(index_type m2 = 0; m2 < r2; ++m2) {
        for(index_type m3 = 0; m3 < r3; ++m3) {
          auto square = [](a_scalar x) {return x*x;};
          res += d1(m1)*d2(m2)*d3(m3)*square(K3operator(tvec, A1.col(m1), A2.col(m2), A3.col(m3), parameter_vector));
        }
      }
    }
    return res;
  }
};




} // namespace CGFs_with_AD
} // namespace saddlepoint













