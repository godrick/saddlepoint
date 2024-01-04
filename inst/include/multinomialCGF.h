#ifndef MULTINOMIALCGF_H_INCLUDED
#define MULTINOMIALCGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
#include "saddlepoint_types.h"
#include "CGF_Defaults.h"
#include "parametric_submodelCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

class MultinomialCGF : public CGF_Defaults<MultinomialCGF> {
public:
  // No data members or constructors needed.
  // All methods can be static methods since there are no data members,
  // and can be accessed by Multinomial_CGF::K(...) etc.
  // Alternatively, objects of type Multinomial_CGF can be created via
  // Multinomial_CGF obj;
  // and methods can be accessed via
  // obj.K(...), obj.K1(...), etc.

  // The following template methods are oriented towards types from the Eigen library.
  // For instance, it is assumed that the vector types have methods such as .array()
  // that allow component-wise multiplication and other operations.
  // Note that Eigen does not always interact well with the auto keyword
  // since it is based on creating expression templates that do not
  // directly correspond to a stored vector or matrix.
  // This is partially overcome by using the .eval() method, which converts many
  // Eigen expression types into concrete objects.
  //-----------------------------------------------------------------------
  // Helper functions to evaluate K, K1, K2, etc in terms of subexpressions
  //-----------------------------------------------------------------------
  template <class t_vector_type, class N_type, class odds_vector_type, class scalar_type>
  static scalar_type K_z1p(const t_vector_type& zm1, const N_type& N,
                           const odds_vector_type& odds_vec,
                           const scalar_type& odds_sum) {
    // Computes K_multinomial as a function of:
    // zm1, the vector whose entries are zm1[i] = exp(tvec[i])-1
    // odds_vec, the vector of odds for the different Multinomial categories
    // This vector must have non-negative entries, not all 0.
    // The Multinomial probabilities are p_i = odds_vec[i] / odds_vec.sum()
    // N, the Multinomial parameter
    // odds_sum, which must equal odds_vec.sum()
    // If odds_vec is known to be a vector of probabilities, odds_sum should equal 1 modulo numerical error
    // and the caller may supply odds_sum = 1.0 without calculating odds_vec.sum()
    // Note:
    // K = N*log(sum_i p_i e^{t_i})
    //   = N*log(sum_i (p_i + p_i (e^{t_i}-1)))
    //   = N*log(1 + sum_i odds_vec[i]*zm1[i]/odds_sum)
    return N*log1p((odds_vec.array() * zm1.array()).sum() / odds_sum);
  }
  template <class t_vector_type, class N_type, class odds_vector_type>
  static auto K_z1p(const t_vector_type& zm1, const N_type& N, const odds_vector_type& odds_vec) {
    // Convenience form that computes odds_sum directly
    return K_z1p(zm1, N, odds_vec, odds_vec.sum());
  }

  template <class z_vector_type, class odds_vector_type>
  static auto q_from_z(const z_vector_type& zvec, const odds_vector_type& odds_vec) {
    // Utility function to compute q from z, where
    // q_i = p_i e^{t_i} = p_i e^{z_i}
    // and p_i = odds_vec[i], t_i = tvec[i], z_i = zvec[i] = e^{t_i}
    // so that v_i = q_i / sum_j(q_j)
    return (zvec.array()*odds_vec.array()).matrix().eval();
  }
  template <class z_vector_type, class odds_vector_type>
  static auto v_from_z(const z_vector_type& zvec, const odds_vector_type& odds_vec) {
    // Utility function to compute v from z, where
    // v_i = p_i e^{t_i} / sum_j(p_j e^{t_j}) = p_i z_i / sum_j(p_j z_j)
    // and p_i = odds_vec[i], t_i = tvec[i], z_i = zvec[i] = e^{t_i}
    auto res = q_from_z(zvec, odds_vec);
    res *= 1/res.sum();
    return res;
  }
  template <class t_vector_type>
  static auto z_from_t(const t_vector_type& tvec) {
    // Utility function to compute v from z, where
    // v_i = p_i e^{t_i} / sum_j(p_j e^{t_j}) = p_i z_i / sum_j(p_j z_j)
    // and p_i = odds_vec[i], t_i = tvec[i], z_i = zvec[i] = e^{t_i}
    return exp(tvec.array()).matrix();
  }
  template <class t_vector_type>
  static auto zm1_from_t(const t_vector_type& tvec) {
    // Utility function to compute v from z, where
    // v_i = p_i e^{t_i} / sum_j(p_j e^{t_j}) = p_i z_i / sum_j(p_j z_j)
    // and p_i = odds_vec[i], t_i = tvec[i], z_i = zvec[i] = e^{t_i}
    // This function should use expm1() but this does not seem to be fully supported
    return (tvec.array().exp() - 1).matrix().eval();
  }
  template <class t_vector_type, class odds_vector_type>
  static auto v_from_t(const t_vector_type& tvec, const odds_vector_type& odds_vec) {
    // Utility function to compute v from z, where
    // v_i = p_i e^{t_i} / sum_j(p_j e^{t_j}) = p_i z_i / sum_j(p_j z_j)
    // and p_i = odds_vec[i], t_i = tvec[i], z_i = zvec[i] = e^{t_i}
    return v_from_z(z_from_t(tvec), odds_vec);
  }

  template <class N_type, class v_vector_type>
  static auto K1_v(const N_type& N, const v_vector_type& v) {
    // v must be the normalised vector of probabilities for the tilted distribution
    // v_i = odds_i * exp(t_i) / sum_j(odds_i * exp(t_i))
    // K1 has a simple formula K1 = N*v but K2 and related quantities have more elaborate dependence on v
    return (N*v).eval();
  }
  template <class N_type, class v_vector_type>
  static auto K2_v(const N_type& N, const v_vector_type& v) {
    // v must be the normalised vector of probabilities for the tilted distribution
    // v_i = odds_i * exp(t_i) / sum_j(odds_i * exp(t_i))
    // K2 = N * ( V - v * v^T) where
    // v is a column vector, so that (v v^T) is a square matrix; and
    // V is the diagonal matrix with entries v
    auto Nv = (N*v).eval();
    auto res = (- Nv * v.transpose()).eval();
    res.diagonal() += Nv;
    return res;
  }
  //-----------------------------------------------------------------------
  template <class t_vector_type, class N_type, class odds_vector_type>
  static auto K(const t_vector_type& tvec, const N_type& N, const odds_vector_type& odds_vec) {
    //return N * log((exp(tvec.array()) * odds_vec.array()).sum() / prob_vec.sum());
    return K_z1p(zm1_from_t(tvec), N, odds_vec);
  }
  // Note: this code does not require prob_vec to be a vector of probabilities.
  // It can be any vector of non-zero numbers not all zero.
  // In this case it can be interpreted as a vector of odds, and if we call it odds_vec
  // then the Multinomial parameters are p_i = odds_vec[i] / odds_vec.sum() for all i.
  // TO DO: make this change and document appropriately
  //-----------------------------------------------------------------------
  template <class t_vector_type, class N_type, class p_vector_type>
  static auto K1(const t_vector_type& tvec, const N_type& N, const p_vector_type& prob_vec) {
    auto v = ((exp(tvec.array()) * prob_vec.array()).matrix()).eval();
    return (v * (N / v.sum())).eval();
  }
  //-----------------------------------------------------------------------
  template <class t_vector_type, class N_type, class p_vector_type>
  static auto K2(t_vector_type tvec, N_type N, p_vector_type prob_vec) {
    auto v = (exp(tvec.array()) * prob_vec.array()).matrix().eval();
    v *= 1 / v.sum();
    auto res = ((-N)*v*v.transpose()).eval();
    res.diagonal() += N*v;
    std::cout << "Multinomial K2 dim: " << res.rows() << " " << res.cols() << std::endl;
    return res;
  }
  //-----------------------------------------------------------------------
  using CGF_Defaults::tilting_exponent;
  using CGF_Defaults::neg_ll;
  using CGF_Defaults::func_T;
  //-----------------------------------------------------------------------
  template <class t_vector_type, class w1_type, class w2_type, class w3_type, class N_type, class odds_vector_type>
  static auto K3operator(const t_vector_type& tvec, const w1_type& w1, const w2_type& w2, const w3_type& w3,
                         const N_type& N, const odds_vector_type& odds_vec) {
    return K3operator_v(w1, w2, w3, N, v_from_t(tvec, odds_vec));
  }
  template <class w1_type, class w2_type, class w3_type, class N_type, class v_vector_type>
  static auto K3operator_v(const w1_type& w1, const w2_type& w2, const w3_type& w3,
                           const N_type& N, const v_vector_type& v) {
    auto vw1 = ( v.array()*w1.array() ).eval();
    auto w2w3 = ( w2.array()*w3.array() ).eval();
    auto vw1s = vw1.sum();
    auto vw2s = ( v.array()*w2.array() ).sum();
    auto vw3s = ( v.array()*w3.array() ).sum();
    return (vw1*w2w3).sum()
      - vw3s*( vw1*w3.array() ).sum() - vw2s*( vw1*w2.array() ).sum() - vw1s*( v.array()*w2w3 ).sum()
      + 2*vw1s*vw2s*vw3s;
  }
  //-----------------------------------------------------------------------
  template <class t_vector_type, class w1_type, class w2_type, class w3_type, class w4_type, class N_type, class odds_vector_type>
  static auto K4operator(const t_vector_type& tvec, const w1_type& w1, const w2_type& w2, const w3_type& w3, const w4_type& w4,
                         const N_type& N, const odds_vector_type& odds_vec) {
    return K4operator_v(w1, w2, w3, w4, N, v_from_t(tvec, odds_vec));
  }
  template <class w1_type, class w2_type, class w3_type, class w4_type, class N_type, class v_vector_type>
  static auto K4operator_v(const w1_type& w1, const w2_type& w2, const w3_type& w3, const w4_type& w4,
                           const N_type& N, const v_vector_type& v) {
    auto vw1 = ( v.array()*w1.array() ).eval();
    auto vw2 = ( v.array()*w2.array() ).eval();
    auto vw3 = ( v.array()*w3.array() ).eval();
    auto vw4 = ( v.array()*w4.array() ).eval();
    auto vw1s = vw1.sum();
    auto vw2s = vw2.sum();
    auto vw3s = vw3.sum();
    auto vw4s = vw4.sum();
    auto w12 = ( w1.array()*w2.array() ).eval();
    auto w34 = ( w3.array()*w4.array() ).eval();
    auto vw12s = ( vw1*w2.array() ).sum();
    auto vw13s = ( vw1*w3.array() ).sum();
    auto vw14s = ( vw1*w4.array() ).sum();
    auto vw23s = ( vw2*w3.array() ).sum();
    auto vw24s = ( vw2*w4.array() ).sum();
    auto vw34s = ( vw3*w4.array() ).sum();
    auto vw123 = ( w12*vw3 ).eval();

    return (vw123*w4.array()).sum() - vw123.sum()*vw4s - (w12*vw4).sum()*vw3s - (w34*vw1).sum()*vw2s - (w34*vw2).sum()*vw1s
    - vw12s*vw34s - vw13s*vw24s - vw14s*vw23s
    + 2*(vw12s*vw3s*vw4s + vw13s*vw2s*vw4s + vw14s*vw2s*vw3s + vw23s*vw1s*vw4s + vw24s*vw1s*vw3s + vw34s*vw1s*vw3s)
    - 6*vw1s*vw2s*vw3s*vw4s;
  }
  //-----------------------------------------------------------------------
  template <class Q_matrix_type, class N_type, class v_vector_type>
  static auto K4operatorAABB_v(const Q_matrix_type& Q1, const Q_matrix_type& Q2, const N_type& N, const v_vector_type& v) {
    typedef decltype( v.sum() ) scalar_type;
    typedef v_vector_type vector_type;

    vector_type Q1v = Q1 * v;
    vector_type Q2v = Q2 * v;
    scalar_type vQ1v = v.transpose() * Q1v;
    scalar_type vQ2v = v.transpose() * Q2v;

    int length_Q1 = Q1.rows();
    scalar_type res_double_indices = 0.0;
    for(int j1 = 0; j1 < length_Q1; ++j1){
      for(int j2 = 0; j2 < length_Q1; ++j2){
        res_double_indices += v(j1)*v(j2)*Q1(j1,j2)*Q2(j1,j2);
      }
    }
    scalar_type tmp = v.transpose() * Q2.diagonal();
    return N * (-2*res_double_indices +
                (v.array()*Q1.diagonal().array() * (Q2.diagonal().array() - 2*Q2v.array() - tmp + 2*vQ2v)).sum() -
                (2*v.array() * Q1v.array() * Q2.diagonal().array()).sum() +
                (8*v.array() * Q1v.array() * Q2v.array()).sum() +
                2*vQ1v*tmp -
                6*vQ1v*vQ2v);
  }
  template <class Q_matrix_type, class t_vector_type, class N_type, class odds_vector_type>
  static auto K4operatorAABB(const t_vector_type& tvec, const Q_matrix_type& Q1, const Q_matrix_type& Q2,
                             const N_type& N, const odds_vector_type& odds_vec) {
    return K4operatorAABB_v(Q1, Q2, N, v_from_t(tvec, odds_vec));
  }

  //----------------------------------------------------------------------------
  template <class Q_matrix_type, class N_type, class v_vector_type>
  static auto K3K3operatorAABBCC_v(const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                                   const N_type& N, const v_vector_type& v) {
    typedef decltype( v.sum() ) scalar_type;
    typedef v_vector_type vector_type;

    vector_type Q1v = Q1 * v;
    vector_type Q2v = Q2 * v;
    vector_type Q3v = Q3 * v;
    scalar_type vQ1v = v.transpose() * Q1v;
    scalar_type vQ2v = v.transpose() * Q2v;
    scalar_type vQ3v = v.transpose() * Q3v;

    int length_Q1 = Q1.rows();
    scalar_type res_double_indices = 0.0;
    for(int j1 = 0; j1 < length_Q1; ++j1){
      for(int j2 = 0; j2 < length_Q1; ++j2){
        res_double_indices += ((v(j1)*v(j2)*Q2(j1,j2)) * (Q1(j1,j1)*Q3(j2,j2) -
          2*Q3(j2,j2)*Q1v(j1) -
          2*Q1(j1,j1)*Q3v(j2) +
          4*Q1v(j1)*Q3v(j2)) );
      }
    }

    return N*N*(res_double_indices +
                (v.array()*Q3.diagonal().array()*Q2v.array()).sum() * (-(v.array()*Q1.diagonal().array()).sum() + 2*vQ1v) +
                (v.array()*Q2v.array()*Q3v.array()).sum() * (2*(v.array()*Q1.diagonal().array()).sum() - 4*vQ1v) +
                (v.array()*Q3.diagonal().array()).sum() * (-(v.array()*Q1.diagonal().array()*Q2v.array()).sum() +
                vQ2v * (v.array()*Q1.diagonal().array()).sum() +
                2*(v.array()*Q1v.array()*Q2v.array()).sum() -
                2*vQ1v*vQ2v) +
                (2*vQ3v) * ((v.array()*Q1.diagonal().array()*Q2v.array()).sum() -
                vQ2v * (v.array()*Q1.diagonal().array()).sum() -
                2*(v.array()*Q1v.array()*Q2v.array()).sum() +
                2*vQ1v*vQ2v) );
  }
  template <class t_vector_type, class Q_matrix_type, class N_type, class odds_vector_type>
  static auto K3K3operatorAABBCC(const t_vector_type& tvec,
                                 const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                                 const N_type& N, const odds_vector_type& odds_vec) {
    return K3K3operatorAABBCC_v(Q1, Q2, Q3, N, v_from_t(tvec, odds_vec));
  }
  //----------------------------------------------------------------------------
  template <class Q_matrix_type, class N_type, class v_vector_type>
  static auto K3K3operatorABCABC_v(const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                                   const N_type& N, const v_vector_type& v) {
    typedef decltype( v.sum() ) scalar_type;
    typedef v_vector_type vector_type;

    scalar_type vQ1v = v.transpose() * Q1 * v;
    scalar_type vQ2v = v.transpose() * Q2 * v;
    scalar_type vQ3v = v.transpose() * Q3 * v;
    vector_type Q1v = Q1 * v;
    vector_type Q2v = Q2 * v;
    vector_type Q3v = Q3 * v;

    int length_Q1 = Q1.rows();
    scalar_type res_double_indices = 0.0;
    for(int j1 = 0; j1 < length_Q1; ++j1){
      for(int j2 = 0; j2 < length_Q1; ++j2){
        res_double_indices += (v(j1)*v(j2) ) * (Q1(j1,j2)*Q2(j1,j2) * (Q3(j1,j2) - Q3v(j1) - Q3v(j2) + vQ3v) +
          Q2(j1,j2)*Q3(j1,j2) * (-Q1v(j1) - Q1v(j2) + vQ1v) +
          Q2(j1,j2)*(Q1v(j1)*Q3v(j2) + Q1v(j2)*Q3v(j1) ) +
          Q1(j1,j2)*Q3(j1,j2)*(-Q2v(j1) - Q2v(j2) ) +
          Q1(j1,j2)*(Q2v(j1)*Q3v(j2) + Q2v(j2)*Q3v(j1) ) +
          Q3(j1,j2)*(Q1v(j1)*Q2v(j2) + Q1v(j2)*Q2v(j1) + Q1(j1,j2)*vQ2v));
      }
    }

    return N*N*(res_double_indices +
                4*(v.array()*Q1v.array()*Q2v.array()*Q3v.array()).sum() -
                4*vQ3v*(v.array()*Q1v.array()*Q2v.array()).sum() -
                4*vQ2v*(v.array()*Q1v.array()*Q3v.array()).sum() -
                4*vQ1v*(v.array()*Q2v.array()*Q3v.array()).sum() +
                4*vQ1v*vQ2v*vQ3v);
  }
  template <class t_vector_type, class Q_matrix_type, class N_type, class odds_vector_type>
  static auto K3K3operatorABCABC(const t_vector_type& tvec,
                                 const Q_matrix_type& Q1, const Q_matrix_type& Q2, const Q_matrix_type& Q3,
                                 const N_type& N, const odds_vector_type& odds_vec) {
    return K3K3operatorABCABC_v(Q1, Q2, Q3, N, v_from_t(tvec, odds_vec));
  }

  // TO DO: add _factored versions; this allows the calculation of some products with v to be computed only once

}; // class MultinomialCGF

} // namespace CGFs_via_templates


namespace CGFs_with_AD {


class MultinomialExtractor {
  // The templated Multinomial class expects its arguments in the order tvec, N, odds_vec
  // The Extractor class assumes that N and odds_vec are stored in a single vector with N in position 0 and odds_vec following
public:
  template <class vector_type>
  auto operator()(const vector_type& parameter_vector) const {
    return std::make_pair(parameter_vector[0], parameter_vector.tail(parameter_vector.size() - 1));
  }
  template <class scalar_type, class vector_type>
  vector_type pack(const scalar_type& N, const vector_type& odds_vec) const {
    vector_type parameter_vector(1 + odds_vec.size());
    parameter_vector << N, odds_vec;
    return parameter_vector;
  }
};

//using Multinomial_with_AD = CGF_with_AD_From_Extractor_And_Template<Multinomial_CGF_template_version, MultinomialExtractor>;
using MultinomialCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::CGFwithExtractor<saddlepoint::CGFs_via_templates::MultinomialCGF, MultinomialExtractor>>;
// This CGF expects its arguments as a single vector, in the order specified by MultinomialExtractor

class MultinomialModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::MultinomialCGF, const ScalarAdaptor*, const VectorAdaptor*> {
private:
  using TemplateBaseCGF = saddlepoint::CGFs_via_templates::MultinomialCGF;
  using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::MultinomialCGF, const ScalarAdaptor*, const VectorAdaptor*>;
public:
  MultinomialModelCGF(const ScalarAdaptor* n_adaptor, const VectorAdaptor* prob_vec_adaptor)
    : Base(TemplateBaseCGF(), n_adaptor, prob_vec_adaptor) {}

  MultinomialModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};
// This CGF is created by specifying two adaptors, one returning a scalar (for N) and one returning a vector (for prob_vec)
// Once created, it takes a single parameter vector of whatever form is appropriate for the specified adaptors

} // namespace CGFs_with_AD

namespace CGFs_via_virtual_functions {

// Scratchpad classes for Multinomial CGF
template <class scalar_type, class vector_type, class matrix_type>
class Multinomial_Scratchpad_Q;

template <class scalar_type, class vector_type, class matrix_type>
class Multinomial_Scratchpad : public CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad {
  // This class is compatible with CGF_with_AD_From_Extractor_And_Template<Multinomial_CGF_template_version, MultinomialExtractor>
private:
  typedef CGF_base<scalar_type, vector_type, matrix_type> CGF_type;
  typedef typename CGF_type::Scratchpad Scratchpad_type;
  typedef typename CGF_type::Scratchpad_Q Scratchpad_Q_type;
  typedef CGFs_via_templates::MultinomialCGF MultinomialCGF;
  typedef CGFs_with_AD::MultinomialExtractor MultinomialExtractor;

  vector_type v;
  scalar_type N;
  // Most Multinomial CGF quantities depend only on v and N.
  // However, the values of K() and tilting_exponent do not.
  // They are computed and stored separately.
  scalar_type k;
  scalar_type te;

  struct K1_finder {
    vector_type operator()(const scalar_type& N, const vector_type& v) {return MultinomialCGF::K1_v(N, v);}
  };
  SavedResult<vector_type, K1_finder> saved_k1;

  struct K2_finder {
    matrix_type operator()(const scalar_type& N, const vector_type& v) {return MultinomialCGF::K2_v(N, v);}
  };
  SavedResult<matrix_type, K2_finder> saved_k2;

  struct neg_ll_finder {
    scalar_type operator()(const scalar_type& te, const matrix_type& k2) {
      // Note that k2.determinant() == 0 modulo numerical error, so this function gives problematic values if called
      return 0.5 * (log(k2.determinant()) + k2.rows() * log(2*M_PI)) - te;
    }
  };
  SavedResult<scalar_type, neg_ll_finder> saved_neg_ll;

  struct func_T_finder {
    scalar_type operator()(const scalar_type& N, const vector_type& v, const matrix_type& k2, Multinomial_Scratchpad* sp) {
      // Note that k2.determinant() == 0 modulo numerical error, so this function gives problematic values if called
      matrix_type mat_Q = k2.inverse();
      std::unique_ptr<Scratchpad_Q_type> p_spq(sp->Multinomial_Scratchpad::scratchpad_q(mat_Q));
      // Extend sp to its Scratchpad_Q equivalent, with Q == mat_Q == K2.inverse()

      return p_spq->K4operatorAABB()/8 - p_spq->K3K3operatorAABBCC()/8 - p_spq->K3K3operatorABCABC()/12;
    }
  };
  SavedResult<scalar_type, func_T_finder> saved_func_T;

  friend class Multinomial_Scratchpad_Q<scalar_type, vector_type, matrix_type>;

public:
  Multinomial_Scratchpad(const vector_type& tvec, const vector_type& parameter_vector) {
    auto p = MultinomialExtractor()(parameter_vector);
    N = p.first;
    const auto& odds_vec = p.second; // alias
    v = MultinomialCGF::v_from_t(tvec, odds_vec);
    k = MultinomialCGF::K(tvec, N, odds_vec);
    te = k - N*(tvec.array() * odds_vec.array()).sum();
  }

  Scratchpad_Q_type* scratchpad_q(const matrix_type& Q) override {
    return new Multinomial_Scratchpad_Q<scalar_type, vector_type, matrix_type>(Q, this);
  }

  const scalar_type& K() override {return k;}

  const vector_type& K1() override {return saved_k1(N, v);}
  const matrix_type& K2() override {return saved_k2(N, v);}

  const scalar_type& tilting_exponent() override {return te;}

  // Note: neg_ll and func_T both return infinite or NaN values, in principle,
  // because the covariance matrix for this CGF is singular
  const scalar_type& neg_ll() override {return saved_neg_ll(te, Multinomial_Scratchpad::K2());}
  const scalar_type& func_T() override {return saved_func_T(N, v, Multinomial_Scratchpad::K2(), this);}
};

template <>
Multinomial_Scratchpad<a_scalar, a_vector, a_matrix>::Multinomial_Scratchpad(const a_vector& tvec, const a_vector& parameter_vector) {
  auto p = MultinomialExtractor()(parameter_vector);
  N = p.first;
  const auto& odds_vec = p.second; // alias
  v = MultinomialCGF::v_from_t(tvec, odds_vec);
  k = MultinomialCGF::K(tvec, N, odds_vec);
  te = k - N*(tvec.array() * v.array()).sum();
}


template <class scalar_type, class vector_type, class matrix_type>
class Multinomial_Scratchpad_Q : public CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q {
private:
  typedef CGF_base<scalar_type, vector_type, matrix_type> CGF_type;
  typedef typename CGF_type::Scratchpad Scratchpad_type;
  typedef typename CGF_type::Scratchpad_Q Scratchpad_Q_type;
  typedef Multinomial_Scratchpad<scalar_type, vector_type, matrix_type> MSP_type;
  typedef CGFs_via_templates::MultinomialCGF MultinomialCGF;

  MSP_type* plain_sp;
  // non-owned pointer to Scratchpad that created this
  // Note: data member, not base class, to avoid copying already calculated values
  // However this means that the virtual functions already implemented by plain_sp have to be explicitly re-implemented to re-route them to plain_sp

  matrix_type Q;

  struct K4operatorAABB_finder {
    scalar_type operator()(const scalar_type& N, const vector_type& v, const matrix_type& q) {
      return MultinomialCGF::K4operatorAABB_v(q, q, N, v);
    }
  };
  SavedResult<scalar_type, K4operatorAABB_finder> saved_k4aabb;

  struct K3K3operatorAABBCC_finder {
    scalar_type operator()(const scalar_type& N, const vector_type& v, const matrix_type& q) {
      return MultinomialCGF::K3K3operatorAABBCC_v(q, q, q, N, v);
    }
  };
  SavedResult<scalar_type, K3K3operatorAABBCC_finder> saved_k3k3aabbcc;

  struct K3K3operatorABCABC_finder {
    scalar_type operator()(const scalar_type& N, const vector_type& v, const matrix_type& q) {
      return MultinomialCGF::K3K3operatorABCABC_v(q, q, q, N, v);
    }
  };
  SavedResult<scalar_type, K3K3operatorABCABC_finder> saved_k3k3abcabc;

public:
  Multinomial_Scratchpad_Q(const matrix_type& q, MSP_type* psp)
    : plain_sp(psp), Q(q) {}
  // Normally this constructor will only be used by Multinomial_Scratchpad::scratchpad_q

  Scratchpad_Q_type* scratchpad_q(const matrix_type& q) override {
    // Note: method creates a new Scratchpad_Q object with a new value q, in addition to the current one
    return plain_sp->scratchpad_q(q);
    // Already implemented by plain_sp
  }

  const scalar_type& K() override {return plain_sp->K();}
  const vector_type& K1() override {return plain_sp->K1();}
  const matrix_type& K2() override {return plain_sp->K2();}
  const scalar_type& tilting_exponent() override {return plain_sp->tilting_exponent();}
  const scalar_type& neg_ll() override {return plain_sp->neg_ll();}
  const scalar_type& func_T() override {return plain_sp->func_T();}
  // Already implemented by plain_sp

  const scalar_type& K4operatorAABB() override {
    return saved_k4aabb(plain_sp->N, plain_sp->v, Q);
  }
  const scalar_type& K3K3operatorAABBCC() override {
    return saved_k3k3aabbcc(plain_sp->N, plain_sp->v, Q);
  }
  const scalar_type& K3K3operatorABCABC() override {
    return saved_k3k3abcabc(plain_sp->N, plain_sp->v, Q);
  }
};

//template <>
//CGF_base<double, vec, mat>::Scratchpad* Multinomial_with_AD::scratchpad(const vec& tvec, const vec& parameter_vector) const {
//    return new Multinomial_Scratchpad<double, vec, mat>(tvec, parameter_vector);
//}
//template <>
//CGF_base<a_scalar, a_vector, a_matrix>::Scratchpad* Multinomial_with_AD::scratchpad(const a_vector& tvec, const a_vector& parameter_vector) const {
//    return new Multinomial_Scratchpad<a_scalar, a_vector, a_matrix>(tvec, parameter_vector);
//}


} // namespace CGFs_via_virtual_functions
} // namespace saddlepoint


#endif // MULTINOMIALCGF_H_INCLUDED
