#ifndef NEGBINCGF_H_INCLUDED
#define NEGBINCGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

// ----------------------------------------------------------------------------
// Negative Binomial vectorised CGF (i.i.d. version)
// ----------------------------------------------------------------------------
class NegBinVectorisedIID_template_version
  : public ScalarVectorisedIID_Defaults<NegBinVectorisedIID_template_version>
{
public:
  
  // ------------------------------------------------------------------------
  // K(tvec, r, p):
  // K(t) = r * log( p / [1 - (1-p)*exp(t)] )
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto K_vectorised_iid(const t_vector_type& tvec,
                               const r_type& r,
                               const p_type& p)
  {
    // alpha(t) = 1 - (1-p)*exp(t)
    // => K(t) = r * [ log(p) - log(alpha(t)) ]
    auto alpha = (1 - (1 - p) * tvec.array().exp());
    return r*(log(p) - log(alpha.array()) );
  }
  
  // ------------------------------------------------------------------------
  // K'(tvec, r, p):
  //   = r * (1-p)*exp(t) / [alpha(t)]
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto K1_vectorised_iid(const t_vector_type& tvec,
                                const r_type& r,
                                const p_type& p)
  {
    auto num_   = (1-p)*tvec.array().exp();
    auto denom_ = 1 - (1 - p)*tvec.array().exp();
    return r * ( num_.array() / denom_.array() );
  }
  
  // ------------------------------------------------------------------------
  // K''(tvec, r, p):
  //   = r*(1-p)*exp(t) / [alpha(t)]^2
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto K2_vectorised_iid(const t_vector_type& tvec,
                                const r_type& r,
                                const p_type& p)
  {
    auto num_   = (1 - p)*tvec.array().exp();
    auto denom_ = 1 - (1 - p)*tvec.array().exp();
    return r*num_.array() / (denom_.array()*denom_.array());
  }
  
  // ------------------------------------------------------------------------
  // K'''(tvec, r, p):
  //   = r*(1-p)*exp(t) * [1 + (1-p)*exp(t)] / [alpha(t)]^3
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto K3_vectorised_iid(const t_vector_type& tvec,
                                const r_type& r,
                                const p_type& p)
  {
    auto e_t      = tvec.array().exp();
    auto alpha_    = 1 - (1 - p)*e_t;
    auto factor_   = (1 - p)*e_t * (1 + (1 - p)*e_t);
    return r*factor_.array() / (alpha_.array()*alpha_.array()*alpha_.array());
  }
  
  // ------------------------------------------------------------------------
  // K''''(tvec, r, p):
  //   = r*(1-p)*exp(t) * [1 + exp(2t) + 4exp(t)
  //                       - 2p exp(2t) - 4p exp(t) + p^2 exp(2t)] / [alpha]^4
  //   [Matches geometric code but multiplied by r]
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto K4_vectorised_iid(const t_vector_type& tvec,
                                const r_type& r,
                                const p_type& p)
  {
    auto e_t    = tvec.array().exp();
    auto e_2t   = (tvec.array()*2.0).exp();
    auto alpha_  = 1.0 - (1.0 - p)*e_t;
    
    // The big bracket from the geometric distribution's 4th derivative:
    // [1 + e^{2t} + 4 e^t - 2 p e^{2t} - 4 p e^t + p^2 e^{2t}]
    auto bracket_ = 1.0
    + e_2t.array()
    + 4.0*e_t.array()
    - 2.0*p*e_2t.array()
    - 4.0*p*e_t.array()
    + p*p*e_2t.array();
    
    auto numerator_ = (1.0 - p)*e_t * bracket_;
    return r*numerator_/(alpha_.array()*alpha_.array()*alpha_.array()*alpha_.array());
  }
  
  // ------------------------------------------------------------------------
  // Inequality constraint:  (1-p)*exp(t) < 1  => t < -log(1-p)
  // Return something that must be < 0. E.g. (1-p)*exp(t) - 1 < 0
  // ------------------------------------------------------------------------
  template <class t_vector_type, class r_type, class p_type>
  static auto ineq_constraint_vectorised_iid(const t_vector_type& tvec,
                                             const r_type&       r,
                                             const p_type&       p)
  {
    // The domain is:  (1-p)*exp(t) < 1.
    return (1.0 - p) * tvec.array().exp() - 1.0;
  }
};

// ----------------------------------------------------------------------------
// A standard CGF type using VectorOfIIDCGF
// ----------------------------------------------------------------------------
using NegBinCGF = VectorOfIIDCGF<NegBinVectorisedIID_template_version>;

} // namespace CGFs_via_templates



// ----------------------------------------------------------------------------
// Some wrappers: NegBinExtractorIID & NegBinModelCGF
// ----------------------------------------------------------------------------
namespace CGFs_with_AD {

class NegBinExtractorIID
{
  // We assume the parameter vector has length=2:
  //    param[0] = r
  //    param[1] = p
public:
  template <class vector_type>
  auto operator()(const vector_type& parameter_vector) const
  {
    // Return (r, p)
    return std::make_pair(parameter_vector[0], parameter_vector[1]);
  }
  
  // Packing convenience
  vec pack(double r, double p) const {
    vec parameter_vector(2);
    parameter_vector << r, p;
    return parameter_vector;
  }
  a_vector pack(a_scalar r, a_scalar p) const {
    a_vector parameter_vector(2);
    parameter_vector << r, p;
    return parameter_vector;
  }
};

using NegBinCGF = CGF_with_AD_from_template<
  saddlepoint::CGFs_via_templates::CGFwithExtractor<
    saddlepoint::CGFs_via_templates::NegBinCGF,
    NegBinExtractorIID
  >
>;



// ----------------------------------------------------------------------------
// NegBinModelCGF
// If you want a "model" that takes adaptors for r & p
// ----------------------------------------------------------------------------
class NegBinModelCGF
  : public SeparateParametersCGF<
    saddlepoint::CGFs_via_templates::NegBinCGF,
    const ScalarAdaptor*,
    const ScalarAdaptor*
  >
{
private:
  using TemplateBaseCGF = saddlepoint::CGFs_via_templates::NegBinCGF;
  using Base = SeparateParametersCGF<TemplateBaseCGF,
                                     const ScalarAdaptor*,
                                     const ScalarAdaptor*>;
public:
  // Constructor that takes two adaptor pointers: r_adaptor, p_adaptor
  NegBinModelCGF(const ScalarAdaptor* r_adaptor,
                 const ScalarAdaptor* p_adaptor)
    : Base(TemplateBaseCGF(), r_adaptor, p_adaptor)
  {}
  
  NegBinModelCGF() = delete; 
};

} // namespace CGFs_with_AD

} // namespace saddlepoint

#endif // NEGBINCGF_H_INCLUDED
