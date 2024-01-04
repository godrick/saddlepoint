#ifndef GAMMACGF_H_INCLUDED
#define GAMMACGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
// #include "saddlepoint_types.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"


namespace saddlepoint {
namespace CGFs_via_templates {

// The code here follows the steps in GeometricCGF.h but we omit the NonIdenticalVersion of the CGF

// CGF of a Gamma distributed r.v with parameter shape and rate
// K(tvec, shape, rate) = -shape*log(1 - tvec/rate); tvec < rate

class GammaVectorisedIID_template_version : public ScalarVectorisedIID_Defaults<GammaVectorisedIID_template_version> {

public:
    // These are vectorised functions of K and its derivatives
    // The functions return a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // shape and rate arguments remain constant across all values of tvec
    template <class t_vector_type, class shape_type, class rate_type>
    static auto K_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        return -shape*log1p(-tvec.array()/rate);
    }
//----------------
    template <class t_vector_type, class shape_type, class rate_type>
    static auto K1_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        return shape / (rate - tvec.array());
    }
//----------------
    template <class t_vector_type, class shape_type, class rate_type>
    static auto K2_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        return shape / ((rate - tvec.array())*(rate - tvec.array()));
    }
//----------------
    template <class t_vector_type, class shape_type, class rate_type>
    static auto K3_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        auto temp = rate - tvec.array();
        return 2*shape/(temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class shape_type, class rate_type>
    static auto K4_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        auto temp = rate - tvec.array();
        return 6*shape/(temp.array()*temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class shape_type, class rate_type>
    static auto ineq_constraint_vectorised_iid(const t_vector_type& tvec, const shape_type& shape, const rate_type& rate) {
        // tvec < rate
        // The value of this function is constrained to be non-positive
        // tvec - rate < 0
        return tvec.array() - rate;
    }
};


// VectorOfIIDCGF class use thentemplate vectorised methods to obtain a standard CGF template version
using GammaCGF = VectorOfIIDCGF<GammaVectorisedIID_template_version>;



} // namespace CGFs_via_templates

namespace CGFs_with_AD {

class GammaExtractorIID {
// The templated GammaCGF class expects only p as its argument
// This Extractor class assumes that p is the only scalar and therefore stored at positon 0
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return std::make_pair(parameter_vector[0], parameter_vector[1]);
    }
};


using GammaCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::CGFwithExtractor<saddlepoint::CGFs_via_templates::GammaCGF, GammaExtractorIID>>;

class GammaModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::GammaCGF, const ScalarAdaptor*, const ScalarAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::GammaCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::GammaCGF, const ScalarAdaptor*, const ScalarAdaptor*>;
public:
    GammaModelCGF(const ScalarAdaptor* shape_adaptor, const ScalarAdaptor* rate_adaptor) : Base(TemplateBaseCGF(), shape_adaptor, rate_adaptor) {}

    GammaModelCGF() = delete;
};


} // namespace CGFs_with_AD

} // namespace saddlepoint


#endif // GAMMACGF_H_INCLUDED
