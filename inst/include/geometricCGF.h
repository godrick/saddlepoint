#ifndef GEOMETRICCGF_H_INCLUDED
#define GEOMETRICCGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
// #include "saddlepoint_types.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"


namespace saddlepoint {
namespace CGFs_via_templates {

// The GeometricVectorisedIID_template_version class is a template specialisation for the case when the Geometric distribution is i.i.d.
// The GeometricVectorisedNonIdentical_template_version class is a template case for non-identical distribution.
// The templates provide vectorised functions to evaluate the CGF and its derivatives for a vector of values tvec and a scalar parameter p
 // for the i.i.d. case and a vector parameter p_vector_type for the non-identical case.
// The code defines K_vectorised_iid, K1_vectorised_iid, K2_vectorised_iid, K3_vectorised_iid, K4_vectorised_iid, and ineq_constraint_vectorised_iid functions for the i.i.d. case.
// The code also defines the same functions for the non-identical case.
 // K, K1, K2, K3, and K4 evaluate the CGF and its derivatives, and ineq_constraint returns the inequality constraint tvec < -log(1-p).

//// CGF of a Geometric distribution with parameter p
// PMF: P(X=k) = (1-p)^k * p; k = 0,1,...
// K(tvec, p) = log(p) - log(1 - exp(tvec) + p * exp(tvec)); tvec < -log(1-p)

class GeometricVectorisedIID_template_version : public ScalarVectorisedIID_Defaults<GeometricVectorisedIID_template_version> {

public:

    // These are vectorised functions of K and its derivatives
    // The functions return a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // p as argument remain constant across all values of tvec
    template <class t_vector_type, class p_type>
    static auto K_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        return log(p)-log(1 - exp(tvec.array()) + p*exp(tvec.array()));
    }
//----------------
    template <class t_vector_type, class p_type>
    static auto K1_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        return (exp(tvec.array()) - p*exp(tvec.array())) / (1 - exp(tvec.array()) + p*exp(tvec.array()));
    }
//----------------
    template <class t_vector_type, class p_type>
    static auto K2_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        auto temp = 1 - exp(tvec.array()) + p*exp(tvec.array());
        return (exp(tvec.array()) - p*exp(tvec.array())) / ( temp.array() * temp.array() );
        // return (exp(tvec.array()) - p*exp(tvec.array())) / vector_pow(1 - exp(tvec.array()) + p*exp(tvec.array()), 2);
    }
//----------------
    template <class t_vector_type, class p_type>
    static auto K3_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        auto temp = 1 - exp(tvec.array()) + p*exp(tvec.array());
        return (exp(tvec.array()) - p*exp(tvec.array())) * (1 + exp(tvec.array()) - p*exp(tvec.array()) ) / (temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class p_type>
    static auto K4_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        auto temp = 1 - exp(tvec.array()) + p*exp(tvec.array());
        return (exp(tvec.array()) - p*exp(tvec.array())) * (1 + exp(2*tvec.array()) + 4*exp(tvec.array()) -
                                                       2*p*exp(2*tvec.array()) - 4*p*exp(tvec.array()) +
                                                       p*p*exp(2*tvec.array())  ) / (temp.array()*temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class p_type>
    static auto ineq_constraint_vectorised_iid(const t_vector_type& tvec, const p_type& p) {
        // tvec < -log(1-p)
        // The value of this function is constrained to be non-positive
        // tvec + log(1 - p) < 0
        // return tvec.array() + log(1 - p);
        return (1-p)*exp(tvec.array()) - 1;
    }
};
//----------------
//----------------
//----------------
class GeometricVectorisedNonIdentical_template_version : public ScalarVectorisedNonIdentical_Defaults<GeometricVectorisedNonIdentical_template_version> {

public:
    // The NonIdentical versions of the K and its derivatives
    // These methods accept p as a vector.
    // Each value in tvec has its own p.
    // tvec.size() == p.size()
    template <class t_vector_type, class p_vector_type>
    static auto K_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        return log(p.array())-log(1 - exp(tvec.array()) + p.array()*exp(tvec.array()));
    }
//----------------
    template <class t_vector_type, class p_vector_type>
    static auto K1_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        return (exp(tvec.array()) - p.array()*exp(tvec.array())) / (1 - exp(tvec.array()) + p.array()*exp(tvec.array()));
    }
//----------------
    template <class t_vector_type, class p_vector_type>
    static auto K2_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        auto temp = 1 - exp(tvec.array()) + p.array()*exp(tvec.array());
        return (exp(tvec.array()) - p.array()*exp(tvec.array())) / (temp.array()*temp.array());
        // return (exp(tvec.array()) - p.array()*exp(tvec.array())) / vector_pow(1 - exp(tvec.array()) + p.array()*exp(tvec.array()), 2);
    }
//----------------
    template <class t_vector_type, class p_vector_type>
    static auto K3_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        auto temp = 1 - exp(tvec.array()) + p.array()*exp(tvec.array());
        return (exp(tvec.array()) - p.array()*exp(tvec.array())) * (1 + exp(tvec.array()) - p.array()*exp(tvec.array()) ) / (temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class p_vector_type>
    static auto K4_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        auto temp = 1 - exp(tvec.array()) + p.array()*exp(tvec.array());
        return (exp(tvec.array()) - p.array()*exp(tvec.array())) * (1 + exp(2*tvec.array()) + 4*exp(tvec.array()) -
                                                       2*p.array()*exp(2*tvec.array()) - 4*p.array()*exp(tvec.array()) +
                                                       p.array()*p.array()*exp(2*tvec.array())  ) / (temp.array()*temp.array()*temp.array()*temp.array());
    }
//----------------
    template <class t_vector_type, class p_vector_type>
    static auto ineq_constraint_vectorisedNonIdentical(const t_vector_type& tvec, const p_vector_type& p) {
        return tvec.array() + log(1 - p.array());
    }
};

// VectorOfIIDCGF and VectorOfNonIdenticalCGF classes uses these template vectorised methods to obtain a standard CGF template version
using GeometricCGF = VectorOfIIDCGF<GeometricVectorisedIID_template_version>;
using GeometricNonIdenticalCGF = VectorOfNonIdenticalCGF<GeometricVectorisedNonIdentical_template_version>;


} // namespace CGFs_via_templates

namespace CGFs_with_AD {

class GeometricExtractorIID {
// The templated GeometricCGF class expects only p as its argument
// This Extractor class assumes that p is the only scalar and therefore stored at positon 0
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return parameter_vector[0];
    }
    vec pack(double p) const {
        vec parameter_vector(1);
        parameter_vector << p;
        return parameter_vector;
    }
    a_vector pack(a_scalar p) const {
        a_vector parameter_vector(1);
        parameter_vector << p;
        return parameter_vector;
    }
};

class GeometricExtractorNonIdentical {
// The template class expects the the arguments in the order (tvec, p1, p2,...)
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return parameter_vector;
    }
    vec pack(vec p) const {
        vec parameter_vector(p.size());
        parameter_vector << p;
        return parameter_vector;
    }
    a_vector pack(a_vector p) const {
        a_vector parameter_vector(p.size());
        parameter_vector << p;
        return parameter_vector;
    }
};

using GeometricCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::GeometricCGF, GeometricExtractorIID>>;
using GeometricNonIdenticalCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::GeometricNonIdenticalCGF, GeometricExtractorNonIdentical>>;

class GeometricModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::GeometricCGF, const ScalarAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::GeometricCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::GeometricCGF, const ScalarAdaptor*>;
public:
    GeometricModelCGF(const ScalarAdaptor* p_adaptor) : Base(TemplateBaseCGF(), p_adaptor) {}

    GeometricModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};

class GeometricNonIdenticalModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::GeometricNonIdenticalCGF, const VectorAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::GeometricNonIdenticalCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::GeometricNonIdenticalCGF, const VectorAdaptor*>;
public:
    GeometricNonIdenticalModelCGF(const VectorAdaptor* p_adaptor) : Base(TemplateBaseCGF(), p_adaptor) {}

    GeometricNonIdenticalModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};

} // namespace CGFs_with_AD

} // namespace saddlepoint


#endif // GEOMETRICCGF_H_INCLUDED
