#ifndef POISSONCGF_H_INCLUDED
#define POISSONCGF_H_INCLUDED

#include "baseCGF.h"
#include "extractor.h"
#include "saddlepoint_types.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

class PoissonVectorisedIID_template_version : public ScalarVectorisedIID_Defaults<PoissonVectorisedIID_template_version> {

public:
    // These are vectorised functions of K and its derivatives
    // The functions returns a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // lambda as argument remain constant across all values of tvec
    template <class t_vector_type, class lambda_type>
    static auto K_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return -lambda + lambda*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_type>
    static auto K1_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return lambda*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_type>
    static auto K2_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return lambda*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_type>
    static auto K3_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return lambda*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_type>
    static auto K4_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda)  {
        return lambda*exp(tvec.array());
    }
};
//----------------
//----------------
//----------------
class PoissonVectorisedNonIdentical_template_version : public ScalarVectorisedNonIdentical_Defaults<PoissonVectorisedNonIdentical_template_version> {

public:
    // The NonIdentical versions of the K and its derivatives take in the parameter lambda as a vector.
    // Each value from the vector tvec has its own lambda.
    // tvec.size() == lambda.size()
    template <class t_vector_type, class lambda_vector_type>
    static auto K_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return -lambda.array() + lambda.array()*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_vector_type>
    static auto K1_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return lambda.array()*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_vector_type>
    static auto K2_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return lambda.array()*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_vector_type>
    static auto K3_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return lambda.array()*exp(tvec.array());
    }
//----------------
    template <class t_vector_type, class lambda_vector_type>
    static auto K4_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return lambda.array()*exp(tvec.array());
    }
//----------------
};

// We use the VectorOfIIDCGF and VectorOfNonIdenticalCGF to converts the Scalar vectorised functions above to obtain a standard CGF template version
using PoissonCGF = VectorOfIIDCGF<PoissonVectorisedIID_template_version>;
using PoissonNonIdenticalCGF = VectorOfNonIdenticalCGF<PoissonVectorisedNonIdentical_template_version>;


}





namespace CGFs_with_AD {

class PoissonExtractorIID {
// The templated Poisson class expects only lambda as its argument
// This Extractor class assumes that lambda is the only scalar and therefore stored at positon 0
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return parameter_vector[0];
    }
    vec pack(double lambda) const {
        vec parameter_vector(1);
        parameter_vector << lambda;
        return parameter_vector;
    }
    a_vector pack(a_scalar lambda) const {
        a_vector parameter_vector(1);
        parameter_vector << lambda;
        return parameter_vector;
    }
};

class PoissonExtractorNonIdentical {
// The template class expects the the arguments in the order (tvec, lambda1, lambda2,...)
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return parameter_vector;
    }
    vec pack(vec lambda) const {
        vec parameter_vector(lambda.size());
        parameter_vector << lambda;
        return parameter_vector;
    }
    a_vector pack(a_vector lambda) const {
        a_vector parameter_vector(lambda.size());
        parameter_vector << lambda;
        return parameter_vector;
    }
};

using PoissonNonIdenticalCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::PoissonNonIdenticalCGF, PoissonExtractorNonIdentical>>;
using PoissonCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::PoissonCGF, PoissonExtractorIID>>;

class PoissonModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::PoissonCGF, const ScalarAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::PoissonCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::PoissonCGF, const ScalarAdaptor*>;
public:
    PoissonModelCGF(const ScalarAdaptor* lambda_adaptor) : Base(TemplateBaseCGF(), lambda_adaptor) {}

    PoissonModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};

} // namespace CGFs_with_AD
} // namespace saddlepoint


#endif // POISSONCGF_H_INCLUDED
