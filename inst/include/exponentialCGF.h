#ifndef BASECGF_CODE_LAYOUT_EXPONENTIALCGF_H
#define BASECGF_CODE_LAYOUT_EXPONENTIALCGF_H

#include "baseCGF.h"
#include "extractor.h"
#include "saddlepoint_types.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

class ExponentialVectorisedIID_template_version : public ScalarVectorisedIID_Defaults<ExponentialVectorisedIID_template_version> {
public:
    template <class t_vector_type, class lambda_type>
    static auto K_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return -log(lambda - tvec.array()) + log(lambda);
    }

    template <class t_vector_type, class lambda_type>
    static auto K1_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return 1 / (lambda - tvec.array());
    }

    template <class t_vector_type, class lambda_type>
    static auto K2_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return 1 / ((lambda - tvec.array()) * (lambda - tvec.array()));
    }

    template <class t_vector_type, class lambda_type>
    static auto K3_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return 2 / ( (lambda - tvec.array()) * (lambda - tvec.array()) * (lambda - tvec.array()) );
    }

    template <class t_vector_type, class lambda_type>
    static auto K4_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return 6 / ((lambda - tvec.array()) * (lambda - tvec.array()) * (lambda - tvec.array()) * (lambda - tvec.array()));
    }

    template <class t_vector_type, class lambda_type>
    static auto ineq_constraint_vectorised_iid(const t_vector_type& tvec, const lambda_type& lambda) {
        return tvec.array() - lambda;
    }
};

class ExponentialVectorisedNonIdentical_template_version : public ScalarVectorisedNonIdentical_Defaults<ExponentialVectorisedNonIdentical_template_version> {

public:
    template <class t_vector_type, class lambda_vector_type>
    static auto K_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return -log(lambda.array() - tvec.array()) + log(lambda.array());
    }

    template <class t_vector_type, class lambda_vector_type>
    static auto K1_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return 1 / (lambda.array() - tvec.array());
    }

    template <class t_vector_type, class lambda_vector_type>
    static auto K2_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return 1 / ((lambda.array() - tvec.array()) * (lambda.array() - tvec.array()));
    }

    template <class t_vector_type, class lambda_vector_type>
    static auto K3_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return 2 / ((lambda.array() - tvec.array()) * (lambda.array() - tvec.array()) * (lambda.array() - tvec.array()));
    }

    template <class t_vector_type, class lambda_vector_type>
    static auto K4_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return 6 / ((lambda.array() - tvec.array()) * (lambda.array() - tvec.array()) * (lambda.array() - tvec.array()) * (lambda.array() - tvec.array()));
    }

    template <class t_vector_type, class lambda_vector_type>
    static auto ineq_constraint_vectorisedNonIdentical(const t_vector_type& tvec, const lambda_vector_type& lambda) {
        return tvec.array() - lambda.array();
    }
};

using ExponentialCGF = VectorOfIIDCGF<ExponentialVectorisedIID_template_version>;
using ExponentialNonIdenticalCGF = VectorOfNonIdenticalCGF<ExponentialVectorisedNonIdentical_template_version>;

} // namespace CGFs_via_templates

namespace CGFs_with_AD {


class ExponentialExtractorIID {
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

class ExponentialExtractorNonIdentical {
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

using ExponentialCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::ExponentialCGF, ExponentialExtractorIID>>;
using ExponentialNonIdenticalCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<saddlepoint::CGFs_via_templates::ExponentialNonIdenticalCGF, ExponentialExtractorNonIdentical>>;

class ExponentialModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::ExponentialCGF, const ScalarAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::ExponentialCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::ExponentialCGF, const ScalarAdaptor*>;
public:
    explicit ExponentialModelCGF(const ScalarAdaptor* lambda_adaptor) : Base(TemplateBaseCGF(), lambda_adaptor) {}

    ExponentialModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};

class ExponentialNonIdenticalModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::ExponentialNonIdenticalCGF, const VectorAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::ExponentialNonIdenticalCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::ExponentialNonIdenticalCGF, const VectorAdaptor*>;
public:
    explicit ExponentialNonIdenticalModelCGF(const VectorAdaptor* lambda_adaptor) : Base(TemplateBaseCGF(), lambda_adaptor) {}

    ExponentialNonIdenticalModelCGF() = delete;
};

}// namespace CGFs_with_AD


} // namespace saddlepoint


#endif //BASECGF_CODE_LAYOUT_EXPONENTIALCGF_H
