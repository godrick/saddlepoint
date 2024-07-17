#ifndef BINOMIALCGF_H_INCLUDED
#define BINOMIALCGF_H_INCLUDED



#include "baseCGF.h"
#include "extractor.h"
#include "saddlepoint_types.h"
#include "parametric_submodelCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"


namespace saddlepoint {
namespace CGFs_via_templates {

class BinomialVectorisedIID_template_version : public ScalarVectorisedIID_Defaults<BinomialVectorisedIID_template_version> {

public:
    // These are vectorised functions of K and its derivatives
    // The functions returns a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // N and p as the Binomial parameters remain constant across all values of tvec
    template <class t_vector_type, class N_type, class p_type>
    static auto K_vectorised_iid(const t_vector_type& tvec, const N_type& N, const p_type& p) {
        return N*log(1 - p + p*exp(tvec.array()) ) ;
    }
//----------------
    template <class t_vector_type, class N_type, class p_type>
    static auto K1_vectorised_iid(const t_vector_type& tvec, const N_type& N, const p_type& p) {
        return N*p*exp(tvec.array()) / (1 - p + p*exp(tvec.array()) );
    }
//----------------
    template <class t_vector_type, class N_type, class p_type>
    static auto K2_vectorised_iid(const t_vector_type& tvec, const N_type& N, const p_type& p) {
        auto temp0 = 1 - p + p*exp(tvec.array());
        return N*(1-p)*p*exp(tvec.array()) / ( temp0.array() * temp0.array() );
    }
//----------------
    template <class t_vector_type, class N_type, class p_type>
    static auto K3_vectorised_iid(const t_vector_type& tvec, const N_type& N, const p_type& p) {
        auto temp0 = 1 - p + p*exp(tvec.array());
        return N*(1-p)*p*exp(tvec.array())*(1 - p - p*exp(tvec.array()) ) / (temp0.array()*temp0.array()*temp0.array());
    }
//----------------
    template <class t_vector_type, class N_type, class p_type>
    static auto K4_vectorised_iid(const t_vector_type& tvec, const N_type& N, const p_type& p)  {
        auto temp0 = p*exp(tvec.array());
        auto temp1 = 1-p;
        auto temp2 = temp1 + p*exp(tvec.array());
        return N*temp1*temp0.array()*(4*p*temp0.array() + (temp0.array()*temp0.array()) - 4*temp0 + (temp1*temp1) ) / (temp2.array()*temp2.array()*temp2.array()*temp2.array());
    }
};
//----------------
//----------------
//----------------
class BinomialVectorisedNonIdentical_template_version : public ScalarVectorisedNonIdentical_Defaults<BinomialVectorisedNonIdentical_template_version> {

public:
    // The NonIdentical versions of the K and its derivatives take in the parameters N and p as vectors.
    // Therefore each scalar value from the vector tvec has its own corresponding N and p.
    // tvec.size() == N.size() == p.size()
    template <class t_vector_type, class N_vector_type, class p_vector_type>
    static auto K_vectorisedNonIdentical(const t_vector_type& tvec, const N_vector_type& N, const p_vector_type& p) {
        return N.array()*log(1 - p.array() + p.array()*exp(tvec.array()) ) ;
    }
//----------------
    template <class t_vector_type, class N_vector_type, class p_vector_type>
    static auto K1_vectorisedNonIdentical(const t_vector_type& tvec, const N_vector_type& N, const p_vector_type& p) {
        return N.array()*p.array()*exp(tvec.array()) / (1 - p.array() + p.array()*exp(tvec.array()) );
    }
//----------------
    template <class t_vector_type, class N_vector_type, class p_vector_type>
    static auto K2_vectorisedNonIdentical(const t_vector_type& tvec, const N_vector_type& N, const p_vector_type& p) {
        auto temp0 = (1 - p.array() + p.array()*exp(tvec.array()));
        return N.array()*(1-p.array())*p.array()*exp(tvec.array()) / (temp0.array()*temp0.array());
    }
//----------------
    template <class t_vector_type, class N_vector_type, class p_vector_type>
    static auto K3_vectorisedNonIdentical(const t_vector_type& tvec, const N_vector_type& N, const p_vector_type& p) {
        auto temp0 = (1 - p.array() + p.array()*exp(tvec.array()));
        return N.array()*(1-p.array())*p.array()*exp(tvec.array())*(1 - p.array() - p.array()*exp(tvec.array()) ) / (temp0.array()*temp0.array()*temp0.array()) ;
    }
//----------------
    template <class t_vector_type, class N_vector_type, class p_vector_type>
    static auto K4_vectorisedNonIdentical(const t_vector_type& tvec, const N_vector_type& N, const p_vector_type& p) {
        auto temp0 = p.array()*exp(tvec.array());
        auto temp1 = 1-p.array();
        auto temp2 = temp1.array() + temp0.array();
        return N.array()*temp1.array()*temp0.array()*(4*p.array()*temp0.array() + (temp0.array()*temp0.array()) - 4*temp0.array() + (temp1.array()*temp1.array()) ) / (temp2.array()*temp2.array()*temp2.array()*temp2.array());
    }
//----------------
};

// We use the VectorOfIIDCGF and VectorOfNonIdenticalCGF to convert our Scalar Binomial vectorised functions defined above to obtain a standard binomialCGF template version
using BinomialCGF = VectorOfIIDCGF<BinomialVectorisedIID_template_version>;
using BinomialNonIdenticalCGF = VectorOfNonIdenticalCGF<BinomialVectorisedNonIdentical_template_version>;


}





namespace CGFs_with_AD {

class BinomialExtractorIID {
// The templated Binomial class expects its arguments in the order tvec, N, p
// The Extractor class assumes that N and p are stored in a single vector with N in position 0 and p in position 1
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        return std::make_pair(parameter_vector[0], parameter_vector[1] );
    }
    vec pack(double N, double p) const {
        vec parameter_vector(2);
        parameter_vector << N, p;
        return parameter_vector;
    }
    a_vector pack(a_scalar N, a_scalar p) const {
        a_vector parameter_vector(2);
        parameter_vector << N, p;
        return parameter_vector;
    }
};

class BinomialExtractorNonIdentical {
// The template class expects the the arguments in the order (tvec, N1, N2, N3,...,p1,p2,p3,...)
// The size of Ni's and pi's should be equal
// Also the size of tvec should be equal to the size of Ni's or pi's
// TO DO: Decide where the error checking for the above conditions are to be performed.
public:
    template <class vector_type>
    auto operator()(const vector_type& parameter_vector) const {
        // TO DO: Error checking - parameter_vector.size() is even
                    // DELETE: The following should work but they give an incorrect result. Why could this be happening? The file will compile without any problem but it causes a problem in "parametric_submodelCGF.h" line 42
                    // decltype(parameter_vector.size()) num_samples = parameter_vector.size() / 2;
                    // vector_type N = parameter_vector.head(num_samples);
                    // vector_type p = parameter_vector.tail(num_samples);
                    // return std::make_pair(N, p);
        return std::make_pair(parameter_vector.head(parameter_vector.size()/2), parameter_vector.tail(parameter_vector.size()/2) );
    }
    vec pack(vec N, vec p) const {
        vec parameter_vector(N.size() + p.size());
        parameter_vector << N, p;
        return parameter_vector;
    }
    a_vector pack(a_vector N, a_vector p) const {
        a_vector parameter_vector(N.size() + p.size());
        parameter_vector << N, p;
        return parameter_vector;
    }
};

using BinomialNonIdenticalCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::CGFwithExtractor<saddlepoint::CGFs_via_templates::BinomialNonIdenticalCGF, BinomialExtractorNonIdentical>>;
using BinomialCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::CGFwithExtractor<saddlepoint::CGFs_via_templates::BinomialCGF, BinomialExtractorIID>>;

class BinomialModelCGF : public SeparateParametersCGF<saddlepoint::CGFs_via_templates::BinomialCGF, const ScalarAdaptor*, const ScalarAdaptor*> {
private:
    using TemplateBaseCGF = saddlepoint::CGFs_via_templates::BinomialCGF;
    using Base = SeparateParametersCGF<saddlepoint::CGFs_via_templates::BinomialCGF, const ScalarAdaptor*, const ScalarAdaptor*>;
public:
    BinomialModelCGF(const ScalarAdaptor* n_adaptor, const ScalarAdaptor* prob_adaptor): Base(TemplateBaseCGF(), n_adaptor, prob_adaptor) {}

    BinomialModelCGF() = delete; // otherwise adaptor pointers are uninitialised
};

} // namespace CGFs_with_AD
} // namespace saddlepoint



#endif // BINOMIALCGF_H_INCLUDED
