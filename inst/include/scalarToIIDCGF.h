#ifndef SCALARTOIIDCGF_H_INCLUDED
#define SCALARTOIIDCGF_H_INCLUDED

#include "baseCGF.h"
#include "cgf_from_scalarCGF.h"
#include "scalarCGF.h"
#include "BaseWrapper.h"

namespace saddlepoint {
namespace CGFs_via_templates {

// This class is used to convert a ScalarCGF into a VectorOfIIDCGF
// The size of 'tvec' for all the methods in the template parameter ScalarCGF is 1
// We first create vectorised_iid methods using ScalarCGF and then convert them to a standard CGF using VectorOfIIDCGF class.

template <class ScalarCGF>
class ScalarCGF2VectorisedIID_template_version : public ScalarVectorisedIID_Defaults<ScalarCGF2VectorisedIID_template_version<ScalarCGF>>, protected BaseWrapper<ScalarCGF> {
private:
    using BaseWrapper<ScalarCGF>::base_cgf;
public:
    template <class ScalarCGF_Type>
    ScalarCGF2VectorisedIID_template_version(ScalarCGF_Type&& scgf) : BaseWrapper<ScalarCGF>(scgf) {}

//    using BaseWrapper<ScalarCGF>::BaseWrapper;
//----------------
    template <class t_vector_type, class... ParamTypes>
    auto K_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        vector_type res(tvec.size());
        typedef decltype(tvec.size()) index_type;
        for (index_type i = 0; i < tvec.size(); i++) res(i) = base_cgf()->K(tvec.segment(i,1), std::forward<ParamTypes>(other_params)...);
        return res;
    }
//----------------
    template <class t_vector_type, class... ParamTypes>
    auto K1_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        vector_type res(tvec.size());
        typedef decltype(tvec.size()) index_type;
        for (index_type i = 0; i < tvec.size(); i++) res(i) = (base_cgf()->K1(tvec.segment(i,1), std::forward<ParamTypes>(other_params)...)).value();
        return res;
    }
//----------------
    template <class t_vector_type, class... ParamTypes>
    auto K2_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        vector_type res(tvec.size());
        typedef decltype(tvec.size()) index_type;
        for (index_type i = 0; i < tvec.size(); i++) res(i) = (base_cgf()->K2(tvec.segment(i,1), std::forward<ParamTypes>(other_params)...)).value();
        return res;
    }
//----------------
    template <class t_vector_type, class... ParamTypes>
    auto K3_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        vector_type res(tvec.size());
        vector_type v1(1); v1 << 1; // v1: vector of size 1
        typedef decltype(tvec.size()) index_type;
        for (index_type i = 0; i < tvec.size(); i++) res(i) = base_cgf()->K3operator(tvec.segment(i,1), v1, v1, v1, std::forward<ParamTypes>(other_params)...);
        return res;
    }
//----------------
    template <class t_vector_type, class... ParamTypes>
    auto K4_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        vector_type res(tvec.size());
        vector_type v1(1); v1 << 1; // v1: vector of size 1
        typedef decltype(tvec.size()) index_type;
        for (index_type i = 0; i < tvec.size(); i++) res(i) = base_cgf()->K4operator(tvec.segment(i,1), v1, v1, v1, v1, std::forward<ParamTypes>(other_params)...);
        return res;
    }
    template <class t_vector_type, class... ParamTypes>
    auto ineq_constraint_vectorised_iid(const t_vector_type& tvec, ParamTypes&&... other_params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;

        typedef decltype(tvec.size()) index_type;

        auto first_ineq_vec = base_cgf()->ineq_constraint(tvec.segment(0,1), std::forward<ParamTypes>(other_params)...);
        auto k = first_ineq_vec.size();
        vector_type res(tvec.size()*k);

        if(k > 0){
            for (index_type i = 0; i < tvec.size(); i++) {
                res.segment(i*k,k) = base_cgf()->ineq_constraint(tvec.segment(i,1), std::forward<ParamTypes>(other_params)...);
            }
        }
        return res;
    }


};


template <class ScalarCGF>
using ScalarToIIDCGF = VectorOfIIDCGF<ScalarCGF2VectorisedIID_template_version<ScalarCGF>>;

} // namespace CGFs_via_templates





namespace CGFs_with_AD {

using ScalarToIIDCGF = CGF_with_AD_from_template< CGFs_via_templates::ScalarToIIDCGF<CGF_with_AD*> >;

} // namespace CGFs_with_AD
} // namespace saddlepoint

#endif // SCALARTOIIDCGF_H_INCLUDED

//// TEST CODE
//// We create an instance of the BinomialCGF class, which is an IID class.
// CGF_with_AD* binomial_cgf = new BinomialCGF();
//
//// We then create an instance of the ScalarToIIDCGF class, passing the binomial_cgf instance as an argument.
// CGF_with_AD* scalar2iid_cgf = new ScalarToIIDCGF(binomial_cgf);
//
//// We create two Eigen vectors: tvec and parameter_vector, to be used as inputs for the K2() function.
// Eigen::VectorXd tvec(3), parameter_vector(2);
// tvec << 0.34, -0.45, 0.02;
// parameter_vector << 10, 0.74;
//
//// We call the K2() function on the binomial_cgf and scalar2iid_cgf instances, passing in the tvec and parameter_vector.
//// We expect to get the same results from both calls since binomial_cgf is an IID class.
// std::cout << "Result from binomial_cgf: " << binomial_cgf->K2(tvec, parameter_vector) << std::endl;
// std::cout << "Result from scalar2iid_cgf: " << scalar2iid_cgf->K2(tvec, parameter_vector) << std::endl;
