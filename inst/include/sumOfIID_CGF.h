
#ifndef SUMOFIID_CGF_H
#define SUMOFIID_CGF_H

#include "CGF_Defaults.h"
#include <utility>
#include "baseCGF.h"
#include "saddlepoint_types.h"
#include "BaseWrapper.h"
#include "parametric_submodelCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {


#include "needs_full_list_of_CGF_methods.h"
template <class BaseCGF>
class SumOfIID_CGF : public CGF_Defaults<SumOfIID_CGF<BaseCGF>>, private BaseWrapper<BaseCGF> {
private:
    using BaseWrapper<BaseCGF>::base_cgf;
public:
    using BaseWrapper<BaseCGF>::BaseWrapper; // Provides access to constructors for BaseCGF
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        // for n = 2, i.e., Y = X1 + X2 
        // K_Y = n * K_X1
        return n * base_cgf()->K(tvec, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        return (n * base_cgf()->K1(tvec, std::forward<ParamTypes>(other_params)...)).eval();
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        return (n * base_cgf()->K2(tvec, std::forward<ParamTypes>(other_params)...)).eval();
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto tilting_exponent(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        //return n * (base_cgf()->K(tvec, std::forward<ParamTypes>(other_params)...) - (tvec.transpose() * base_cgf()->K1(tvec, std::forward<ParamTypes>(other_params)...)) );
        return n * base_cgf()->tilting_exponent(tvec, std::forward<ParamTypes>(other_params)...) ;
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto neg_ll(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        auto K2_val = K2(tvec, n, std::forward<ParamTypes>(other_params)...);
        return (0.5 * (log(K2_val.determinant()) + tvec.size()*log(2*M_PI))) - tilting_exponent(tvec, n, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class n_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, 
                    const n_type& n, ParamTypes&&... other_params) const {
        return n * base_cgf()->K3operator(tvec, v1, v2, v3, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class n_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const v4_type& v4, 
                    const n_type& n, ParamTypes&&... other_params) const {
        return n * base_cgf()->K4operator(tvec, v1, v2, v3, v4, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class x_type, class y_type, class n_type, class... ParamTypes>
    auto K2operator(const t_vector_type& tvec, const x_type& x, const y_type& y,
                    const n_type& n, ParamTypes&&... other_params) const {
        return n * base_cgf()->K2operator(tvec, x, y, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class B_matrix_type, class n_type, class... ParamTypes>
    auto K2operatorAK2AT(const t_vector_type& tvec, const B_matrix_type& B,
                         const n_type& n, ParamTypes&&... other_params) const {
        return (n * base_cgf()->K2operatorAK2AT(tvec, B, std::forward<ParamTypes>(other_params)...)).eval();
    }
    //---------------------
    // For the factored forms where Q = B D B^T and D has diagonal vector d. Here Q = K2_Y^{-1} = (n * K2_X1)^{-1}
    // The n has an effect on the vector d, but the matrix B can remain unchanged.
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class B3_type, class d3_type, class n_type, class... ParamTypes>
    auto K3K3operatorAABBCC_factored(const t_vector_type& tvec,
                                     const B1_type& B1, d1_type&& d1,
                                     const B2_type& B2, d2_type&& d2,
                                     const B3_type& B3, d3_type&& d3,
                                     const n_type& n, ParamTypes&&... other_params) const {

        return (base_cgf()->K3K3operatorAABBCC_factored(tvec,
                                                        B1, n*std::forward<d1_type>(d1),
                                                        B2, n*std::forward<d2_type>(d2),
                                                        B3, n*std::forward<d3_type>(d3),
                                                        std::forward<ParamTypes>(other_params)...) )/n ;
    }
    //---------------------
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class B3_type, class d3_type, class n_type, class... ParamTypes>
    auto K3K3operatorABCABC_factored(const t_vector_type& tvec,
                                     const B1_type& B1, d1_type&& d1,
                                     const B2_type& B2, d2_type&& d2,
                                     const B3_type& B3, d3_type&& d3,
                                     const n_type& n, ParamTypes&&... other_params) const {
        return (base_cgf()->K3K3operatorABCABC_factored(tvec,
                                                        B1, n*std::forward<d1_type>(d1),
                                                        B2, n*std::forward<d2_type>(d2),
                                                        B3, n*std::forward<d3_type>(d3),
                                                        std::forward<ParamTypes>(other_params)...) )/n;
    }
    //---------------------
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class n_type, class... ParamTypes>
    auto K4operatorAABB_factored(const t_vector_type& tvec,
                                 const B1_type& B1, d1_type&& d1,
                                 const B2_type& B2, d2_type&& d2,
                                 const n_type& n, ParamTypes&&... other_params) const {
        return (base_cgf()->K4operatorAABB_factored(tvec,
                                                   B1, n*std::forward<d1_type>(d1),
                                                   B2, n*std::forward<d2_type>(d2),
                                                   std::forward<ParamTypes>(other_params)...) )/n;
    }
    //---------------------
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class n_type, class... ParamTypes>
    auto K3K3operatorAABBCC(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2, const Q3_type& Q3,
                            const n_type& n, ParamTypes&&... other_params) const {
        return n * n * base_cgf()->K3K3operatorAABBCC(tvec, Q1, Q2, Q3, std::forward<ParamTypes>(other_params)...);
    }
    //---------------------
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class n_type, class... ParamTypes>
    auto K3K3operatorABCABC(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2, const Q3_type& Q3,
                            const n_type& n, ParamTypes&&... other_params) const {
        return n * n * base_cgf()->K3K3operatorABCABC(tvec, Q1, Q2, Q3, std::forward<ParamTypes>(other_params)...) ;
    }
    //---------------------
    template <class t_vector_type, class Q1_type, class Q2_type, class n_type, class... ParamTypes>
    auto K4operatorAABB(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2,
                        const n_type& n, ParamTypes&&... other_params) const {
        return n * base_cgf()->K4operatorAABB(tvec, Q1, Q2, std::forward<ParamTypes>(other_params)...) ;
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto func_T(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        return (base_cgf()->func_T(tvec, std::forward<ParamTypes>(other_params)...))/n;
    }
    //---------------------
    template <class t_vector_type, class n_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const n_type& n, ParamTypes&&... other_params) const {
        // ineq_constraint is independent of n, therefore, we directly call the ineq_constraint
        // function of base_cgf() with tvec and other_params.
        // i.e., Y = X1 + X2 + X3, the inequality constraints on Y are the same as those on X1, X2, X3.
        return base_cgf()->ineq_constraint(tvec, std::forward<ParamTypes>(other_params)...);
    }


};

} // namespace CGFs_via_templates

namespace CGFs_with_AD {


// Function object to provide a saved (n) either as double or AD version
class Provide_n_1CFOwAD {
private:
    double n;
public:
    explicit Provide_n_1CFOwAD(const double& n_) : n(n_) {}
    
    Provide_n_1CFOwAD() = delete;

    
    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const vec& parameter_vector) const {
        return co(n, parameter_vector);
    }
    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const a_vector& parameter_vector) const {
        a_scalar a_n = n;
        return co(a_n, parameter_vector);
    }

};


using SumOfIID_CGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::AdaptedParametersCallbackCGF<saddlepoint::CGFs_via_templates::SumOfIID_CGF<CGF_with_AD*>, Provide_n_1CFOwAD>>;

} // namespace CGFs_with_AD



} // namespace saddlepoint

#endif // SUMOFIID_CGF_H
