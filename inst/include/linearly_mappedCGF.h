#ifndef LINEARLY_MAPPEDCGF_H_INCLUDED
#define LINEARLY_MAPPEDCGF_H_INCLUDED


#include "CGF_Defaults.h"
#include <utility>
#include "baseCGF.h"
#include "saddlepoint_types.h"
#include "BaseWrapper.h"
#include "parametric_submodelCGF.h"

namespace saddlepoint {

namespace CGFs_via_virtual_functions {
template <class scalar_type, class vector_type, class matrix_type, class A_matrix_type>
class Linearly_Mapped_Scratchpad;
} // namespace CGFs_via_virtual_functions

namespace CGFs_with_AD {

class linearly_mappedCGF_with_AD : public CGF_with_AD {
private:
    CGF_with_AD* base_CGF;
    mat A;
    a_matrix a_A;

    typedef CGF_base<double, vec, mat> D_CGF;
    typedef CGF_base<a_scalar, a_vector, a_matrix> AD_CGF;
    // convenience typedefs

public:
    linearly_mappedCGF_with_AD(CGF_with_AD* b_CGF, const mat& matrix_A) : base_CGF(b_CGF), A(matrix_A), a_A(matrix_A.rows(), matrix_A.cols()) {
        for(int i = 0; i < A.rows(); ++i){
            for(int j = 0; j < A.cols(); ++j){
                a_A(i,j) = A(i,j);
            }
        } // a_A = A.cast<a_scalar>();
    }
//------------------------------------------------------------------------------------------------
// Distribution-specific functions
// These compute transform tvec into t_inner = A.transpose() * tvec and pass t_inner to base_CGF
    // along with appropriate changes to output values based on the matrix A
//------------------------------------------------------------------------------------------------
    double K(const vec& tvec, const vec& parameter_vector) const override {
        return base_CGF->K( (A.transpose() * tvec).eval(), parameter_vector);
    }
    a_scalar K(const a_vector& tvec, const a_vector& parameter_vector) const override {
        return base_CGF->K( (a_A.transpose() * tvec).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    vec K1(const vec& tvec, const vec& parameter_vector) const override {
        return (A * base_CGF->K1( (A.transpose() * tvec).eval(), parameter_vector) ).eval();
    }
    a_vector K1(const a_vector& tvec, const a_vector& parameter_vector) const override {
        return (a_A * base_CGF->K1( (a_A.transpose() * tvec).eval(), parameter_vector)).eval();
    }
//------------------------------------------------------------------------------------------------
    mat K2(const vec& tvec, const vec& parameter_vector) const override {
        return (A * base_CGF->K2( (A.transpose() * tvec).eval(), parameter_vector) * A.transpose() ).eval();
    }
    a_matrix K2(const a_vector& tvec, const a_vector& parameter_vector) const override {
        return (a_A * base_CGF->K2( (a_A.transpose() * tvec).eval(), parameter_vector)  * a_A.transpose()).eval();
    }
//------------------------------------------------------------------------------------------------
    double K2operator(const vec& tvec, const vec& x, const vec& y, const vec& parameter_vector) const override {
        // Since K2 = A base_K2 A^T, this reduces to x^T K2 x = (A^T x)^T base_K2 (A^T x)
        return base_CGF->K2operator( A.transpose() * tvec, A.transpose() * x, A.transpose() * y, parameter_vector);
    }
    a_scalar K2operator(const a_vector& tvec, const a_vector& x, const a_vector& y, const a_vector& parameter_vector) const override {
        return base_CGF->K2operator( a_A.transpose() * tvec, a_A.transpose() * x, a_A.transpose() * y, parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    mat K2operatorAK2AT(const vec& tvec, const mat& B, const vec& parameter_vector) const override {
        // Returns B K2 B^T
        // Since K2 = A base_K2 A^T, this reduces to (B A) base_K2 (B A)^T
        return base_CGF->K2operatorAK2AT( A.transpose() * tvec, B*A, parameter_vector);
    }
    a_matrix K2operatorAK2AT(const a_vector& tvec, const a_matrix& B, const a_vector& parameter_vector) const override {
        return base_CGF->K2operatorAK2AT( a_A.transpose() * tvec, B*a_A, parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    double K3operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& parameter_vector) const override {
        return base_CGF->K3operator( A.transpose() * tvec, A.transpose()*v1, A.transpose()*v2, A.transpose()*v3, parameter_vector);
    }
    a_scalar K3operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& parameter_vector) const override {
        return base_CGF->K3operator( a_A.transpose() * tvec, a_A.transpose()*v1, a_A.transpose()*v2, a_A.transpose()*v3, parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    double K4operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& v4, const vec& parameter_vector) const override {
        return base_CGF->K4operator( A.transpose() * tvec, A.transpose()*v1, A.transpose()*v2, A.transpose()*v3, A.transpose()*v4, parameter_vector);
    }
    a_scalar K4operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& v4, const a_vector& parameter_vector) const override {
        return base_CGF->K4operator( a_A.transpose() * tvec, a_A.transpose()*v1, a_A.transpose()*v2, a_A.transpose()*v3, a_A.transpose()*v4, parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    double K4operatorAABB(const vec& tvec,
                          const mat& Q1,
                          const mat& Q2,
                          const vec& parameter_vector) const override {
        return base_CGF->K4operatorAABB((A.transpose()*tvec).eval(), (A.transpose() * Q1 * A).eval(), (A.transpose() * Q2 * A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    a_scalar K4operatorAABB(const a_vector& tvec,
                            const a_matrix& Q1,
                            const a_matrix& Q2,
                            const a_vector& parameter_vector) const override {
        return base_CGF->K4operatorAABB((a_A.transpose()*tvec).eval(), (a_A.transpose() * Q1 * a_A).eval(), (a_A.transpose() * Q2 * a_A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    double K3K3operatorAABBCC(const vec& tvec,
                              const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {
        return base_CGF->K3K3operatorAABBCC((A.transpose()*tvec).eval(), (A.transpose() * Q1 * A).eval(), (A.transpose() * Q2 * A).eval(), (A.transpose() * Q3 * A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    a_scalar K3K3operatorAABBCC(const a_vector& tvec,
                                const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {
        return base_CGF->K3K3operatorAABBCC((a_A.transpose()*tvec).eval(), (a_A.transpose() * Q1 * a_A).eval(), (a_A.transpose() * Q2 * a_A).eval(), (a_A.transpose() * Q3 * a_A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    double K3K3operatorABCABC(const vec& tvec,
                              const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {
        return base_CGF->K3K3operatorABCABC((A.transpose()*tvec).eval(), (A.transpose() * Q1 * A).eval(), (A.transpose() * Q2 * A).eval(), (A.transpose() * Q3 * A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
    a_scalar K3K3operatorABCABC(const a_vector& tvec,
                                const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {
        return base_CGF->K3K3operatorABCABC((a_A.transpose()*tvec).eval(), (a_A.transpose() * Q1 * a_A).eval(), (a_A.transpose() * Q2 * a_A).eval(), (a_A.transpose() * Q3 * a_A).eval(), parameter_vector);
    }
//------------------------------------------------------------------------------------------------
// Optional overrides
// tilting_exponent(...), and therefore neg_ll(...), can be implemented more directly by noting that with K(t, theta) = K_inner(A.transpose()*t, theta),
    // K(s, theta) - s K'(s, theta) == K_inner(s A, theta) - s A K'_inner(s A, theta)
//------------------------------------------------------------------------------------------------
    double tilting_exponent(const vec& tvec, const vec& parameter_vector) const override {
        return base_CGF->tilting_exponent(A.transpose()*tvec, parameter_vector);
    }
    double neg_ll(const vec& tvec, const vec& parameter_vector) const override {
        double res = linearly_mappedCGF_with_AD::tilting_exponent(tvec, parameter_vector);
        mat K2_val = linearly_mappedCGF_with_AD::K2(tvec, parameter_vector);
        // non-virtual function calls

        res -= (0.5 * (log(K2_val.determinant()) + tvec.size()*log(2*M_PI)));
        return -res;
    }
//------------------------------------------------------------------------------------------------
// Scratchpad methods, implementation below
//------------------------------------------------------------------------------------------------
    CGF_base<double, vec, mat>::Scratchpad* scratchpad(const vec& tvec, const vec& parameter_vector) const override;
    CGF_base<a_scalar, a_vector, a_matrix>::Scratchpad* scratchpad(const a_vector& tvec, const a_vector& parameter_vector) const override;
};

} // namespace CGFs_with_AD

//------------------------------------------------------------------------------------------------
// Scratchpad classes for linearly mapped CGFs
// The value A.transpose()*tvec is computed once and passed to the scratchpad for the inner CGF, but not otherwise needed
// The matrix A (stored by reference) is retained in order to transform the values of K1, K2
    // and the matrices Q in the operator-style functions
//------------------------------------------------------------------------------------------------
namespace CGFs_via_virtual_functions {

template <class scalar_type, class vector_type, class matrix_type, class A_matrix_type>
class Linearly_Mapped_Scratchpad_Q;

template <class scalar_type, class vector_type, class matrix_type, class A_matrix_type>
class Linearly_Mapped_Scratchpad : public CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad {
private:
    typedef CGF_base<scalar_type, vector_type, matrix_type> CGF_type;
    typedef typename CGF_type::Scratchpad Scratchpad_type;
    typedef typename CGF_type::Scratchpad_Q Scratchpad_Q_type;

    const A_matrix_type& A;
    std::unique_ptr<Scratchpad_type> base_scratchpad;

    struct K1_finder {
        vector_type operator()(const A_matrix_type& a, const std::unique_ptr<Scratchpad_type>& bsp) {return a * bsp->K1();}
    };
    SavedResult<vector_type, K1_finder> saved_k1;

    struct K2_finder {
        matrix_type operator()(const A_matrix_type& a, const std::unique_ptr<Scratchpad_type>& bsp) {return a * bsp->K2() * a.transpose();}
    };
    SavedResult<matrix_type, K2_finder> saved_k2;

    struct neg_ll_finder {
        scalar_type operator()(Linearly_Mapped_Scratchpad* sp, const A_matrix_type& a) {
            const scalar_type& te = sp->Linearly_Mapped_Scratchpad::tilting_exponent();
            const matrix_type& k2 = sp->Linearly_Mapped_Scratchpad::K2();
            // Note: bypasses virtual function dispatch mechanism and calls the Linearly_Mapped_Scratchpad versions of these functions,
                // even if *sp is actually a derived type with overridden versions of these virtual functions

            return 0.5 * (log(k2.determinant()) + a.rows() * log(2*M_PI)) - te;
            // size of observed vector can be recovered from A as number of rows
            // Note: the determinant could be more easily computed if K2_val had an LU or related decomposition.
            // This would also share work in computing the matric mat_Q in the funt_T() method
            // TO DO: adapt this method and the related func_T method to share this work.
        }
    };
    SavedResult<scalar_type, neg_ll_finder> saved_neg_ll;

    struct func_T_finder {
        scalar_type operator()(Linearly_Mapped_Scratchpad* sp, const A_matrix_type& a) {
            matrix_type mat_Q = sp->Linearly_Mapped_Scratchpad::K2().inverse();
            std::unique_ptr<Scratchpad_Q_type> p_spq(sp->Linearly_Mapped_Scratchpad::scratchpad_q(mat_Q));
            // Extend sp to its Scratchpad_Q equivalent, with Q == mat_Q == K2.inverse()
            // Note: bypasses virtual function dispatch mechanism and calls the Linearly_Mapped_Scratchpad versions of these functions,
                // even if *sp is actually a derived type with overridden versions of these virtual functions

            return p_spq->K4operatorAABB()/8 - p_spq->K3K3operatorAABBCC()/8 - p_spq->K3K3operatorABCABC()/12;
            // TO DO: See also comments about neg_ll and sharing work with mat_Q
        }
    };
    SavedResult<scalar_type, func_T_finder> saved_func_T;

    friend class Linearly_Mapped_Scratchpad_Q<scalar_type, vector_type, matrix_type, A_matrix_type>;

public:
    Linearly_Mapped_Scratchpad(const vector_type& tvec, const vector_type& parameter_vector, const A_matrix_type& a, const CGF_type* base_cgf)
        : A(a),
        base_scratchpad(base_cgf->scratchpad(a.transpose() * tvec, parameter_vector)) {}
        // base_scratchpad is constructed using the vector s_inner = s A as its CGF argument (for s == tvec.transpose(), in row vector notation)

    Scratchpad_Q_type* scratchpad_q(const matrix_type& Q) override {
        return new Linearly_Mapped_Scratchpad_Q<scalar_type, vector_type, matrix_type, A_matrix_type>(Q, this);
    }

    const scalar_type& K() override {
        return base_scratchpad->K();
    }
    const vector_type& K1() override {return saved_k1(A, base_scratchpad);}
    const matrix_type& K2() override {return saved_k2(A, base_scratchpad);}

    const scalar_type& tilting_exponent() override {
        // Delegates to base_scratchpad because of the identity K(s, theta) - s K'(s, theta) == K_inner(s A, theta) - s A K'_inner(s A, theta)
        return base_scratchpad->tilting_exponent();
    }
    const scalar_type& neg_ll() override {return saved_neg_ll(this, A);}
    const scalar_type& func_T() override {return saved_func_T(this, A);}
};
//------------------------------------------------------------------------------------------------
template <class scalar_type, class vector_type, class matrix_type, class A_matrix_type>
class Linearly_Mapped_Scratchpad_Q : public CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q {
private:
    typedef CGF_base<scalar_type, vector_type, matrix_type> CGF_type;
    typedef typename CGF_type::Scratchpad Scratchpad_type;
    typedef typename CGF_type::Scratchpad_Q Scratchpad_Q_type;
    typedef Linearly_Mapped_Scratchpad<scalar_type, vector_type, matrix_type, A_matrix_type> LMSP_type;

    LMSP_type* plain_sp;
    // non-owned pointer to Scratchpad that created this
    // Note: data member, not base class, to avoid copying already calculated values
    // However this means that the virtual functions already implemented by plain_sp have to be explicitly re-implemented to re-route them to plain_sp

    std::unique_ptr<Scratchpad_Q_type> base_sp_q;
    // owned pointer to Scratchpad_Q for inner CGF

    static matrix_type conjugate(const matrix_type& Q, const A_matrix_type& A) {return A.transpose() * Q * A;}
    // Helper function for use in constructor

public:
    Linearly_Mapped_Scratchpad_Q(const matrix_type& Q, LMSP_type* psp)
        : plain_sp(psp),
          base_sp_q(psp->base_scratchpad->scratchpad_q(conjugate(Q, psp->A))) {}
          // base_sp_q is constructed using the matrix A^T Q A as its square matrix argument
    // Normally this constructor will only be used by Linearly_Mapped_Scratchpad::scratchpad_q

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

    const scalar_type& K4operatorAABB() override {return base_sp_q->K4operatorAABB();}
    const scalar_type& K3K3operatorAABBCC() override {return base_sp_q->K3K3operatorAABBCC();}
    const scalar_type& K3K3operatorABCABC() override {return base_sp_q->K3K3operatorABCABC();}
    // Equivalent because of relations between K3, K4; K3_inner, K4_inner; and the matrix A
};

} // namespace CGFs_via_virtual_functions

namespace CGFs_with_AD {
//------------------------------------------------------------------------------------------------
// Scratchpad methods from linearly_mappedCGF_with_AD
// The type of the scratchpad object returned is determined by the types of the arguments
// The value A.transpose()*tvec is computed once and used to construct a scratchpad for the inner CGF, but not otherwise needed
//------------------------------------------------------------------------------------------------
CGFs_via_virtual_functions::CGF_base<double, vec, mat>::Scratchpad* linearly_mappedCGF_with_AD::scratchpad(const vec& tvec, const vec& parameter_vector) const {
    return new CGFs_via_virtual_functions::Linearly_Mapped_Scratchpad<double, vec, mat, mat>(tvec, parameter_vector, A, static_cast<const D_CGF*>(base_CGF));
}
CGFs_via_virtual_functions::CGF_base<a_scalar, a_vector, a_matrix>::Scratchpad* linearly_mappedCGF_with_AD::scratchpad(const a_vector& tvec, const a_vector& parameter_vector) const {
    return new CGFs_via_virtual_functions::Linearly_Mapped_Scratchpad<a_scalar, a_vector, a_matrix, a_matrix>(tvec, parameter_vector, a_A, static_cast<const AD_CGF*>(base_CGF));
}

} // namespace CGFs_with_AD



namespace CGFs_via_templates {

template <class BaseCGF>
class LinearlyMappedCGF : public CGF_Defaults<LinearlyMappedCGF<BaseCGF>>, private BaseWrapper<BaseCGF> {
    // This class creates a CGF for the random vector Y=AX, where X is the random vector specified in BaseCGF
    // The matrix A is interpreted as a parameter to be supplied to each method
    //
    // Methods take the form K(tvec, A, other_params...) or K3operator(tvec, v1, v2, v3, A, other_params...)
    // where other_params... are the parameters passed to BaseCGF
private:
    using BaseWrapper<BaseCGF>::base_cgf;
    typedef CGF_Defaults<LinearlyMappedCGF<BaseCGF>> Defaults;
public:
    using BaseWrapper<BaseCGF>::BaseWrapper; // Provides access to constructors for BaseCGF

    template <class t_vector_type, class A_matrix_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const A_matrix_type& A, ParamTypes&&... other_params) const {
        // Key identity: K_Y(t) = K_X(A^T t) (with t assumed to be a column vector)
        return base_cgf()->K(A.transpose() * tvec, std::forward<ParamTypes>(other_params)...);
    }

    template <class t_vector_type, class A_matrix_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const A_matrix_type& A, ParamTypes&&... other_params) const {
        // Key identity: K_Y' = A K_X'
        return (A * base_cgf()->K1(A.transpose() * tvec, std::forward<ParamTypes>(other_params)...)).eval();
    }

    template <class t_vector_type, class A_matrix_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const A_matrix_type& A, ParamTypes&&... other_params) const {
        //// Key identity: K_Y'' = A K_X'' A^T
        //return base_cgf()->K2operatorAK2AT((A.transpose() * tvec).eval(), A, std::forward<ParamTypes>(other_params)...).eval();
        return (A * base_cgf()->K2( (A.transpose() * tvec).eval(), std::forward<ParamTypes>(other_params)...) * A.transpose() ).eval();
    }

    template <class t_vector_type, class A_matrix_type, class... ParamTypes>
    auto tilting_exponent(const t_vector_type& tvec, const A_matrix_type& A, ParamTypes&&... other_params) const {
        // Key identity: K_Y(t) - t^T K_Y'(t) = K_X(A^T t) - t^T A K_X'(A^T t) = K_X(A^T t) - (A^T t)^T K_X'(A^T t)
        return base_cgf()->tilting_exponent(A.transpose() * tvec, std::forward<ParamTypes>(other_params)...);
    }

    using Defaults::neg_ll; // optional; this line documents explicitly that we use default behaviour
    using Defaults::func_T;

    template <class t_vector_type, class x_type, class y_type, class A_matrix_type, class... ParamTypes>
    auto K2operator(const t_vector_type& tvec, const x_type& x, const y_type& y,
                    const A_matrix_type& A, ParamTypes&&... other_params) const {
        // Key identity: x^T K_Y'' y = x^T A K_X'' A^T y = (A^T x)^T K_X'' A^T y
        return base_cgf()->K2operator(A.transpose() * tvec, A.transpose()*x, A.transpose()*y, std::forward<ParamTypes>(other_params)...);
    }

    template <class t_vector_type, class B_matrix_type, class A_matrix_type, class... ParamTypes>
    auto K2operatorAK2AT(const t_vector_type& tvec, const B_matrix_type& B,
                         const A_matrix_type& A, ParamTypes&&... other_params) const {
        // Returns B K_Y'' B^T as a function of the supplied (non-parameter) argument B
        // Key identity: B K_Y'' B^T = B A K_X'' A^T B^T = (B A) K_X'' (B A)^T
        return base_cgf()->K2operatorAK2AT(A.transpose() * tvec, B * A, std::forward<ParamTypes>(other_params)...);
    }

    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class A_matrix_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3,
                    const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K3operator(A.transpose() * tvec, A.transpose()*v1, A.transpose()*v2, A.transpose()*v3, std::forward<ParamTypes>(other_params)...);
    }

    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class A_matrix_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const v4_type& v4,
                    const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K4operator(A.transpose() * tvec,
                                  A.transpose()*v1, A.transpose()*v2, A.transpose()*v3, A.transpose()*v4,
                                  std::forward<ParamTypes>(other_params)...);
    }

    // All the operator forms involving matrices Q are equivalent to applying the same method for BaseCGF with Q_inner = A^T Q A
    template <class t_vector_type, class Q1_type, class Q2_type, class A_matrix_type, class... ParamTypes>
    auto K4operatorAABB(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2,
                        const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K4operatorAABB(A.transpose() * tvec, A.transpose()*Q1*A, A.transpose()*Q2*A,
                                      std::forward<ParamTypes>(other_params)...);
    }
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class A_matrix_type, class... ParamTypes>
    auto K3K3operatorAABBCC(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2, const Q3_type& Q3,
                            const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K3K3operatorAABBCC(A.transpose() * tvec, A.transpose()*Q1*A, A.transpose()*Q2*A, A.transpose()*Q3*A,
                                        std::forward<ParamTypes>(other_params)...);
    }
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class A_matrix_type, class... ParamTypes>
    auto K3K3operatorABCABC(const t_vector_type& tvec, const Q1_type& Q1, const Q2_type& Q2, const Q3_type& Q3,
                            const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K3K3operatorABCABC(A.transpose() * tvec, A.transpose()*Q1*A, A.transpose()*Q2*A, A.transpose()*Q3*A,
                                          std::forward<ParamTypes>(other_params)...);
    }

    // For the factored forms where Q = B D B^T and D has diagonal vector d, note that Q_inner = A^T Q A = (A^T B) D (A^T B)^T
    // Note about sizes: if A is n-by-m then B is n-by-r for some r, and A^T B is m-by-r
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class A_matrix_type, class... ParamTypes>
    auto K4operatorAABB_factored(const t_vector_type& tvec, const B1_type& B1, d1_type&& d1, const B2_type& B2, d2_type&& d2,
                                 const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K4operatorAABB_factored(A.transpose() * tvec,
                                                   A.transpose() * B1, std::forward<d1_type>(d1),
                                                   A.transpose() * B2, std::forward<d2_type>(d2),
                                               std::forward<ParamTypes>(other_params)...);
    }
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class B3_type, class d3_type,
                class A_matrix_type, class... ParamTypes>
    auto K3K3operatorAABBCC_factored(const t_vector_type& tvec,
                                     const B1_type& B1, d1_type&& d1, const B2_type& B2, d2_type&& d2, const B3_type& B3, d3_type&& d3,
                                     const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K3K3operatorAABBCC_factored(A.transpose() * tvec,
                                                       A.transpose() * B1, std::forward<d1_type>(d1),
                                                       A.transpose() * B2, std::forward<d2_type>(d2),
                                                       A.transpose() * B3, std::forward<d3_type>(d3),
                                               std::forward<ParamTypes>(other_params)...);
    }
    template <class t_vector_type, class B1_type, class d1_type, class B2_type, class d2_type, class B3_type, class d3_type,
                class A_matrix_type, class... ParamTypes>
    auto K3K3operatorABCABC_factored(const t_vector_type& tvec,
                                     const B1_type& B1, d1_type&& d1, const B2_type& B2, d2_type&& d2, const B3_type& B3, d3_type&& d3,
                                     const A_matrix_type& A, ParamTypes&&... other_params) const {
        return base_cgf()->K3K3operatorABCABC_factored(A.transpose() * tvec,
                                                       A.transpose() * B1, std::forward<d1_type>(d1),
                                                       A.transpose() * B2, std::forward<d2_type>(d2),
                                                       A.transpose() * B3, std::forward<d3_type>(d3),
                                               std::forward<ParamTypes>(other_params)...);
    }
    template <class t_vector_type, class A_matrix_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const A_matrix_type& A, ParamTypes&&... other_params) const {
        // inequality constraints for the transformed variable Y = A * X are the same as those
        // for the original variable X, evaluated at the transformed input A.transpose() * tvec.
        return base_cgf()->ineq_constraint(A.transpose() * tvec, std::forward<ParamTypes>(other_params)...);
    }

}; // class LinearlyMappedCGF

} // namespace CGFs_via_templates

namespace CGFs_with_AD {


// Function object to provide a saved matrix either as double or AD version
// For use with AdaptedParametersCallbackCGF (see above) in conjunction with CGF_with_AD* as the BaseCGF class
class ProvideMatrix1CFOwAD {
private:
    mat A;
    Eigen::SparseMatrix<double> sparse_A;

public:
    explicit ProvideMatrix1CFOwAD(const mat& a) : A(a), sparse_A(a.sparseView()) {}

    ProvideMatrix1CFOwAD() = delete;
    // No default constructor as otherwise A will be initialized as an empty matrix, which is undesirable

    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const vec& parameter_vector) const {
        return co(sparse_A, parameter_vector);
    }
    // If parameter_vector is provided as type vec, pass the sparse matrix A as Eigen::SparseMatrix<double> type

    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const a_vector& parameter_vector) const {
        return co(sparse_A.cast<a_scalar>(), parameter_vector);
    }
    // If parameter_vector is provided as type a_vector, convert sparse matrix A to be of an AD type
};

using LinearlyMappedCGF = CGF_with_AD_from_template<CGFs_via_templates::AdaptedParametersCallbackCGF<CGFs_via_templates::LinearlyMappedCGF<CGF_with_AD*>, ProvideMatrix1CFOwAD>>;
// Note: constructor of the form LinearlyMappedCGF(CGF_with_AD*, const mat&) is provided (via AdaptedParametersCallbackCGF)

} // namespace CGFs_with_AD
} // namespace saddlepoint




#endif // LINEARLY_MAPPEDCGF_H_INCLUDED
