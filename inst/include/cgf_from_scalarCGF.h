#ifndef CGF_FROM_SCALARCGF_H_INCLUDED
#define CGF_FROM_SCALARCGF_H_INCLUDED

#include <utility>
#include "CGF_Defaults.h"
#include "BaseWrapper.h"

namespace saddlepoint {
namespace CGFs_via_templates {

// Converts a VectorisedScalarCGF into the CGF for the random vector whose entries are i.i.d. copies of the scalar distribution
// The template parameter VectorisedScalarCGF should follow conventions similar to the ScalarCGF_baseIID abstract base class
// As implemented here, VectorisedScalarCGF may in general accept multiple parameter arguments
template <class VectorisedScalarCGF>
class VectorOfIIDCGF : protected BaseWrapper<VectorisedScalarCGF>, public CGF_Defaults<VectorOfIIDCGF<VectorisedScalarCGF>> {
private:
    typedef CGF_Defaults<VectorOfIIDCGF<VectorisedScalarCGF>> Defaults;
    using BaseWrapper<VectorisedScalarCGF>::base_cgf;
public:
    using BaseWrapper<VectorisedScalarCGF>::BaseWrapper;

    template <class t_vector_type, class... ParameterTypes>
    auto K(t_vector_type&& tvec, ParameterTypes&&... params) const {
        return base_cgf()->K_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    // Note: throughout, tvec and params are forwarded
    template <class t_vector_type, class... ParameterTypes>
    auto K1(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->K1_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto K2(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
                    // Eigen::SparseMatrix<scalar_type> diagonal_sparse_matrix(tvec.size(), tvec.size());
                    // diagonal_sparse_matrix.diagonal() = base_cgf()->K2_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
        return (base_cgf()->K2_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...)).matrix().asDiagonal();
    }


    template <class t_vector_type, class... ParameterTypes>
    auto tilting_exponent(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->tilting_exponent_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto neg_ll(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->neg_ll_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto func_T(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->func_T_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }

    using Defaults::K2operator;
    using Defaults::K2operatorAK2AT;
    // These should already be efficient since K2 is returned as a diagonal-typed matrix
    // Previous version:
//    template <class t_vector_type, class vector_type, class... ParameterTypes>
//    auto K2operator(t_vector_type&& tvec, const vector_type& x, const vector_type& y, ParameterTypes&&... params) const {
//        // Since K2 is diagonal, the quadratic form x^T K2 y is the sum of the entrywise products
//        return ( base_cgf()->K2_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * x.array() * y.array() ).sum();
//    }
//    template <class t_vector_type, class... ParameterTypes, class matrix_type>
//    auto K2operatorAK2AT(t_vector_type&& tvec, const matrix_type& A, ParameterTypes&&... params) const {
//        // Since K2 should return its result as a diagonal matrix class, the product A K2 A^T should already be efficiently computed
//        return A * K2(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...) * A.transpose();
//    }

    template <class t_vector_type, class vector_type, class... ParameterTypes>
    auto K3operator(t_vector_type&& tvec, const vector_type& v1, const vector_type& v2, const vector_type& v3,
                    ParameterTypes&&... params) const
    {
        // K3 is a diagonal tensor so the result is a sum of entrywise products
        return ( base_cgf()->K3_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * v1.array() * v2.array() * v3.array() ).sum();
    }
    template <class t_vector_type, class vector_type, class... ParameterTypes>
    auto K4operator(t_vector_type&& tvec, const vector_type& v1, const vector_type& v2, const vector_type& v3, const vector_type& v4,
                    ParameterTypes&&... params) const
    {
        // K4 is a diagonal tensor so the result is a sum of entrywise products
        return ( base_cgf()->K4_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * v1.array() * v2.array() * v3.array() * v4.array() ).sum();
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K4operatorAABB(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, ParameterTypes&&... params) const
    {
        // K4 is a diagonal tensor so the result only uses the diagonals of Q1,Q2
        return ( base_cgf()->K4_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * Q1.diagonal().array() * Q2.diagonal().array() ).sum();
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K3K3operatorAABBCC(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                            ParameterTypes&&... params) const
    {
        // K3(i,j,k) vanishes except when i==j==k so
        // sum_{i1,...,i6} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
        // reduces to sum_{i,j} K3(i,i,i)*K3(j,j,j)*Q1(i,i)*Q2(i,j)*Q3(j,j)
        // The dependence on Q1 and Q3 is through diagonals only
        auto k3val = base_cgf()->K3_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).eval();
        return ( (Q1.diagonal().array()*k3val.array() ).matrix().asDiagonal() * Q2 * (Q3.diagonal().array()*k3val.array() ).matrix().asDiagonal() ).sum() ;
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K3K3operatorABCABC(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                            ParameterTypes&&... params) const
    {
        // K3(i,j,k) vanishes except when i==j==k so
        // sum_{i1,...,i6} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
        // reduces to sum_{i,j} K3(i,i,i)*K3(j,j,j)*Q1(i,j)*Q2(i,j)*Q3(i,j)
        auto k3val = base_cgf()->K3_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).eval();
        return ( k3val.matrix().asDiagonal() * (Q1.array()*Q2.array()*Q3.array()).matrix() * k3val.matrix().asDiagonal() ).sum() ;
    }
    // Not currently implementing other variants of K4/K3K3 operators; default versions used
    template <class t_vector_type, class... ParameterTypes>
    auto ineq_constraint(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->ineq_constraint_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
    }

};

template <class VectorisedScalarCGF>
class VectorOfNonIdenticalCGF : protected BaseWrapper<VectorisedScalarCGF>, public CGF_Defaults<VectorOfNonIdenticalCGF<VectorisedScalarCGF>> {
private:
    typedef CGF_Defaults<VectorOfNonIdenticalCGF<VectorisedScalarCGF>> Defaults;
    using BaseWrapper<VectorisedScalarCGF>::base_cgf;
public:
    using BaseWrapper<VectorisedScalarCGF>::BaseWrapper;

    template <class t_vector_type, class... ParameterTypes>
    auto K(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->K_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    // Note: throughout, tvec and params are forwarded
    template <class t_vector_type, class... ParameterTypes>
    auto K1(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->K1_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto K2(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return (base_cgf()->K2_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...)).matrix().asDiagonal();
    }

    template <class t_vector_type, class... ParameterTypes>
    auto tilting_exponent(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->tilting_exponent_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto neg_ll(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->neg_ll_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto func_T(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->func_T_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).sum();
    }

    using Defaults::K2operator;
    using Defaults::K2operatorAK2AT;
    // These should already be efficient since K2 is returned as a diagonal-typed matrix
    // Previous version:
//    template <class t_vector_type, class vector_type, class... ParameterTypes>
//    auto K2operator(t_vector_type&& tvec, const vector_type& x, const vector_type& y, ParameterTypes&&... params) const {
//        // Since K2 is diagonal, the quadratic form x^T K2 y is the sum of the entrywise products
//        return ( base_cgf()->K2_vectorised_iid(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * x.array() * y.array() ).sum();
//    }
//    template <class t_vector_type, class... ParameterTypes, class matrix_type>
//    auto K2operatorAK2AT(t_vector_type&& tvec, const matrix_type& A, ParameterTypes&&... params) const {
//        // Since K2 should return its result as a diagonal matrix class, the product A K2 A^T should already be efficiently computed
//        return A * K2(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...) * A.transpose();
//    }

    template <class t_vector_type, class vector_type, class... ParameterTypes>
    auto K3operator(t_vector_type&& tvec, const vector_type& v1, const vector_type& v2, const vector_type& v3,
                    ParameterTypes&&... params) const
    {
        // K3 is a diagonal tensor so the result is a sum of entrywise products
        return ( base_cgf()->K3_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * v1.array() * v2.array() * v3.array() ).sum();
    }
    template <class t_vector_type, class vector_type, class... ParameterTypes>
    auto K4operator(t_vector_type&& tvec, const vector_type& v1, const vector_type& v2, const vector_type& v3, const vector_type& v4,
                    ParameterTypes&&... params) const
    {
        // K4 is a diagonal tensor so the result is a sum of entrywise products
        return ( base_cgf()->K4_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * v1.array() * v2.array() * v3.array() * v4.array() ).sum();
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K4operatorAABB(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, ParameterTypes&&... params) const
    {
        // K4 is a diagonal tensor so the result only uses the diagonals of Q1,Q2
        return ( base_cgf()->K4_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).array() * Q1.diagonal().array() * Q2.diagonal().array() ).sum();
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K3K3operatorAABBCC(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                            ParameterTypes&&... params) const
    {
        // K3(i,j,k) vanishes except when i==j==k so
        // sum_{i1,...,i6} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
        // reduces to sum_{i,j} K3(i,i,i)*K3(j,j,j)*Q1(i,i)*Q2(i,j)*Q3(j,j)
        // The dependence on Q1 and Q3 is through diagonals only
        auto k3val = base_cgf()->K3_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...).eval();
        return ( (Q1.diagonal().array()*k3val.array() ).matrix().asDiagonal() * Q2 * (Q3.diagonal().array()*k3val.array() ).matrix().asDiagonal() ).sum() ;
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K3K3operatorABCABC(t_vector_type&& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                            ParameterTypes&&... params) const
    {
        // K3(i,j,k) vanishes except when i==j==k so
        // sum_{i1,...,i6} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
        // reduces to sum_{i,j} K3(i,i,i)*K3(j,j,j)*Q1(i,j)*Q2(i,j)*Q3(i,j)
        auto k3val = base_cgf()->K3_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
        return ( k3val.matrix().asDiagonal() * (Q1.array()*Q2.array()*Q3.array()).matrix() * k3val.matrix().asDiagonal() ).sum() ;
    }
    template <class t_vector_type, class... ParameterTypes>
    auto ineq_constraint(t_vector_type&& tvec, ParameterTypes&&... params) const
    {
        return base_cgf()->ineq_constraint_vectorisedNonIdentical(std::forward<t_vector_type>(tvec), std::forward<ParameterTypes>(params)...);
    }
};


} // namespace CGFs_via_templates
} // namespace saddlepoint

#endif // CGF_FROM_SCALARCGF_H_INCLUDED
