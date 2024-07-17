#ifndef BASECGF_H_INCLUDED
#define BASECGF_H_INCLUDED
// A header file that defines CGF_base //
// -------------------------------------
// The notation "CGF" appears throughout the code in different forms.
// We use it in a general sense to refer to a number of functions.
// For instance, CGF_base is defined as a class containing the cumulant
// generating function of a specific distribution/model, first and
// second tvec-gradients of the cumulant generating function and a
// number of operators useful for error computation.

#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>
// #include "saddlepoint_types.h"

#include "needs_full_list_of_CGF_methods.h"

// #include "atomic_funcs.h"

namespace saddlepoint {
namespace CGFs_via_virtual_functions {

template<class scalar_type, class vector_type, class matrix_type>
class CGF_base{
    // CGF_base defines a virtual function interface for the cumulant generating function (CGF) denoted as K,
    // first and second-order tvec-gradients of K, labelled as K1 and K2,
    // and a number of operators ***including a better (faster) implementation of K2 (not yet implemented)***

    // Usage: The user implements specific distributions of interest (or parametric families of distributions)
       // as derived classes inheriting from the CGF_base class. i.e. CGF_base is a base
       // See for instance Multinomial_CGF or Gamma_CGF or ....

    // The arguments to the following functions are
       // tvec : saddlepoint value as a vector of length m, the length of the random variable whose distribution is under consideration.
       // parameter_vector : parameter(s) of the underlying distribution of interest, stored in a single vector
    // Template parameters:
        // scalar_type : the return type for the scalar-valued functions
        // vector_type : a vector type, typically with entries of type scalar_type; the return type for the tvec-gradient K1
        // matrix_type : a matrix type, typically with entries of type scalar_type; the return type for the Hessian K2
        // See class CGF_with_AD below for two examples.
    // Therefore,
       //K(tvec, parameter_vector): is the value of the cumulant generating function at saddlepoint value "tvec" and parameter value "parameter_vector"
       //K1 and K2 are also similar functions whose values are the first and second tvec-gradients of cumulant generating function.
       //K4operatorAABB as an operator takes additional arguments Q1 and Q2 (symmetric) matrices of dimension (m by m)
            // and returns $\sum_{j_1,j_2,j_3,j_4}$ from 1 to m of Q1(j1,j2)*Q2(j3,j4)*K4(j1,j2,j3,j4)
            // where K4 is the fourth-order partial derivative of K wrt. tvec (i.e. forth tvec-gradients of the cumulant generating function)
public:
    virtual ~CGF_base() = default;

//-----------------------------------------------------------
    // Distribution-specific functions must be overridden by derived classes
    virtual scalar_type K(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K1(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual matrix_type K2(const vector_type& tvec, const vector_type& parameter_vector) const = 0;

//-----------------------------------------------------------
    // Default implementations based on calling distribution-specific virtual functions
    // These functions may be overridden if derived classes can implement them more efficiently
    // For instance, in neg_ll, the determinant can be calculated more efficiently if K2 is known to be a diagonal matrix
    virtual scalar_type tilting_exponent(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns the log of the numerator in the saddlepoint approximation
        scalar_type res = -tvec.transpose()*K1(tvec, parameter_vector);
        res += K(tvec, parameter_vector);
        return res;
    }

    virtual scalar_type neg_ll(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns the saddlepoint negative log-likelihood
        scalar_type K_val = K(tvec, parameter_vector);
        vector_type K1_val = K1(tvec, parameter_vector);
        matrix_type K2_val = K2(tvec, parameter_vector);
        //---------------------
        scalar_type res = -tvec.transpose()*K1_val;
        res += K_val;
        res -= (0.5 * (atomic::logdet(matrix<scalar_type>(K2_val)) + tvec.size()*log(2*M_PI)));
        return -res;
    }
//-----------------------------------------------------------
    virtual scalar_type func_T(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns the first higher-order correction term in the saddlepoint expansion

        // Specify Q = K2^{-1} by computing an LDLT decomposition for K2
        // and using it to compute the A and d arguments to the _factored form of the operator methods
        // To save repeated inversion, this is performed once rather than delegated to the _factored_inverse methods
        // Note that with K2 = P^T L D L^T P, we have Q = K2^{-1} = P^{-1} L^{-1}^T D^{-1} L^{-1} P^{-1}^T
        // and since P is a permutation matrix P^T = P^{-1} and this reduces to
        // Q = (L^{-1} P)^T D^{-1} L^{-1} P, i.e., A=(L^{-1}P)^T = P^T (L^T)^{-1}
                                            // auto k2ldlt = K2(tvec, parameter_vector).ldlt();
                                            //
                                            // //matrix_type A = (k2ldlt.matrixL().inverse()*k2ldlt.transpositionsP()).transpose(); // this does not seem to work
                                            // auto dim = tvec.size();
                                            // matrix_type A = matrix_type::Identity(dim, dim) * k2ldlt.transpositionsP().transpose();
                                            // k2ldlt.matrixU().template solveInPlace<Eigen::OnTheRight>(A); // matrixU() provides U=L^T and solveInPlace changes A into A U^{-1}

                                            
        auto inverse_k2ldlt = K2(tvec, parameter_vector).inverse().ldlt();
        matrix_type inverse_k2ldlt_L = inverse_k2ldlt.matrixL();
        matrix_type A = inverse_k2ldlt.transpositionsP().transpose() * inverse_k2ldlt_L;

        vector_type d = inverse_k2ldlt.vectorD();

        scalar_type K4operatorAABB_val = K4operatorAABB_factored(tvec, A, d, A, d, parameter_vector);
        scalar_type K3K3operatorAABBCC_val = K3K3operatorAABBCC_factored(tvec, A, d, A, d, A, d, parameter_vector);
        scalar_type K3K3operatorABCABC_val = K3K3operatorABCABC_factored(tvec, A, d, A, d, A, d, parameter_vector);

//        // Previous version:
//        matrix_type mat_Q = (K2(tvec, parameter_vector)).inverse();
//        scalar_type K4operatorAABB_val = K4operatorAABB(tvec, mat_Q, mat_Q, parameter_vector);
//        scalar_type K3K3operatorAABBCC_val = K3K3operatorAABBCC(tvec, mat_Q, mat_Q, mat_Q, parameter_vector);
//        scalar_type K3K3operatorABCABC_val = K3K3operatorABCABC(tvec, mat_Q, mat_Q, mat_Q, parameter_vector);

        return K4operatorAABB_val/8 - K3K3operatorAABBCC_val/8 - K3K3operatorABCABC_val/12;
    }

//-----------------------------------------------------------
    // "Operator" forms of K2 that does not directly compute K2 but uses it in a quadratic form.
    // Default implementations of these methods call the method K2 above.
    // Overrides of these methods may be implemented more efficiently for specific distributions
    // For instance, if K2 is known to be diagonal, the quadratic form can be computed more efficiently
    // without needing to compute the many 0 entries of K2
    virtual scalar_type K2operator(const vector_type& tvec, const vector_type& x, const vector_type& y,
                                   const vector_type& parameter_vector) const {
        // Returns x^T K2 y, the quadratic form as a function of two (column) vectors x, y
        // x and y must be vectors of the same size as tvec
        return x.transpose() * K2(tvec, parameter_vector) * y;
    }
    virtual matrix_type K2operatorAK2AT(const vector_type& tvec, const matrix_type& A, const vector_type& parameter_vector) const {
        // tvec is a vector of size d = tvec.size()
        // A is a k-by-d matrix for some k
        // Returns the matrix A K2 A^T, which is k-by-k symmetric (and positive semi-definite if K2 is)
        return A * K2(tvec, parameter_vector) * A.transpose();
    }

    // "Operator" forms of K3 and K4, considered as tensors that act on three or four supplied vectors
    // These are the analogues of the quadratic form v1^T K2 v2, as in K2operator
    virtual scalar_type K3operator(const vector_type& tvec, const vector_type& v1, const vector_type& v2, const vector_type& v3,
                                   const vector_type& parameter_vector) const = 0;
    virtual scalar_type K4operator(const vector_type& tvec, const vector_type& v1, const vector_type& v2,
                                   const vector_type& v3, const vector_type& v4, const vector_type& parameter_vector) const = 0;
    // Returns the value of the tensors K3 or K4 applied to the vectors v1,...,v3 or v1,...,v4
    // i.e., the sums K4(i1,i2,i3,i4)*v1[i1]*v2[i2]*v3[i3]*v4[i4]
    // where K4(i1,i2,i3,i4) means the fourth-order mixed partial derivative of K
    // and the sum is over all choices of indices i1,...,i4
    // Note that this is equivalent to the fourth-order mixed partial derivative, with respect to scalars x1, x2, x3, x4,
    // of the function K(tvec + x1*v1 + x2*v2 + x3*v3 + x4*v4, parameter_vector), evaluated at x1=x2=x3=x4=0
    // Similarly for K3operator

//-----------------------------------------------------------
    // Operator forms of K3 and K4 involving matrices instead of vectors
    // Currently, the matrices Q are required to be symmetric
    // (In practice they will be symmetric positive semi-definite as well)
    virtual scalar_type K4operatorAABB(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2,
                                       const vector_type& parameter_vector) const;
        // Returns the sum of K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4) over all indices i1,...,i4
        // where K4(i,j,k,l) denotes the mixed partial derivative of K with respect to the i,j,k,l indices of tvec
    virtual scalar_type K3K3operatorAABBCC(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2,
                                           const matrix_type& Q3, const vector_type& parameter_vector) const;
        // Returns the sum of K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6) over all indices i1,...,i6
        // where K3(i,j,k) denotes the mixed partial derivative of K with respect to the i,j,k indices of tvec
        // Note that Q2 plays a distinguished role compared to Q1 and Q3
    virtual scalar_type K3K3operatorABCABC(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2,
                                           const matrix_type& Q3, const vector_type& parameter_vector) const;
        // Returns the sum of K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6) over all indices i1,...,i6
        // where K3(i,j,k) denotes the mixed partial derivative of K with respect to the i,j,k indices of tvec

//-----------------------------------------------------------
    // Alternative forms of the above operators where the matrices Qk (k=1,2 or k=1,2,3) are specified in factored form

    // Forms where the d-by-d symmetric matrix Qk is specified as Qk = Ak Dk Ak^T, where
    // Ak is d-by-r
    // Dk is r-by-r diagonal with diagonal entries given by the vector dk
    // This factorisation may arise as, for instance,
        // an LDLT factorisation, with Ak square and lower triangular (or permuted lower triangular)
        // an eigenvalue/eigenvector decomposition, with Ak the square orthogonal matrix whose columns are (right) eigenvectors of Qk
            // and dk the vector of eigenvalues of Qk
    // or in any other way.
    // In particular, d and r may be different, i.e., Ak need not be square (and the different Ak's may have different r values)
    // Note that this factorisation implies Qk(i,j) = sum_{m=1,...,r} dk(m)*Ak(i,m)*Ak(j,m)
    virtual scalar_type K4operatorAABB_factored(const vector_type& tvec, const matrix_type& A1, const vector_type& d1,
                                                const matrix_type& A2, const vector_type& d2, const vector_type& parameter_vector) const {
        // Key identity:
        // sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*Q1(i1,i2)*Q2(i3,i4)
        // = sum_{m1=1,...,r1} sum_{m2=1,...,r2} d1(m1)*d2(m2)*( sum_{i1,...i4=1,...,d} K4(i1,i2,i3,i4)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)*A2(i4,m2) )
        // Complexity: r1*r2 calls to K4operator
        typedef decltype(d1.size()) index_type;
        index_type r1 = d1.size();
        index_type r2 = d2.size();

        scalar_type res = 0;
        for(index_type m1 = 0; m1 < r1; ++m1) {
            for (index_type m2 = 0; m2 < r2; ++m2) {
                res += d1(m1) * d2(m2) * K4operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), A2.col(m2), parameter_vector);
            }
        }
        return res;
    }
    virtual scalar_type K3K3operatorAABBCC_factored(const vector_type& tvec,
                                                    const matrix_type& A1, const vector_type& d1,
                                                    const matrix_type& A2, const vector_type& d2,
                                                    const matrix_type& A3, const vector_type& d3,
                                                    const vector_type& parameter_vector) const {
        // Key identity:
        // sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i2)*Q2(i3,i4)*Q3(i5,i6)
        // = sum_{m2=1,...,r2} d2(m2)*( sum_{m1=1,...,r1} d1(m1)*(sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A1(i2,m1)*A2(i3,m2)) )
        //                           *( sum_{m3=1,...,r3} d3(m3)*(sum_{i4,...i6=1,...,d} K3(i4,i5,i6)*A2(i4,m2)*A3(i5,m3)*A3(i6,m3)) )
        // Complexity: r2*(r1+r3) calls to K3operator
        typedef decltype(d1.size()) index_type;
        index_type r1 = d1.size();
        index_type r2 = d2.size();
        index_type r3 = d3.size();

        scalar_type res = 0;
        for(index_type m2 = 0; m2 < r2; ++r2) {
            scalar_type factor1 = 0;
            for(index_type m1 = 0; m1 < r1; ++m1) {
                factor1 += d1(m1) * K3operator(tvec, A1.col(m1), A1.col(m1), A2.col(m2), parameter_vector);
            }
            scalar_type factor2 = 0;
            for(index_type m3 = 0; m3 < r3; ++m3) {
                factor2 += d3(m3) * K3operator(tvec, A2.col(m2), A3.col(m3), A3.col(m3), parameter_vector);
            }
            res += d2(m2)*factor1*factor2;
        }
        return res;
    }
    virtual scalar_type K3K3operatorABCABC_factored(const vector_type& tvec,
                                                    const matrix_type& A1, const vector_type& d1,
                                                    const matrix_type& A2, const vector_type& d2,
                                                    const matrix_type& A3, const vector_type& d3,
                                                    const vector_type& parameter_vector) const {
        // Key identity:
        // sum_{i1,...i6=1,...,d} K3(i1,i2,i3)*K3(i4,i5,i6)*Q1(i1,i4)*Q2(i2,i5)*Q3(i3,i6)
        // = sum_{m1=1,...,r1} sum_{m2=1,...,r2} sum_{m3=1,...,r3} d1(m1)*d2(m2)*d3(m3)*
        //          ( sum_{i1,...i3=1,...,d} K3(i1,i2,i3)*A1(i1,m1)*A2(i2,m2)*A3(i3,m3) )^2
        // since the sum involving i4,i5,i6 is identical to the sum involving i1,i2,i3 except for the labelling of the variables of summation
        // Complexity: r1*r2*r3 calls to K3operator
        typedef decltype(d1.size()) index_type;
        index_type r1 = d1.size();
        index_type r2 = d2.size();
        index_type r3 = d3.size();

        scalar_type res = 0;
        for(index_type m1 = 0; m1 < r1; ++m1) {
            for(index_type m2 = 0; m2 < r2; ++m2) {
                for(index_type m3 = 0; m3 < r3; ++m3) {
                    auto square = [](scalar_type x) {return x*x;};
                    res += d1(m1)*d2(m2)*d3(m3)*square(K3operator(tvec, A1.col(m1), A2.col(m2), A3.col(m3), parameter_vector));
                }
            }
        }
        return res;
    }
//-----------------------------------------------------------
    virtual vector_type ineq_constraint(const vector_type& tvec, const vector_type& parameter_vector) const {
        vector_type res;
        return res;
    }

//-----------------------------------------------------------
    // Scratchpad class and accessor for multiple evaluations
    // In many cases, several CGF functions will be called with the same values tvec and theta
    // The CGF_base class is designed to promote modularity, with many functions implementing
        // various specific and self-contained mathematical building blocks.
    // However, this modularity means that shared sub-expressions may be calculated multiple times,
        // one for each function call.
    // The Scratchpad abstract base class provides an interface for saving a particular tvec, theta pair
        // and evaluating multiple CGF functions with the same function arguments
    // Specific implementations deriving from Scratchpad should save repeated sub-expressions for later reuse
        // as appropriate in particular implementations
    // The related Scratchpad_Q abstract base class saves in addition a square symmetric matrix Q
        // which is supplied as all of the Q1, Q2, Q3 arguments in the operator-style functions
//-----------------------------------------------------------
    class Scratchpad;
    class Scratchpad_Q;

    // Each CGF_base object provides creator methods returning a newly constructed Scratchpad object by pointer.
    // To keep the Scratchpad feature optional, this method is not pure virtual, and the default returns a basic Scratchpad implementation.
    virtual Scratchpad* scratchpad(const vector_type& tvec, const vector_type& parameter_vector) const;
    // The objects returned by the default implementation requires the CGF_base object that created it to remain accessible
    // throughout the duration of the returned Scratchpad object's lifetime.
    // Other implementations might not impose this restriction.
    // It is the user's responsibility to manage the pointer returned by this method.
};

//-----------------------------------------------------------
// Implementation of default K4operatorAABB, K3operatorAABBCC, K3operatorABCABC
//-----------------------------------------------------------
// All three implementations compute LDLT decompositions of Q1,Q2 or Q1,Q2,Q3 and pass to the corresponding _factored methods
    // Namely, Qk = Pk^T Lk Dk Lk^T Pk where
    // Pk is a permutation matrix (accessed implicitly via .transpositionsP() method)
    // Lk is lower triangular
    // Dk is diagonal
// Note: *this code assumes that matrix_type supports the Eigen library's format for LDLT*

// Note that if K3 and K4 are "diagonal" (i.e., if K3(i,j,k) is zero unless i==j==k, and similarly for K4)
// then this implementation is inefficient.
// For this case, a more efficient implementation can simply extract the diagonals of Q1,Q2,Q3

template <class scalar_type, class vector_type, class matrix_type>
scalar_type CGF_base<scalar_type, vector_type, matrix_type>::
    K4operatorAABB(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2, const vector_type& parameter_vector) const {
        auto ldlt1 = Q1.ldlt();
        auto ldlt2 = Q2.ldlt();
        auto extract_A = [](const decltype(ldlt1)& ldlt) {
            matrix_type res = ldlt.matrixL();
            res = ldlt.transpositionsP().transpose() * res;
            return res;
            };
        return CGF_base::K4operatorAABB_factored(tvec,
                                                 extract_A(ldlt1), ldlt1.vectorD(),
                                                 extract_A(ldlt2), ldlt2.vectorD(), parameter_vector);
}
template <class scalar_type, class vector_type, class matrix_type>
scalar_type CGF_base<scalar_type, vector_type, matrix_type>::
    K3K3operatorAABBCC(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                       const vector_type& parameter_vector) const {
        auto ldlt1 = Q1.ldlt();
        auto ldlt2 = Q2.ldlt();
        auto ldlt3 = Q3.ldlt();
        auto extract_A = [](const decltype(ldlt1)& ldlt) {
            matrix_type res = ldlt.matrixL();
            res = ldlt.transpositionsP().transpose() * res;
            return res;
            };
        return CGF_base::K3K3operatorAABBCC_factored(tvec,
                                                     extract_A(ldlt1), ldlt1.vectorD(),
                                                     extract_A(ldlt2), ldlt2.vectorD(),
                                                     extract_A(ldlt3), ldlt3.vectorD(), parameter_vector);
}
template <class scalar_type, class vector_type, class matrix_type>
scalar_type CGF_base<scalar_type, vector_type, matrix_type>::
    K3K3operatorABCABC(const vector_type& tvec, const matrix_type& Q1, const matrix_type& Q2, const matrix_type& Q3,
                       const vector_type& parameter_vector) const {
        auto ldlt1 = Q1.ldlt();
        auto ldlt2 = Q2.ldlt();
        auto ldlt3 = Q3.ldlt();
        auto extract_A = [](const decltype(ldlt1)& ldlt) {
            matrix_type res = ldlt.matrixL();
            res = ldlt.transpositionsP().transpose() * res;
            return res;
            };
        return CGF_base::K3K3operatorABCABC_factored(tvec,
                                                     extract_A(ldlt1), ldlt1.vectorD(),
                                                     extract_A(ldlt2), ldlt2.vectorD(),
                                                     extract_A(ldlt3), ldlt3.vectorD(), parameter_vector);
}

} // namespace CGFs_via_virtual_functions
} // namespace saddlepoint

#include "scratchpad.h"

namespace saddlepoint {
namespace CGFs_with_AD {

#include "needs_full_list_of_CGF_methods.h"
//-----------------------------------------------------------
// CGF_with_AD: the class we will use, which is a join class implementing the CGF_base interface with two sets of template parameters
class CGF_with_AD : public CGFs_via_virtual_functions::CGF_base<double, vec, mat>,
                    public CGFs_via_virtual_functions::CGF_base<a_scalar, a_vector, a_matrix> {
public:
    // To allow overload resolution to operate within CGF_with_AD,
        // explicit using-declarations bring the normal and AD versions into the same scope
    using CGF_base<double, vec, mat>::K;
    using CGF_base<double, vec, mat>::K1;
    using CGF_base<double, vec, mat>::K2;
    using CGF_base<double, vec, mat>::tilting_exponent;
    using CGF_base<double, vec, mat>::neg_ll;
    using CGF_base<double, vec, mat>::func_T;
    using CGF_base<double, vec, mat>::K2operator;
    using CGF_base<double, vec, mat>::K2operatorAK2AT;
    using CGF_base<double, vec, mat>::K3operator;
    using CGF_base<double, vec, mat>::K4operator;
    using CGF_base<double, vec, mat>::K4operatorAABB;
    using CGF_base<double, vec, mat>::K3K3operatorAABBCC;
    using CGF_base<double, vec, mat>::K3K3operatorABCABC;
    using CGF_base<double, vec, mat>::K4operatorAABB_factored;
    using CGF_base<double, vec, mat>::K3K3operatorAABBCC_factored;
    using CGF_base<double, vec, mat>::K3K3operatorABCABC_factored;
    using CGF_base<double, vec, mat>::scratchpad;
    using CGF_base<double, vec, mat>::ineq_constraint;

    using CGF_base<a_scalar, a_vector, a_matrix>::K;
    using CGF_base<a_scalar, a_vector, a_matrix>::K1;
    using CGF_base<a_scalar, a_vector, a_matrix>::K2;
    using CGF_base<a_scalar, a_vector, a_matrix>::tilting_exponent;
    using CGF_base<a_scalar, a_vector, a_matrix>::neg_ll;
    using CGF_base<a_scalar, a_vector, a_matrix>::func_T;
    using CGF_base<a_scalar, a_vector, a_matrix>::K2operator;
    using CGF_base<a_scalar, a_vector, a_matrix>::K2operatorAK2AT;
    using CGF_base<a_scalar, a_vector, a_matrix>::K3operator;
    using CGF_base<a_scalar, a_vector, a_matrix>::K4operator;
    using CGF_base<a_scalar, a_vector, a_matrix>::K4operatorAABB;
    using CGF_base<a_scalar, a_vector, a_matrix>::K3K3operatorAABBCC;
    using CGF_base<a_scalar, a_vector, a_matrix>::K3K3operatorABCABC;
    using CGF_base<a_scalar, a_vector, a_matrix>::K4operatorAABB_factored;
    using CGF_base<a_scalar, a_vector, a_matrix>::K3K3operatorAABBCC_factored;
    using CGF_base<a_scalar, a_vector, a_matrix>::K3K3operatorABCABC_factored;
    using CGF_base<a_scalar, a_vector, a_matrix>::scratchpad;
    using CGF_base<a_scalar, a_vector, a_matrix>::ineq_constraint;
};
//-----------------------------------------------------------

#include "needs_full_list_of_CGF_methods.h"
template <class templateCGF>
class CGF_with_AD_from_template : public CGF_with_AD, protected templateCGF {
public:
    using templateCGF::templateCGF;

    vec ineq_constraint(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::ineq_constraint(tvec, parameter_vector);}
    a_vector ineq_constraint(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::ineq_constraint(tvec, parameter_vector);}



    double K(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K(tvec, parameter_vector);}
    vec K1(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K1(tvec, parameter_vector);}
    mat K2(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K2(tvec, parameter_vector);}
    double tilting_exponent(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::tilting_exponent(tvec, parameter_vector);}
    double neg_ll(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::neg_ll(tvec, parameter_vector);}
    double func_T(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::func_T(tvec, parameter_vector);}
    double K2operator(const vec& tvec, const vec& x, const vec& y, const vec& parameter_vector) const override {return templateCGF::K2operator(tvec, x, y, parameter_vector);}
    mat K2operatorAK2AT(const vec& tvec, const mat& A, const vec& parameter_vector) const override {return templateCGF::K2operatorAK2AT(tvec, A, parameter_vector);}
    double K3operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& parameter_vector) const override {return templateCGF::K3operator(tvec, v1, v2, v3, parameter_vector);}
    double K4operator(const vec& tvec, const vec& v1, const vec& v2, const vec& v3, const vec& v4, const vec& parameter_vector) const override {return templateCGF::K4operator(tvec, v1, v2, v3, v4, parameter_vector);}
    double K4operatorAABB(const vec& tvec, const mat& Q1, const mat& Q2, const vec& parameter_vector) const override {return templateCGF::K4operatorAABB(tvec, Q1, Q2, parameter_vector);}
    double K3K3operatorAABBCC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {return templateCGF::K3K3operatorAABBCC(tvec, Q1, Q2, Q3, parameter_vector);}
    double K3K3operatorABCABC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {return templateCGF::K3K3operatorABCABC(tvec, Q1, Q2, Q3, parameter_vector);}
    double K4operatorAABB_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const vec& parameter_vector) const override {return templateCGF::K4operatorAABB_factored(tvec, A1, d1, A2, d2, parameter_vector);}
    double K3K3operatorAABBCC_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const mat& A3, const vec& d3, const vec& parameter_vector) const override {return templateCGF::K3K3operatorAABBCC_factored(tvec, A1, d1, A2, d2, A3, d3, parameter_vector);}
    double K3K3operatorABCABC_factored(const vec& tvec, const mat& A1, const vec& d1, const mat& A2, const vec& d2, const mat& A3, const vec& d3, const vec& parameter_vector) const override {return templateCGF::K3K3operatorABCABC_factored(tvec, A1, d1, A2, d2, A3, d3, parameter_vector);}
    // currently scratchpad() is not supported by template versions so use default from CGF_base

    a_scalar K(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K(tvec, parameter_vector);}
    a_vector K1(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K1(tvec, parameter_vector);}
    a_matrix K2(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K2(tvec, parameter_vector);}
    a_scalar tilting_exponent(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::tilting_exponent(tvec, parameter_vector);}
    a_scalar neg_ll(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::neg_ll(tvec, parameter_vector);}
    a_scalar func_T(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::func_T(tvec, parameter_vector);}
    a_scalar K2operator(const a_vector& tvec, const a_vector& x, const a_vector& y, const a_vector& parameter_vector) const override {return templateCGF::K2operator(tvec, x, y, parameter_vector);}
    a_matrix K2operatorAK2AT(const a_vector& tvec, const a_matrix& A, const a_vector& parameter_vector) const override {return templateCGF::K2operatorAK2AT(tvec, A, parameter_vector);}
    a_scalar K3operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& parameter_vector) const override {return templateCGF::K3operator(tvec, v1, v2, v3, parameter_vector);}
    a_scalar K4operator(const a_vector& tvec, const a_vector& v1, const a_vector& v2, const a_vector& v3, const a_vector& v4, const a_vector& parameter_vector) const override {return templateCGF::K4operator(tvec, v1, v2, v3, v4, parameter_vector);}
    a_scalar K4operatorAABB(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_vector& parameter_vector) const override {return templateCGF::K4operatorAABB(tvec, Q1, Q2, parameter_vector);}
    a_scalar K3K3operatorAABBCC(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {return templateCGF::K3K3operatorAABBCC(tvec, Q1, Q2, Q3, parameter_vector);}
    a_scalar K3K3operatorABCABC(const a_vector& tvec, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {return templateCGF::K3K3operatorABCABC(tvec, Q1, Q2, Q3, parameter_vector);}
    a_scalar K4operatorAABB_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_vector& parameter_vector) const override {return templateCGF::K4operatorAABB_factored(tvec, A1, d1, A2, d2, parameter_vector);}
    a_scalar K3K3operatorAABBCC_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_matrix& A3, const a_vector& d3, const a_vector& parameter_vector) const override {return templateCGF::K3K3operatorAABBCC_factored(tvec, A1, d1, A2, d2, A3, d3, parameter_vector);}
    a_scalar K3K3operatorABCABC_factored(const a_vector& tvec, const a_matrix& A1, const a_vector& d1, const a_matrix& A2, const a_vector& d2, const a_matrix& A3, const a_vector& d3, const a_vector& parameter_vector) const override {return templateCGF::K3K3operatorABCABC_factored(tvec, A1, d1, A2, d2, A3, d3, parameter_vector);}
};

} // namespace CGFs_with_AD
} // namespace saddlepoint

#endif // BASECGF_H_INCLUDED
