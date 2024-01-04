#ifndef ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED
#define ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED

#include <utility>
#include <cppad/cppad.hpp>
#include <Eigen/Eigen>

// Assumptions:
// This atomic function implements the derivative behaviour pertaining to an implicitly defined function g(v), where
// u = g(v) is the solution of f(u, v) = constant.
// The behaviour of the function f is specified via a FunctionObject.
// It is required that, at the points of interest, the partial derivaive of f with respect to u is non-zero.
// Both u and v may be vectors, say of dimensions n and m, in which case f must have a vector result of dimension n,
// and the n-by-n matrix of partial derivatives of f with respect to u is required to be non-singular.

// Calculus underpinnings:
// Given vectors u0, v0 such that df/du(u0,v0) is non-singular, there is a (mathematical) function hat_u(v)
// defined in a neighbourhood of v0 such that f(hat_u(v),v) = constant and hat_u(v0)=u0.
// The constant is precisely the value f(u0,v0).
// (Note that the (mathematical) function hat_u depends on the choice of points u0,v0, in general,
// although within a local neighbourhood the (mathematical) function depends only on the value of f(u0,v0).)
// Differentiating leads to 0 = df/du(hat_u(v),v) * dhat_u/dv(v) + df/dv(hat_u(v),v), so that
// dhat_u/dv = - (df/du)^{-1} * df/dv

class atomic_implicit_function_common_base : public CppAD::atomic_three<double> {
    // Helper class to implement for_type, which does not vary across the template parameter of atomic_implicit_function
public:
    using atomic_three<double>::atomic_three;
private:
    bool for_type(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        CppAD::vector<CppAD::ad_type_enum>&        type_y      ) override
        {
          size_t n = type_y.size(); // size of u
          bool ok = (type_x.size() > 2*n); // x.size() should equal 2*n+m where m==v.size()
          if (!ok) {return ok;}
          size_t m = type_x.size() - 2*n;

          // Convenience mappings to Eigen types
          using Eigen::Map; using Eigen::Matrix; using Eigen::Dynamic;
          typedef Matrix<CppAD::ad_type_enum, Dynamic, 1> Matrix_type;
          Map<const Matrix_type> type_x_mapped(type_x.data(), type_x.size(), 1);
          Map<Matrix_type> type_y_mapped(type_y.data(), type_y.size(), 1);

          // Fill type_y with max of type_x entries that correspond to v
          // i.e., the max of the size-m segment of type_x starting just after the first n entries
          type_y_mapped.setConstant(type_x_mapped.segment(n,m).maxCoeff());

          return ok;
        }
};

template <class FunctionObject>
class atomic_implicit_function : public atomic_implicit_function_common_base, protected FunctionObject {
    // Usage: the atomic function is not intended to be called directly by the user.
    // In order to invoke
        // atomic_implicit_function aif(a_x,a_y);
    // the caller must ensure that, in terms of the mathematical vectors u,v described above,
        // a_y.size() == u.size();
        // a_x.size() == 2*u.size() + v.size()
        // a_x is the concatenation of the vectors u, v, and f(u,v), in that order
    // It is irrelevant how the AD vector u has been created as a function of v,
    // and indeed u may have been coerced to appear constant.
    // No matter what the dependence of u on underlying variables, atomic_implicit_function will
    // produce a result that appears to depend correctly on v only.
    // However, it is essential that the third block of a_x should be
    // evaluated with the same vectors u,v passed to atomic_implicit_function.
    // See the template function convert_to_implicit below.

    // Implementation:
    // FunctionObject must provide methods compatible with:
        // Eigen::VectorXd y = FunctionObject::dfdu_solve(const Eigen::VectorXd& u,
        //                                                const Eigen::VectorXd& v,
        //                                                const Eigen::VectorXd& w);
        // where y and w are (column) vectors of the same size as u.
        // dfdu_solve should return (df/du)^{-1} * w, where df/du is evaluated at (u,v)
    // and
        // Eigen::VectorXd yT = FunctionObject::dfdu_solve_row(const Eigen::VectorXd& u,
        //                                                     const Eigen::VectorXd& v,
        //                                                     const Eigen::VectorXd& wT);
        // where the column vectors yT and wT are the transposes of row vectors of the same size as u.
        // dfdu_solve_row should return the column vector (w * (df/du)^{-1})^T, where df/du is evaluated at (u,v)
        // or in other words yT = ((df/du)^T)^{-1} * wT.
        // Mathematically, this is simply right-multiplication by (df/du)^{-1},
        // except that all vectors are changed into column vectors in order to keep to a single vector type.
public:
    template <class... FOArgTypes>
    atomic_implicit_function(const std::string& name, FOArgTypes&&... args)
        : atomic_implicit_function_common_base(name), FunctionObject(std::forward<FOArgTypes>(args)...) {}
private:
    bool forward(
        const CppAD::vector<double>&              parameter_x  ,
        const CppAD::vector<CppAD::ad_type_enum>& type_x       ,
        size_t                                    need_y       ,
        size_t                                    order_low    ,
        size_t                                    order_up     ,
        const CppAD::vector<double>&              taylor_x     ,
        CppAD::vector<double>&                    taylor_y     ) override
    {
        // General preconditions:
        // Call the notional input vector x and the notional output vector y.
        // With q_1 = order_up+1, the entries of x appear in entries taylor_x[i*q1 + 0], 0 <= i < x.size()
        // When order_up >= 1, the k^th-order Taylor coefficients of the i^th entry of x
        // appear in entries taylor_x[i*q1 + k], 0 <= i < x.size(), 0 <= k < q_1
        // where the Taylor series should be understood as being with respect to an unspecified underlying scalar variable, x=x(t) for real t.

        // General postconditions:
        // Every entry taylor_y[j*q1 + k], 0 <= j < y.size(), 0 <= k < q_1, contains the k^th-order Taylor coefficient of the j^th entry of y

        // Here we only implement up to order 1, and 0 <= order_low <= order_up <= 1

        // Preconditions:
        // x is the concatenation of n entries of u, m entries of v, and n entries of f(u,v)
        // df/du(u,v) is a non-singular square matrix
        // 0 <= order_low <= order_up <= 1

        // Postconditions:
        // The order-0 entries of taylor_y are determined consistent with the assumption y = hat_u = u
        // If order_up >= 1, the order-1 entries of taylor_y are determined so that,
        // in combination with the assumption that the third block of x is f(u,v),
        // the overall result is that dy/du = 0
        // and dy/dv = (df/du)^{-1}*df/dv, where the partial derivatives are evaluated at (u,v)

        // In more detail for order 1:
        // Call the input values u0, v0, c0 = f(u0,v0).
        // The AD assumption for forward mode is that there are unspecified functions u(t), v(t), c(t)
        // such that u(0)==u0, v(0)==v0, c(0)==c0, and such that the derivatives u'(0), v'(0), c'(0)
        // are stored in taylor_x.
        // Moreover, because we assume that the input c has been calculated as c = f(u,v) in the sense of AD values,
        // we can assume that that c(t) = f(u(t),v(t)).
        // The function we are looking for is hat_u(t) such that f(hat_u(t), v(t)) == c0 for all t,
        // and in forward mode of order 1 we must store the value hat_u'(0) into taylor_y.
        // Note that hat_u(0)=u0.
        // Differentiating w.r.t. t and setting t=0, we find
            // c'(0) == df/du * u'(0) + df/dv * v'(0)
            // 0 == df/du * hat_u'(0) + df/dv * v'(0)
        // and therefore
            // hat_u'(0) = u'(0) - (df/du)^{-1} c'(0)
        // (Note that this identity is valid no matter how u(t) depends on t.)

        bool ok = (order_up <= 1);
        // Order 0 and 1 implemented only
        if (!ok) {return ok;}

        // Infer and check sizes
        size_t q1 = order_up + 1;
        size_t n = taylor_y.size() / q1; // size of u
        ok &= (type_x.size() > 2*n); // x.size() should equal 2*n+m where m==v.size()
        if (!ok) {return ok;}
        size_t m = type_x.size() - 2*n;

        // Spacing (stride) q1 for input vector taylor_x and output vector taylor_y
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_u(taylor_x.data(), n, Eigen::InnerStride<>(q1));
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_y(taylor_y.data(), n, Eigen::InnerStride<>(q1));

        // Order zero forward mode
        if( order_low <= 0 ){
            // Function value is the part of the input vector supplied by the user
            mapped_y = mapped_u;
        }
        if( order_up <= 0 ) {
            return ok;
        }

        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_v(taylor_x.data() + n*q1, m, Eigen::InnerStride<>(q1));
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_u_prime(taylor_x.data() + 1, n, Eigen::InnerStride<>(q1));
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_c_prime(taylor_x.data() + (n+m)*q1 + 1, n, Eigen::InnerStride<>(q1));
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>>
            mapped_y_prime(taylor_y.data() + 1, n, Eigen::InnerStride<>(q1));

        // Order one forward mode.
        if( order_low <= 1 ) {
            mapped_y_prime = mapped_u_prime - FunctionObject::dfdu_solve(mapped_u, mapped_v, mapped_c_prime);
        }
        // We only implement forward mode up to order 1
        return ok;
    }

    bool reverse(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        size_t                              order_up    ,
        const CppAD::vector<double>&               taylor_x    ,
        const CppAD::vector<double>&               taylor_y    ,
        CppAD::vector<double>&                     partial_x   ,
        const CppAD::vector<double>&               partial_y   ) override
    {
        // General preconditions:
        // Call the notional input vector x.
        // With q_1 = order_up+1, the entries of x appear in entries taylor_x[i*q1 + 0], 0 <= i < x.size()
        // When order_up >= 1, the k^th-order Taylor coefficients of the i^th entry of x
        // appear in entries taylor_x[i*q1 + k], 0 <= i < x.size(), 0 <= k < q_1
        // [where the Taylor series should be understood as being with respect to an unspecified underlying scalar variable]
        // taylor_y has the same structure relative to the notional output vector y,
        // and in particular the output values taylor_y[j*q1 + 0], 0 <= j < y.size() are already computed
        // partial_y has the same arrangement as taylor_y
        // and each entry partial_y[k] should be interpreted as the partial derivative with respect to taylor_y[k]
        // of a notional scalar function G(taylor_y) applied to the output vector taylor_y

        // General postcondition:
        // partial_x[m] represents the partial derivative with respect to taylor_x[m]
        // of the composite scalar function H(F(taylor_x)),
        // where F is the function that expresses taylor_y in terms of taylor_x
        // as determined according to the underlying function h represented by the atomic function being implemented

        // Here we only implement order_up == 0, q_1 == 1,
        // so that taylor_x and taylor_y can be identified with x and y
        // Thus partial_y[k] represents dG/dy_k and we are to compute dH/dx_m, which are related by
        // dH/dx_m = sum_k dG/dy_k * dy_k/dx_m

        // Preconditions:
        // x is the concatenation of n entries of u, m entries of v, and n entries of f(u,v)
        // df/du(u,v) is a non-singular square matrix
        // order_up == 0

        // Postconditions:
        // The entries of partial_x are determined so that,
        // in combination with the assumption that the third block of x is f(u, v),
        // the overall result is that dy/du = 0 and dy/dv = - (df/du)^{-1} * df/dv

        // In more detail:
        // Call the input values u0, v0, c0 = f(u0,v0).
        // The AD assumption for reverse mode is that there is an unspecified scalar-valued function G(y)
        // such that the entries of partial_y are the partial derivatives dG/dy_i, i=1,...,n, evaluated at y = u_hat.
        // Since y = y(u,v,c) is considered to be a function of u, v, c, we can define a function H(u,v,c) = G(y(u,v,c))
        // and we should store partial derivative values for H into partial_x
        // based on knowledge of the function y(u,v,c) we are trying to represent.
        // Now add the additional assumption that c = f(u1,v1), where for the moment u1 and v1 are formally distinct
        // copies of u and v.
        // We know that later reverse-mode operations for the mapping (u, v, u1, v1) -> (u, v, f(u1, v1))
        // will produce a function Q(u,v,u1,v1) = H(G(y(u,v,f(u1,v1))))
        // and will map the partial derivatives of H into
            // dQ/du = dH/du and dQ/dv = dH/dv,
            // dQ/du1 = dH/dc * df/du1 and dQ/dv1 = dH/dc * df/dv1.
        // We wish to match this with the function R(v) = G(hat_u(v)) of v alone,
        // under the additional assumption u1 = u, v1 = v.
        // So we require
            // dQ/du + dQ/du1 = dR/du = 0
            // dQ/dv + dQ/dv1 = dR/dv = dG/dy dhat_u/dv = - dG/dy * (df/du)^{-1} df/dv
        // which reduces to
            // dH/du + dH/dc * df/du = 0
            // dH/dv + dH/dc * df/dv = - dG/dy * (df/du)^{-1} df/dv
        // Eliminating dH/dc = - dH/du * (df/du)^{-1} from the second equation leads to
            // dH/dv = ( dH/du - dG/dy ) * (df/du)^{-1} * df/dv
        // This is underdetermined, but we can choose a solution that does not require knowledge of df/dv by taking
            // dH/du = dG/dy
            // dH/dv = 0
            // dH/dc = - dG/dy * (df/du)^{-1}
        // Note the duality with the forward case discussed above: if we are given
        // u(t), v(t) and we set c(t) = f(u(t),v(t)) and y(t) = hat_u(v(t)),
        // and H(t) = H(u(t), v(t), c(t)) and G(t) = G(y(t)), then these partial derivatives yield
            // H'(0) = dH/du u'(0) + dH/dv v'(0) + dH/dc c'(0)
            // = dG/dy u'(0) + 0 - dG/dy * (df/du)^{-1} (df/du u'(0) + df/dv v'(0))
            // = dG/dy ( u'(0) - (df/du)^{-1} (df/du) u'(0) - (df/du)^{-1} df/dv v'(0))
            // = - dG/dy (df/du)^{-1} df/dv v'(0)
        // which matches with
            // G'(0) = dG/dy y'(0) = dG/dy dhat_u/dv = - dG/dy (df/du)^{-1 df/dv
        // as required.

        // First-order reverse mode only.
        // Note that the "order" for reverse mode is order_up+1.
        size_t q1 = order_up + 1;
        bool ok = (q1 == 1);
        if (!ok) {return ok;}
        // The rest of this function makes explicit use of the assumption q1 == 1.

        // Infer sizes (no checking because this will have been done by forward() on a previous occasion)
        size_t n = taylor_y.size(); // size of u
        size_t m = type_x.size() - 2*n;

        // dH/du = dG/dy
        Eigen::Map<Eigen::VectorXd> partial_u(partial_x.data(), n);
        partial_u = Eigen::Map<const Eigen::VectorXd>(partial_y.data(), n);
        // dH/dv = 0
        Eigen::Map<Eigen::VectorXd> partial_v(partial_x.data()+n, m);
        partial_v.setConstant(0);

        Eigen::Map<Eigen::VectorXd> partial_c(partial_x.data()+n+m, n);
        Eigen::Map<const Eigen::VectorXd> mapped_u(taylor_y.data(), n);
        Eigen::Map<const Eigen::VectorXd> mapped_v(taylor_x.data()+n, m);
        // Note: the value of u can be read either from taylor_y or from the first part of taylor_x.
        // This makes no difference here.
        // TO DO: consider whether the AD<double> version of this method should use one or the other.

        // dH/dc = - dG/dy * (df/du)^{-1}
        partial_c = - FunctionObject::dfdu_solve_row(mapped_u,
                                                     mapped_v,
                                                     Eigen::Map<const Eigen::VectorXd>(partial_y.data(), n));
        return ok;
    }
}; // End of class atomic_implicit_function

// Helper class, see convert_to_implicit below
template <class dfdu_FO>
struct convert_to_implicit_helper_FO : public dfdu_FO {
public:
    template <class u_type, class v_type, class w_type>
    auto dfdu_solve(const u_type& u, const v_type& v, const w_type& w) {
        return dfdu_FO::operator()(u, v).solve(w);
    }
    template <class u_type, class v_type, class wT_type>
    auto dfdu_solve_row(const u_type& u, const v_type& v, const wT_type& wT) {
        return (dfdu_FO::operator()(u, v)).transpose().solve(wT);
    }
};


template <class f_FO, class dfdu_FO>
CppAD::vector<CppAD::AD<double>>
convert_to_implicit(const CppAD::vector<CppAD::AD<double>>& a_u, const CppAD::vector<CppAD::AD<double>>& a_v)
// f_FO must provide f_FO::operator()(a_u, a_v) returning a CppAD::vector<CppAD::AD<double>>
// dfdu_FO must provide dfdu_FO::operator()(u, v) (non-AD arguments) returning an Eigen::MatrixXd or compatible type
// These operators should compute a function f(u,v) and its partial derivative matrix df/du, respectively
// (with the i,j entry of df/du being the df_i/du_j)
// f_FO and dfdu_FO must be default-constructible
// This function could also be adapted to accept specified functions but this is not currently done.
{
    static atomic_implicit_function<convert_to_implicit_helper_FO<dfdu_FO>> aif;
    static f_FO f;

    typedef Eigen::Map<Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1>> MappedADVector;
    typedef Eigen::Map<const Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1>> ConstMappedADVector;

    size_t n = a_u.size();
    size_t m = a_v.size();
    CppAD::vector<CppAD::AD<double>> a_x(2*n+m);
    MappedADVector mapped_a_x(a_x.data(), a_x.size());
    mapped_a_x.head(n) = ConstMappedADVector(a_u.data(), n);
    mapped_a_x.segment(n, m) = ConstMappedADVector(a_v.data(), m);

    CppAD::vector<CppAD::AD<double>> fuv = f(a_u, a_v);
    if (fuv.size() != n) {
        throw std::range_error("Function object must return vector of same size as first argument in convert_to_implicit()");
    }
    mapped_a_x.tail(n) = ConstMappedADVector(fuv, n);
    CppAD::vector<CppAD::AD<double>> a_y(n);
    aif(a_x, a_y);
    return a_y;
}

#endif // ATOMIC_IMPLICIT_FUNCTION_H_INCLUDED
