#ifndef RANDOMLYSTOPPEDSUMCGF_H
#define RANDOMLYSTOPPEDSUMCGF_H

#include "CGF_Defaults.h"
#include "BaseWrapper.h"
#include "baseCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

template <class CountCGF, class SummandCGF>
class RandomlyStoppedSumCGF : public CGF_Defaults<RandomlyStoppedSumCGF<CountCGF,SummandCGF>>, private BaseWrapper<CountCGF, 0>, BaseWrapper<SummandCGF, 1> {
    // CGF for Y = sum_{i=1}^N X_i, where
    // N is a scalar given by CountCGF,
    // X_i are i.i.d. given by SummandCGF,
    // and N and the X_i's are independent.
private:
    auto count_cgf() {return BaseWrapper<CountCGF, 0>::base_cgf();}
    auto count_cgf() const {return BaseWrapper<CountCGF, 0>::base_cgf();}
    auto summand_cgf() {return BaseWrapper<SummandCGF, 1>::base_cgf();}
    auto summand_cgf() const {return BaseWrapper<SummandCGF, 1>::base_cgf();}
    // CountCGF and SummandCGF may be the same type; the integer argument to BaseWrapper allows them to be distinguished

public:
    template <class CountCGF_Type, class SummandCGF_Type>
    RandomlyStoppedSumCGF(CountCGF_Type&& c_cgf, SummandCGF_Type&& s_cgf)
            : BaseWrapper<CountCGF, 0>(std::forward<CountCGF_Type>(c_cgf)), BaseWrapper<SummandCGF, 1>(std::forward<SummandCGF_Type>(s_cgf)) {}

    template <class t_vector_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const ParamTypes&... params) const {
        // K_Y = K_N( K_X(tvec, params), params)
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));
        return count_cgf()->K(count_tvec, params...);
    }
//---------------
    template <class t_vector_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const ParamTypes&... params) const {
        // K1_Y = K1_N( K_X(tvec, params), params) K1_X(tvec, params)
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));
        return ((count_cgf()->K1(count_tvec, params...)).value() * (summand_cgf()->K1(tvec, params...)).array()).eval();
    }
//---------------
    template <class t_vector_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const ParamTypes&... params) const {
        // K2_Y = (K1_X) (##) (K1_X)^T + (#) K2_X
        // #  : K1_N(K_X(tvec, params), params) // this will be a scalar
        // ## : K2_N(K_X(tvec, params), params) // this is also a scalar
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));
        auto K1_N = count_cgf()->K1(count_tvec, params...).eval(); // this is still of a vector type, but of size 1
        auto K2_N = count_cgf()->K2(count_tvec, params...).eval(); // (1 by 1) matrix
        auto K1_X = summand_cgf()->K1(tvec, params...).eval();
        auto K2_X = summand_cgf()->K2(tvec, params...).eval();
        return (K1_X * K2_N * K1_X.transpose() + K1_N.value() * K2_X).eval();
    }
//---------------
    template <class t_vector_type, class w1_type, class w2_type, class w3_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const w1_type& w1, const w2_type& w2, const w3_type& w3, const ParamTypes&... params) const  {
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        //----------
        //----------
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));
        auto K1_N = count_cgf()->K1(count_tvec, params...).eval();
        auto K1_X = summand_cgf()->K1(tvec, params...).eval();
        //----------
        //----------
        return K1_N.value()*summand_cgf()->K3operator(tvec, w1, w2, w3, params...) +
               count_cgf()->K2operator(count_tvec, w1.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w2, w3, params...)), params...) +
               count_cgf()->K2operator(count_tvec, w2.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w1, w3, params...)), params...) +
               count_cgf()->K2operator(count_tvec, w3.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w1, w2, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w1.transpose()*K1_X, w2.transpose()*K1_X, w3.transpose()*K1_X, params...);
    }
//---------------
    template <class t_vector_type, class w1_type, class w2_type, class w3_type, class w4_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const w1_type& w1, const w2_type& w2, const w3_type& w3, const w4_type& w4, const ParamTypes&... params) const  {
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));

        auto K1_N = count_cgf()->K1(count_tvec, params...).eval();
        auto K1_X = summand_cgf()->K1(tvec, params...).eval();
        //----------
        //----------
        return K1_N.value()*summand_cgf()->K4operator(tvec, w1, w2, w3, w4, params...) +
               count_cgf()->K2operator(count_tvec, w1.transpose()*K1_X, s_to_v1(summand_cgf()->K3operator(tvec, w2, w3, w4, params...)), params...) +
               count_cgf()->K2operator(count_tvec, w2.transpose()*K1_X, s_to_v1(summand_cgf()->K3operator(tvec, w1, w3, w4, params...)), params...) +
               count_cgf()->K2operator(count_tvec, w3.transpose()*K1_X, s_to_v1(summand_cgf()->K3operator(tvec, w1, w2, w4, params...)), params...) +
               count_cgf()->K2operator(count_tvec, w4.transpose()*K1_X, s_to_v1(summand_cgf()->K3operator(tvec, w1, w2, w3, params...)), params...) +
               count_cgf()->K2operator(count_tvec, s_to_v1(summand_cgf()->K2operator(tvec, w1, w2, params...)), s_to_v1(summand_cgf()->K2operator(tvec, w3, w4, params...)), params...) +
               count_cgf()->K2operator(count_tvec, s_to_v1(summand_cgf()->K2operator(tvec, w1, w4, params...)), s_to_v1(summand_cgf()->K2operator(tvec, w2, w3, params...)), params...) +
               count_cgf()->K2operator(count_tvec, s_to_v1(summand_cgf()->K2operator(tvec, w1, w3, params...)), s_to_v1(summand_cgf()->K2operator(tvec, w2, w4, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w1.transpose()*K1_X, w2.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w3, w4, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w1.transpose()*K1_X, w3.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w2, w4, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w1.transpose()*K1_X, w4.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w2, w3, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w2.transpose()*K1_X, w3.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w1, w4, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w2.transpose()*K1_X, w4.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w1, w3, params...)), params...) +
               count_cgf()->K3operator(count_tvec, w3.transpose()*K1_X, w4.transpose()*K1_X, s_to_v1(summand_cgf()->K2operator(tvec, w1, w2, params...)), params...) +
               count_cgf()->K4operator(count_tvec, w1.transpose()*K1_X, w2.transpose()*K1_X, w3.transpose()*K1_X, w4.transpose()*K1_X, params...);
    }
//---------------
    template <class t_vector_type, class x_type, class y_type, class... ParamTypes>
    auto K2operator(const t_vector_type& tvec, const x_type& x, const y_type& y, const ParamTypes&... params) const  {
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));
        auto K1_N = count_cgf()->K1(count_tvec, params...).eval();
        auto K1_X = summand_cgf()->K1(tvec, params...).eval();
        //----------
        //----------
        return K1_N.value()*summand_cgf()->K2operator(tvec, x.transpose(), y, params...) +
               count_cgf()->K2operator(count_tvec, K1_X.transpose()*x, K1_X.transpose()*y, params...);
    }
//---------------
    template <class t_vector_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef decltype(summand_cgf()->K(tvec, params...)) scalar_type;
        auto s_to_v1 = [](auto x){return Eigen::Matrix<scalar_type, 1, 1>{x};}; // scalar to vector of size 1 conversion
        Eigen::Matrix<scalar_type, 1, 1> count_tvec = s_to_v1(summand_cgf()->K(tvec, params...));


        auto summand_ineq_vec = summand_cgf()->ineq_constraint(tvec, params...);
        auto count_ineq_v1 = count_cgf()->ineq_constraint(count_tvec, params...);
        
        Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> ineq_constraint_vec(summand_ineq_vec.size()+count_ineq_v1.size());

        ineq_constraint_vec << count_ineq_v1, summand_ineq_vec;

        return ineq_constraint_vec;
    }



};

} //namespace CGFs_via_templates


namespace CGFs_with_AD {

using RandomlyStoppedSumCGF = CGF_with_AD_from_template< CGFs_via_templates::RandomlyStoppedSumCGF<CGF_with_AD*, CGF_with_AD*> >;

} // namespace CGFs_with_AD
} // namespace saddlepoint

//////TEST CODE
//// This code creates an instance of the RandomlyStoppedSumCGF.
//// The count variable N follows a Poisson distribution with parameter lambda = 9,
//   // and the summand variable follows Binomial(n = 10, p = 0.64)
//// The Poisson and Binomial CGFs are created first, and then an instance of the RSS CGF is created by passing in these CGFs.
// CGF_with_AD* countCGF = new PoissonCGF();
// CGF_with_AD* summandCGF = new BinomialCGF();
// Adaptor* countAdaptor = new SubvectorAdaptor(0, 1);
// CGF_with_AD* countModelCGF = new parametric_submodelCGF(countCGF, countAdaptor);
// Adaptor* summandAdaptor = new SubvectorAdaptor(1, 2);
// CGF_with_AD* summandModelCGF = new parametric_submodelCGF(summandCGF, summandAdaptor);
// CGF_with_AD* rssCGF = new RandomlyStoppedSumCGF(countModelCGF, summandModelCGF);
//
//// The CGF is then evaluated at various points using the K, K1, K2, K2operator, K3operator, and K4operator functions.
// Eigen::VectorXd tvec(1), parameter_vector(3), vec_of_one(1);
// tvec << 0.34; parameter_vector << 9, 10, 0.64; vec_of_one << 1;
// std::cout << "K(tvec, parameter_vector) = " << rssCGF->K(tvec, parameter_vector) << std::endl;
//// K_Y = -lambda + lambda*(1-p+p*exp(tvec))^n = 81.17247
// std::cout << "K1(tvec, parameter_vector) = " << rssCGF->K1(tvec, parameter_vector) << std::endl;
//// K1_Y = lambda * n * p * (1 + p * (exp(tvec) - 1))^(n - 1) * exp(tvec) = 643.9185
// std::cout << "K2(tvec, parameter_vector) = " << rssCGF->K2(tvec, parameter_vector) << std::endl; // K2_Y = 4782.299
// std::cout << "K2operator(tvec, vec_of_one, vec_of_one, parameter_vector) = " << rssCGF->K2operator(tvec, vec_of_one, vec_of_one, parameter_vector) << std::endl;
//// K2operator is the same as K2_Y
// std::cout << "K3operator(tvec, vec_of_one, vec_of_one, vec_of_one, parameter_vector) = " << rssCGF->K3operator(tvec, vec_of_one, vec_of_one, vec_of_one, parameter_vector) << std::endl;
//// K3_Y = 36700.68
// std::cout << "K4operator(tvec, vec_of_one, vec_of_one, vec_of_one, vec_of_one, parameter_vector) = " << rssCGF->K4operator(tvec, vec_of_one, vec_of_one, vec_of_one, vec_of_one, parameter_vector) << std::endl;
//// K4_Y = 289639.2
//
//// Note: Adaptor created here is used to make each CGF aware of the positions of its parameters within the overall parameter_vector.

#endif // RANDOMLYSTOPPEDSUMCGF_H
