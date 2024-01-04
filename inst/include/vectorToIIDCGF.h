#ifndef VECTORTOIIDCGF_H_INCLUDED
#define VECTORTOIIDCGF_H_INCLUDED

#include "CGF_Defaults.h"
#include "baseCGF.h"
#include "saddlepoint_types.h"
#include "BaseWrapper.h"
#include "parametric_submodelCGF.h"

namespace saddlepoint {
namespace CGFs_via_templates {

// IIDReplicates_CGF class takes a vector base CGF as input and creates an (iid) based CGF.
// The block_size argument specifies the size of the vector to be replicated
// based on the available data (observations).
template <class BaseCGF>
class IIDReplicates_CGF : public CGF_Defaults<IIDReplicates_CGF<BaseCGF>>, private BaseWrapper<BaseCGF> {
private:
    using BaseWrapper<BaseCGF>::base_cgf;
public:
    using BaseWrapper<BaseCGF>::BaseWrapper; // Provides access to constructors for BaseCGF
    //---------------------
    template <class t_vector_type, class block_size_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const block_size_type& block_size, ParamTypes&&... other_params) const {
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));
        const auto m = tvec.size() / block_size;
        typedef decltype(base_cgf()->K(tvec.segment(0, block_size), std::forward<ParamTypes>(other_params)...)) scalar_type;
        scalar_type res = 0;
        for (auto i = 0; i < m; i++) {
            res += base_cgf()->K(tvec.segment(i*block_size, block_size), std::forward<ParamTypes>(other_params)...);
        }
        return res;
    }
    //---------------------
    template <class t_vector_type, class block_size_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const block_size_type& block_size, ParamTypes&&... other_params) const {
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));
        const auto m = tvec.size() / block_size;
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(base_cgf()->K1(tvec.segment(0, block_size), std::forward<ParamTypes>(other_params)...).eval())>> vector_type;
        vector_type res = vector_type::Zero(m * block_size);
        for (auto i = 0; i < m; i++) {
            res.segment(i*block_size, block_size) = base_cgf()->K1(tvec.segment(i*block_size, block_size), std::forward<ParamTypes>(other_params)...);
        }
        return res;
    }
    //---------------------
    template <class t_vector_type, class block_size_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const block_size_type& block_size, ParamTypes&&... other_params) const {
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));

        const auto m = tvec.size() / block_size;
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(base_cgf()->K2(tvec.segment(0, block_size), std::forward<ParamTypes>(other_params)...).eval())>> matrix_type;
        matrix_type res = matrix_type::Zero(m*block_size, m*block_size);
        for (auto i = 0; i < m; ++i) {
                res.block(i*block_size, i*block_size, block_size, block_size) = base_cgf()->K2(tvec.segment(i*block_size, block_size), std::forward<ParamTypes>(other_params)...);
        }
        return res;
    }
    //---------------------
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class block_size_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const block_size_type& block_size, ParamTypes&&... other_params) const {
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));

        const auto m = tvec.size() / block_size;
        typedef decltype(base_cgf()->K3operator(tvec.segment(0, block_size), v1.segment(0, block_size), v2.segment(0, block_size), v3.segment(0, block_size), std::forward<ParamTypes>(other_params)...)) scalar_type;
        scalar_type res = 0;
        for (auto i = 0; i < m; i++) res += base_cgf()->K3operator(tvec.segment(i * block_size, block_size), v1.segment(i * block_size, block_size), v2.segment(i * block_size, block_size), v3.segment(i * block_size, block_size), std::forward<ParamTypes>(other_params)...);
        return res;
    }
    //---------------------
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class block_size_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const v4_type& v4, const block_size_type& block_size, ParamTypes&&... other_params) const {
        const auto m = tvec.size() / block_size;
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));
        typedef decltype(base_cgf()->K4operator(tvec.segment(0,block_size), v1.segment(0,block_size), v2.segment(0,block_size), v3.segment(0,block_size), v4.segment(0,block_size), std::forward<ParamTypes>(other_params)...)) scalar_type;
        scalar_type res = 0;
        for (auto i = 0; i < m; i++) res += base_cgf()->K4operator(tvec.segment(i*block_size,block_size), v1.segment(i*block_size,block_size), v2.segment(i*block_size,block_size), v3.segment(i*block_size,block_size), v4.segment(i*block_size,block_size), std::forward<ParamTypes>(other_params)...);
        return res;
    }
    //---------------------
    template <class t_vector_type, class block_size_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const block_size_type& block_size, ParamTypes&&... other_params) const {
        const auto m = tvec.size() / block_size;
        if (tvec.size() % block_size != 0) throw std::invalid_argument("Dimension mismatch: The length of data is not divisible by block_size. Expected size is a multiple of " + std::to_string(block_size));
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;

        auto first_ineq_vec = base_cgf()->ineq_constraint(tvec.segment(0,block_size), std::forward<ParamTypes>(other_params)...);
        auto j = first_ineq_vec.size();
        vector_type res(m*j);

        if(j > 0){
            for (auto i = 0; i < m; i++) {
                res.segment(i*j,j) = base_cgf()->ineq_constraint(tvec.segment(i*block_size,block_size), std::forward<ParamTypes>(other_params)...);
            }
        }
        return res;
    }
    //---------------------
    //---------------------

};
} // namespace CGFs_via_templates

namespace CGFs_with_AD {



class Provide_block_size {
private:
    int block_size;
public:
    explicit Provide_block_size(const int& block_size_) : block_size(block_size_) {
        if (block_size <= 0) throw std::invalid_argument("Invalid arguments: block_size must be a positive integer.");
    }

    Provide_block_size() = delete;


    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const vec& parameter_vector) const {
        return co(block_size, parameter_vector);
    }
    template <class CallbackObject>
    auto operator()(const CallbackObject& co, const a_vector& parameter_vector) const {
        return co(block_size, parameter_vector);
    }

};


using IIDReplicates_CGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::AdaptedParametersCallbackCGF<saddlepoint::CGFs_via_templates::IIDReplicates_CGF<CGF_with_AD*>, Provide_block_size>>;

} // namespace CGFs_with_AD



} // namespace saddlepoint


#endif // VECTORTOIIDCGF_H_INCLUDED



//// TEST CODE
//// The simple MVPoisson code Y = X + AZ_0 // the cgf of Y corresponds to sum_independent_cgf below
//// Create a Poisson CGF
//CGF_with_AD* pois_cgf = new PoissonCGF();
//// Create an adaptor for the first CGF, where lambda is the first parameter in the model parameter theta = c(14,7)
//Adaptor* first_adaptor = new SubvectorAdaptor(0, 1);
//// Create a parametric submodel CGF for the first Poisson model CGF with the given adaptor
//CGF_with_AD* first_poisson_model_cgf = new parametric_submodelCGF(pois_cgf, first_adaptor);
//// Create a 3x1 Eigen Matrix with all values set to 1
//Eigen::MatrixXd A(3, 1);
//A.setOnes();
//// Create a second Poisson model CGF with the lambda parameter in the second position
//CGF_with_AD* second_poisson_model_cgf = new parametric_submodelCGF(new PoissonCGF(), new SubvectorAdaptor(1, 1));
//// Create a linearly mapped CGF from the second Poisson model CGF using the A matrix
//CGF_with_AD* linear_mapped_cgf = new LinearlyMappedCGF(second_poisson_model_cgf, A);
//// Create a vector of WrapAsCGF objects that wrap the CGFs for the two Poisson models
//std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
//cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(first_poisson_model_cgf));
//cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(linear_mapped_cgf));
//// Create a sum of independent CGF from the vector of wrapped CGFs
//CGF_with_AD* sum_independent_cgf = new SumOfIndependentCGF(cgfs);
//
//// set vector Y to be of length 3 and make a cgf for 2-iid copies of Y
//CGF_with_AD* vector_to_iid_cgf = new IIDReplicates_CGF(sum_independent_cgf, 2);
//Eigen::VectorXd tvec(6), parameter_vector(2);
//tvec << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
//parameter_vector << 14, 7;
//std::cout << "K1 = " << vector_to_iid_cgf->K1(tvec, parameter_vector) << std::endl; // Output: 28.22722 29.85447 31.65285 52.25737 54.45392 56.88149
//Eigen::VectorXd v1(6), v2(6), v3(6), v4(6);
//v1 << 0.45, 0.65, 0.3, 0.45, 2.4, 0.45;
//v2 << 0.44, 0.1, 0.2, 0.3, 0.67, 0.1;
//v3 << -0.45, 0.123, 1.5, 0.68, 0.004, 1.2;
//v4 << 0.32, -0.65, -0.006, -1.5, 0.4, 0.45;
//std::cout << "K4operator " << vector_to_iid_cgf->K4operator(tvec, v1, v2, v3, v4, parameter_vector) << std::endl; // Output: -143.599







//// Using the following as an alternative should yield the same CGF
//CGF_with_AD* pois_cgf = new PoissonCGF();
//Adaptor* first_adaptor = new SubvectorAdaptor(0, 1);
//CGF_with_AD* first_poisson_model_cgf = new parametric_submodelCGF(pois_cgf, first_adaptor);
//int d = 3, m = 2;
//Eigen::MatrixXd A(m * d, m);
//A.setZero();
//for (int i = 0; i < m; ++i) {
//    A.block(i * d, i, d, 1).setOnes();
//}
//CGF_with_AD* second_poisson_model_cgf = new parametric_submodelCGF(new PoissonCGF(), new SubvectorAdaptor(1, 1));
//CGF_with_AD* linear_mapped_cgf = new LinearlyMappedCGF(second_poisson_model_cgf, A);
//std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
//cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(first_poisson_model_cgf));
//cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(linear_mapped_cgf));
//CGF_with_AD* sum_independent_cgf = new SumOfIndependentCGF(cgfs);
//
//
//Eigen::VectorXd tvec(6), parameter_vector(2);
//tvec << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
//parameter_vector << 14, 7;
//Eigen::VectorXd v1(6), v2(6), v3(6), v4(6);
//v1 << 0.45, 0.65, 0.3, 0.45, 2.4, 0.45;
//v2 << 0.44, 0.1, 0.2, 0.3, 0.67, 0.1;
//v3 << -0.45, 0.123, 1.5, 0.68, 0.004, 1.2;
//v4 << 0.32, -0.65, -0.006, -1.5, 0.4, 0.45;
//std::cout << "K4operator " << sum_independent_cgf->K4operator(tvec, v1, v2, v3, v4, parameter_vector) << std::endl;
