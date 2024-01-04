#ifndef SUMOFINDEPENDENTCGF_H
#define SUMOFINDEPENDENTCGF_H

#include "CGF_Defaults.h"
// other include directives needed for CGF_with_AD version only
#include "baseCGF.h"
#include "parametric_submodelCGF.h"
#include <vector>
#include <type_traits>

namespace saddlepoint {
namespace CGFs_via_templates {

#include "needs_full_list_of_CGF_methods.h"
template <class Container>
class SumOfIndependentCGF : public CGF_Defaults<SumOfIndependentCGF<Container>> {
    // The Container class should provide an interface similar to the standard container classes
    // with elements of a suitable CGF type.
    // To use CGF_with_AD as the CGF type, use the WrapAsCGF class, as in
    // SumOfIndependentCGF<std::vector<WrapAsCGF<CGF_with_AD*>>>
private:
    Container cgfs;
public:
    template <class... ArgTypes>
    explicit SumOfIndependentCGF(ArgTypes&&... args) : cgfs(std::forward<ArgTypes>(args)...) {}

    //typedef Container container_type;

    template <class t_vector_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const ParamTypes&... params) const {
        // All arguments are passed as const reference, not using std::forward, since they will be passed to every CGF in cgfs
        typedef decltype(cgfs.begin()->K(tvec, params...)) scalar_type;
        scalar_type res = 0;
        // If cgfs is empty, the return value will be 0
        // This correctly corresponds to the CGF of an empty sum
        // (i.e., the CGF of the constant random variable whose value is the 0 vector of size matching tvec)
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {res += i->K(tvec, params...);}
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const ParamTypes&... params) const {
        // Extract the return type (stripped of const and reference qualifiers)
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(cgfs.begin()->K1(tvec, params...).eval())>> vector_type;

        // We assume that vector_type follows the Eigen-style naming convention
        // with a static method Zero(n) returning a zero vector of size n
        vector_type res = vector_type::Zero(tvec.size());
        // If cgfs is empty, the return value will be the 0 vector of the size matching tvec
        // This correctly corresponds to the mean of an empty sum, interpreted as the 0 vector of size matching tvec
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {res += i->K1(tvec, params...);}
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const ParamTypes&... params) const {
        // Extract the return type (stripped of const and reference qualifiers)
        typedef std::remove_cv_t<typename std::remove_reference_t<decltype(cgfs.begin()->K2(tvec, params...).eval())>> matrix_type;
        // We assume that matrix_type follows the Eigen-style naming convention
        // with a static method Zero(n, m) returning a zero matrix of size n-by-m
        matrix_type res = matrix_type::Zero(tvec.size(), tvec.size());
        // If cgfs is empty, the return value will be the 0 vector of the size matching tvec
        // This correctly corresponds to the mean of an empty sum, interpreted as the 0 vector of size matching tvec
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {res += i->K2(tvec, params...);}
        return res;
    }

    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->K3operator(tvec, v1, v2, v3, params...)) scalar_type;
        scalar_type res = 0;
        // If cgfs is empty, the return value will be 0
        // This correctly corresponds to the CGF of an empty sum
        // (i.e., the CGF of the constant random variable whose value is the 0 vector of size matching tvec)
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {res += i->K3operator(tvec, v1, v2, v3, params...);}
        return res;
    }
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const v4_type& v4,
                    const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->K4operator(tvec, v1, v2, v3, v4, params...)) scalar_type;
        scalar_type res = 0;
        // If cgfs is empty, the return value will be 0
        // This correctly corresponds to the CGF of an empty sum
        // (i.e., the CGF of the constant random variable whose value is the 0 vector of size matching tvec)
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {res += i->K4operator(tvec, v1, v2, v3, v4, params...);}
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const ParamTypes&... params) const {
        // Define vector_type and inner_vector_type based on the input type and the return type of internal ineq_constraints
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        typedef decltype(cgfs.begin()->ineq_constraint(tvec, params...)) inner_vector_type;

        // Calculate the total size needed for the result
        int total_size = 0;
        for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {
            inner_vector_type ineq_vec = i->ineq_constraint(tvec, params...);
            total_size += ineq_vec.size();
        }

        // Allocate the size for res
        vector_type res(total_size);

        // If total_size is zero, it means that all ineq_constraint calls returned empty vectors, so we can return an empty vector.
        // If it's greater than zero, we need to copy the results into res
        if (total_size > 0) {
            // Copy the results into res by iterating through cgfs
            // and extracting the ineq_constraint from each element
            int current_index = 0;
            for (auto i = cgfs.begin(); i != cgfs.end(); ++i) {
                inner_vector_type ineq_vec = i->ineq_constraint(tvec, params...);
                res.segment(current_index, ineq_vec.size()) = ineq_vec;
                current_index += ineq_vec.size();
            }
        }

        return res;
    }
    // TO DO: implement other methods.
    // Currently these are being implemented by CGF_Defaults
};

} // namespace CGFs_via_templates

namespace CGFs_with_AD {

using SumOfIndependentCGF = CGF_with_AD_from_template<CGFs_via_templates::SumOfIndependentCGF<std::vector<CGFs_via_templates::WrapAsCGF<CGF_with_AD*>>>>;
} // namespace CGFs_with_AD
} // namespace saddlepoint

#endif // SUMOFINDEPENDENTCGF_H



//// TEST CODE
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
//CGF_with_AD* final_cgf = new SumOfIndependentCGF(cgfs);
//// Create a 3x1 Eigen Vector with given values
//Eigen::VectorXd tvec(3), parameter_vector(2);
//tvec << 0.1, 0.2, 0.3;
//parameter_vector << 14, 7;
//// Compute and print K for the final CGF with the given tvec and parameter_vector
//std::cout << "K " << final_cgf->K(tvec, parameter_vector) << std::endl; // Output: 15.22489

