#ifndef CONCATENATIONCGF_H
#define CONCATENATIONCGF_H

#include "CGF_Defaults.h"
// other include directives needed for CGF_with_AD version only
#include "baseCGF.h"
#include "parametric_submodelCGF.h"
#include <vector>
#include <type_traits>

namespace saddlepoint {
namespace CGFs_via_templates {

template <class Container>
class ConcatenationCGF : public CGF_Defaults<ConcatenationCGF<Container>> {
private:
    Container cgfs;
    Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> vecLengths;

public:
    template <class ContainerType>
    explicit ConcatenationCGF(ContainerType&& container, const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& vecLengths_)
            : cgfs(std::forward<ContainerType>(container)), vecLengths(vecLengths_) {}

    template <class t_vector_type, class... ParamTypes>
    auto K(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->K(tvec.segment(0, vecLengths[0]), params...)) scalar_type;
        scalar_type res = 0;
        //Eigen::Index current_start = 0;
        //if (vecLengths.sum() != tvec.size()) throw std::invalid_argument("Mismatch between vector length and data size in ConcatenationCGF::K");
        //for(Eigen::Index i = 0; i < cgfs.size(); ++i) {
        //    res += cgfs[i].K(tvec.segment(current_start, vecLengths[i]), params...);
        //    current_start += vecLengths[i];
        //}
        // if (vecLengths.sum() != tvec.size()) throw std::invalid_argument("Mismatch between vector length and data size in ConcatenationCGF");
        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res += it->K(tvec.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }

        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto K1(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(cgfs.begin()->K1(tvec.segment(0, vecLengths[0]), params...).eval())>> vector_type;
        vector_type res = vector_type::Zero(tvec.size());

        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res.segment(current_start, vecLengths[i]) += it->K1(tvec.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto K2(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef std::remove_cv_t<typename std::remove_reference_t<decltype(cgfs.begin()->K2(tvec.segment(0, vecLengths[0]), params...).eval())>> matrix_type;
        matrix_type res = matrix_type::Zero(tvec.size(), tvec.size());

        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res.block(current_start, current_start, vecLengths[i], vecLengths[i]) += it->K2(tvec.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }
        return res;
    }

    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class... ParamTypes>
    auto K3operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->K3operator(tvec.segment(0, vecLengths[0]), v1.segment(0, vecLengths[0]), v2.segment(0, vecLengths[0]), v3.segment(0, vecLengths[0]), params...)) scalar_type;
        scalar_type res = 0;

        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res += it->K3operator(tvec.segment(current_start, vecLengths[i]),
                                  v1.segment(current_start, vecLengths[i]),
                                  v2.segment(current_start, vecLengths[i]),
                                  v3.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }
        return res;
    }

    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class... ParamTypes>
    auto K4operator(const t_vector_type& tvec, const v1_type& v1, const v2_type& v2, const v3_type& v3, const v4_type& v4,
                    const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->K4operator(tvec.segment(0, vecLengths[0]), v1.segment(0, vecLengths[0]), v2.segment(0, vecLengths[0]),
                                    v3.segment(0, vecLengths[0]), v4.segment(0, vecLengths[0]), params...)) scalar_type;
        scalar_type res = 0;

        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res += it->K4operator(tvec.segment(current_start, vecLengths[i]),
                                  v1.segment(current_start, vecLengths[i]),
                                  v2.segment(current_start, vecLengths[i]),
                                  v3.segment(current_start, vecLengths[i]),
                                  v4.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto func_T(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef decltype(cgfs.begin()->func_T(tvec.segment(0, vecLengths[0]), params...)) scalar_type;
        scalar_type res = 0;

        Eigen::Index current_start = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            res += it->func_T(tvec.segment(current_start, vecLengths[i]), params...);
            current_start += vecLengths[i];
        }
        return res;
    }

    template <class t_vector_type, class... ParamTypes>
    auto ineq_constraint(const t_vector_type& tvec, const ParamTypes&... params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
        typedef decltype(cgfs.begin()->ineq_constraint(tvec.segment(0, vecLengths[0]), params...)) inner_vector_type;

        std::vector<typename inner_vector_type::Scalar> vals;
        Eigen::Index current_index = 0, i = 0;
        for(auto it = cgfs.begin(); it != cgfs.end(); ++it, ++i) {
            inner_vector_type ineq_vec = it->ineq_constraint(tvec.segment(current_index, vecLengths[i]), params...);
            vals.insert(vals.end(), ineq_vec.data(), ineq_vec.data() + ineq_vec.size());
            current_index += vecLengths[i];
        }

        vector_type res = Eigen::Map<vector_type>(vals.data(), vals.size());
        return res;
    }

//            template <class t_vector_type, class... ParamTypes>
//            auto ineq_constraint(const t_vector_type& tvec, const ParamTypes&... params) const {
//                typedef typename std::remove_cv_t<typename std::remove_reference_t<decltype(tvec.eval())>> vector_type;
//                typedef decltype(cgfs.begin()->ineq_constraint(tvec.segment(0, vecLengths[0]), params...)) inner_vector_type;
//
//                int total_size = 0;
//                int current_index0 = 0;
//                for(Eigen::Index i = 0; i < cgfs.size(); ++i) {
//                    inner_vector_type ineq_vec = cgfs[i]->ineq_constraint(tvec.segment(current_index0, vecLengths[i]), params...);
//                    current_index0 += vecLengths[i];
//                    total_size += ineq_vec.size();
//                }
//
//                vector_type res(total_size);
//
//                if (total_size > 0) {
//                    int current_index = 0;
//                    for(Eigen::Index i = 0; i < cgfs.size(); ++i) {
//                        inner_vector_type ineq_vec = cgfs[i]->ineq_constraint(tvec.segment(current_index, vecLengths[i]), params...);
//                        res.segment(current_index, ineq_vec.size()) = ineq_vec;
//                        current_index += ineq_vec.size();
//                    }
//                }
//
//                return res;
//            }


};

} // namespace CGFs_via_templates

namespace CGFs_with_AD {
using ConcatenationCGF = CGF_with_AD_from_template<CGFs_via_templates::ConcatenationCGF< std::vector<CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> >>;
// To create an instance of ConcatenationCGF, use the following syntax:
    // new ConcatenationCGF(ContainerType, Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>);
    // where ContainerType is a container type such as std::vector - see example below
    // and the Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> is a vector of the same length as the container - see example below
}// namespace CGFs_with_AD



} // namespace saddlepoint


////  TEST CODE
////  Instantiate a Poisson model CGF. An adaptor to point to the position of the distributional paramemter in the model parameter.
// Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices1(1); indices1 << 0;  // Set the index to 0
// CGF_with_AD* pois_cgf = new PoissonModelCGF( new ScalarAdaptorFromVectorAdaptor(new SubsetVectorByIndicesAdaptor(indices1)) );
////  Instantiate a Binomial model CGF with appropriate adaptors for 'n' and 'p'
// Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices_n(1), indices_p(1);
// indices_n << 1; indices_p << 2;
// CGF_with_AD* binom_cgf = new BinomialModelCGF(new ScalarAdaptorFromVectorAdaptor(new SubsetVectorByIndicesAdaptor(indices_n)),
//                                               new ScalarAdaptorFromVectorAdaptor(new SubsetVectorByIndicesAdaptor(indices_p)));
////  A container for CGFs containing the Poisson and Binomial CGFs
// std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> cgfs;
// cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(pois_cgf));
// cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(binom_cgf));

////  vector lengths
// Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> vecLengths(2);
// vecLengths << 3, 2;
//// vecLengths.sum() must correspond to the length of tvec passed to the K, K1, K2, etc. methods of ConcatenationCGF

//// An instance of ConcatenationCGF with the CGFs container and the vector lengths
//CGF_with_AD* concat_cgf = new ConcatenationCGF(cgfs, vecLengths);

//// Define test vectors
//Eigen::VectorXd tvec(5), parameter_vector(3);
//tvec << 0.1, 0.2, 0.3, 0.4, 0.5;
//parameter_vector << 12, 10, 0.64;
//std::cout << "ConcatenationCGF - K: " << concat_cgf->K(tvec, parameter_vector) << std::endl;
//
//
//
//
//// TEST CODE: the following code is equivalent to the above (should give the same result)
//Eigen::MatrixXd A0(5, 3);
//A0 << Eigen::MatrixXd::Identity(3, 3), // A 3x3 identity matrix
//      Eigen::MatrixXd::Zero(2, 3); // A 2x3 zero matrix
//Eigen::MatrixXd A1(5, 2);
//A1 << Eigen::MatrixXd::Zero(3, 2), // A 3x2 zero matrix
//      Eigen::MatrixXd::Identity(2, 2); // A 2x2 identity matrix
//
//CGF_with_AD* cgf1 = new LinearlyMappedCGF(new PoissonCGF(), A0);
//CGF_with_AD* cgf2 = new LinearlyMappedCGF(new BinomialCGF(), A1);
//
//Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> cgf1_indices(1), cgf2_indices(2);
//cgf1_indices << 0; cgf2_indices << 1,2;
//CGF_with_AD* model_poisson_cgf = new parametric_submodelCGF( cgf1, new SubsetVectorByIndicesAdaptor(cgf1_indices) );
//CGF_with_AD* model_binomi_cgf = new parametric_submodelCGF( cgf2, new SubsetVectorByIndicesAdaptor(cgf2_indices) );
//
//std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> sumOfIndpendnent_cgf;
//sumOfIndpendnent_cgf.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(model_poisson_cgf));
//sumOfIndpendnent_cgf.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(model_binomi_cgf));
//
//CGF_with_AD* sumOfIndpendnent_CGF = new SumOfIndependentCGF(sumOfIndpendnent_cgf);
//std::cout << "sumOfIndpendnent - K: " << sumOfIndpendnent_CGF->K(tvec, parameter_vector) << std::endl;

//// Equivalence check for the other methods
//// K1
//assert(concat_cgf->K1(tvec, parameter_vector).isApprox(sumOfIndpendnent_CGF->K1(tvec, parameter_vector), 1e-6));
//// K2
//assert(concat_cgf->K2(tvec, parameter_vector).isApprox(sumOfIndpendnent_CGF->K2(tvec, parameter_vector), 1e-6));
//// K3operator and K4operator
//Eigen::VectorXd vec1(5), vec2(5), vec3(5), vec4(5);
//vec1 << 1, 2, 3, 4, 5; vec2 << 2.6, 3.2, 4.1, 5.6, 6.3; vec3 << 3.1, 4, 5, 6, 7.2; vec4 << 8.0, 9.0, 10, 4.2, 1.3;
//assert(std::abs(concat_cgf->K3operator(tvec, vec1, vec2, vec3, parameter_vector) -
//                sumOfIndpendnent_CGF->K3operator(tvec, vec1, vec2, vec3, parameter_vector)) < 1e-6);
//assert(std::abs(concat_cgf->K4operator(tvec, vec1, vec2, vec3, vec4, parameter_vector) -
//                sumOfIndpendnent_CGF->K4operator(tvec, vec1, vec2, vec3, vec4, parameter_vector)) < 1e-6);


////ineq_constraint
//Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> position4_for_a_new_parameter(1); position4_for_a_new_parameter << 3;
//CGF_with_AD* geom_modelCGF = new GeometricModelCGF(new ScalarAdaptorFromVectorAdaptor(new SubsetVectorByIndicesAdaptor(position4_for_a_new_parameter)));
//// A container for CGFs containing the Poisson and Binomial CGFs
//std::vector<saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>> updated_cgfs;
//updated_cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(concat_cgf));
//updated_cgfs.push_back(saddlepoint::CGFs_via_templates::WrapAsCGF<CGF_with_AD*>(geom_modelCGF));
//
//Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> updatedvecLengths(2);
//updatedvecLengths << 5, 4;
//
//CGF_with_AD* updated_concat_cgf = new ConcatenationCGF(updated_cgfs, updatedvecLengths);
//
//Eigen::VectorXd updated_tvec(9), updated_parameter_vector(4);
//updated_tvec << tvec, -5.6, -0.3, -0.1, -0.54;
//updated_parameter_vector << parameter_vector, 0.03;
//
//std::cout << "UpdatedConcatenationCGF - ineq_constraint: " << updated_concat_cgf->ineq_constraint(updated_tvec, updated_parameter_vector) << std::endl;
////results from the linearly mapped approach: -5.6304592 -0.3304592 -0.1304592 -0.5704592

#endif // CONCATENATIONCGF_H



