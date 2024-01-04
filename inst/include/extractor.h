#ifndef EXTRACTOR_H_INCLUDED
#define EXTRACTOR_H_INCLUDED

#include "baseCGF.h"

namespace saddlepoint {
namespace CGFs_with_AD {

template <class Templated_CGF, class Extractor>
class CGF_with_AD_From_Extractor_And_Template : public CGF_with_AD, protected Templated_CGF, Extractor {
// Templated CGFs often expect their parameters passed separately, whereas CGF_with_AD passes a single parameter vector
// The Extractor class takes a single parameter vector and returns a pair or tuple of parameters
// The CGF_with_AD_From_Extractor_And_Template class takes a parameter vector, passes it to Extractor, then passes 
// Currently this is implemented only for Extractor classes that return a std::pair
    // Note: it would be natural to use a std::tuple of arbitrary size, but some of the relevant C++ features are in 
    // versions of C++ that are not currently supported via Rcpp
public:
    // Templated_CGF is a protected base class to avoid accidental clashes of the names K, K1, etc
    // but can be accessed via base() method
    typedef Templated_CGF base_CGF_type;
    base_CGF_type& base() {return static_cast<base_CGF_type&>(*this);}
    const base_CGF_type& base() const {return static_cast<const base_CGF_type&>(*this);}
    
    typedef Extractor extractor_type;
    extractor_type& extractor() {return static_cast<extractor_type&>(*this);}
    const extractor_type& extractor() const {return static_cast<const extractor_type&>(*this);}
    
    // Constructors: up to two arguments may be supplied, the first for the templated CGF and the second for the extractor
    // If no argument provided, the corresponding base class is default-initialised
    CGF_with_AD_From_Extractor_And_Template() = default;
    template <class single_arg_type>
    explicit CGF_with_AD_From_Extractor_And_Template(single_arg_type&& arg) : base_CGF_type(std::forward<single_arg_type>(arg)) {}
    template <class single_arg_CGF_type, class single_arg_Extractor_type>
    CGF_with_AD_From_Extractor_And_Template(single_arg_CGF_type&& cgf_arg, single_arg_Extractor_type&& extractor_arg)
        : base_CGF_type(std::forward<single_arg_CGF_type>(cgf_arg)), extractor_type(std::forward<single_arg_Extractor_type>(extractor_arg)) {}
    
protected:
    // Helper functions use overload resolution to process the return type of extractor()::operator()(parameter_vector)
    // Currently only std::pair is supported but it would be natural to support std::tuple<...> as well
    template <class extracted_parameters_type>
    double helper_K(const vec& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    vec helper_K1(const vec& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    mat helper_K2(const vec& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    double helper_K4operatorAABB(const vec& tvec, const extracted_parameters_type& eps, const mat& Q1, const mat& Q2) const;
    template <class extracted_parameters_type>
    double helper_K3K3operatorAABBCC(const vec& tvec, const extracted_parameters_type& eps, const mat& Q1, const mat& Q2, const mat& Q3) const;
    template <class extracted_parameters_type>
    double helper_K3K3operatorABCABC(const vec& tvec, const extracted_parameters_type& eps, const mat& Q1, const mat& Q2, const mat& Q3) const;
    template <class extracted_parameters_type>
    double helper_tilting_exponent(const vec& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    double helper_neg_ll(const vec& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    double helper_func_T(const vec& tvec, const extracted_parameters_type& eps) const;
    // Currently only instantiated for std::pair
    template <class p1_type, class p2_type>
    double helper_K(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    vec helper_K1(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K1(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    mat helper_K2(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K2(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    double helper_K4operatorAABB(const vec& tvec, const std::pair<p1_type, p2_type>& p, const mat& Q1, const mat& Q2) const {return base().K4operatorAABB(tvec, p.first, p.second, Q1, Q2);}
    template <class p1_type, class p2_type>
    double helper_K3K3operatorAABBCC(const vec& tvec, const std::pair<p1_type, p2_type>& p, const mat& Q1, const mat& Q2, const mat& Q3) const {return base().K3K3operatorAABBCC(tvec, p.first, p.second, Q1, Q2, Q3);}
    template <class p1_type, class p2_type>
    double helper_K3K3operatorABCABC(const vec& tvec, const std::pair<p1_type, p2_type>& p, const mat& Q1, const mat& Q2, const mat& Q3) const {return base().K3K3operatorABCABC(tvec, p.first, p.second, Q1, Q2, Q3);}
    template <class p1_type, class p2_type>
    double helper_tilting_exponent(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().tilting_exponent(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    double helper_neg_ll(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().neg_ll(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    double helper_func_T(const vec& tvec, const std::pair<p1_type, p2_type>& p) const {return base().func_T(tvec, p.first, p.second);}
    
    template <class extracted_parameters_type>
    a_scalar helper_K(const a_vector& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    a_vector helper_K1(const a_vector& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    a_matrix helper_K2(const a_vector& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    a_scalar helper_K4operatorAABB(const a_vector& tvec, const extracted_parameters_type& eps, const a_matrix& Q1, const a_matrix& Q2) const;
    template <class extracted_parameters_type>
    a_scalar helper_K3K3operatorAABBCC(const a_vector& tvec, const extracted_parameters_type& eps, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3) const;
    template <class extracted_parameters_type>
    a_scalar helper_K3K3operatorABCABC(const a_vector& tvec, const extracted_parameters_type& eps, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3) const;
    template <class extracted_parameters_type>
    a_scalar helper_tilting_exponent(const a_vector& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    a_scalar helper_neg_ll(const a_vector& tvec, const extracted_parameters_type& eps) const;
    template <class extracted_parameters_type>
    a_scalar helper_func_T(const a_vector& tvec, const extracted_parameters_type& eps) const;
    // Currently only instantiated for std::pair
    template <class p1_type, class p2_type>
    a_scalar helper_K(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    a_vector helper_K1(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K1(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    a_matrix helper_K2(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().K2(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    a_scalar helper_K4operatorAABB(const a_vector& tvec, const std::pair<p1_type, p2_type>& p, const a_matrix& Q1, const a_matrix& Q2) const {return base().K4operatorAABB(tvec, p.first, p.second, Q1, Q2);}
    template <class p1_type, class p2_type>
    a_scalar helper_K3K3operatorAABBCC(const a_vector& tvec, const std::pair<p1_type, p2_type>& p, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3) const {return base().K3K3operatorAABBCC(tvec, p.first, p.second, Q1, Q2, Q3);}
    template <class p1_type, class p2_type>
    a_scalar helper_K3K3operatorABCABC(const a_vector& tvec, const std::pair<p1_type, p2_type>& p, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3) const {return base().K3K3operatorABCABC(tvec, p.first, p.second, Q1, Q2, Q3);}
    template <class p1_type, class p2_type>
    a_scalar helper_tilting_exponent(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().tilting_exponent(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    a_scalar helper_neg_ll(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().neg_ll(tvec, p.first, p.second);}
    template <class p1_type, class p2_type>
    a_scalar helper_func_T(const a_vector& tvec, const std::pair<p1_type, p2_type>& p) const {return base().func_T(tvec, p.first, p.second);}
    
public:
    double K(const vec& tvec, const vec& parameter_vector) const override {return helper_K(tvec, extractor()(parameter_vector));}
    vec K1(const vec& tvec, const vec& parameter_vector) const override {return helper_K1(tvec, extractor()(parameter_vector));}
    mat K2(const vec& tvec, const vec& parameter_vector) const override {return helper_K2(tvec, extractor()(parameter_vector));}
    double K4operatorAABB(const vec& tvec, const mat& Q1, const mat& Q2, const vec& parameter_vector) const override {return helper_K4operatorAABB(tvec, extractor()(parameter_vector), Q1, Q2);}
    double K3K3operatorAABBCC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {return helper_K3K3operatorAABBCC(tvec, extractor()(parameter_vector), Q1, Q2, Q3);}
    double K3K3operatorABCABC(const vec& tvec, const mat& Q1, const mat& Q2, const mat& Q3, const vec& parameter_vector) const override {return helper_K3K3operatorABCABC(tvec, extractor()(parameter_vector), Q1, Q2, Q3);}
    double tilting_exponent(const vec& tvec, const vec& parameter_vector) const override {return helper_tilting_exponent(tvec, extractor()(parameter_vector));}
    double neg_ll(const vec& tvec, const vec& parameter_vector) const override {return helper_neg_ll(tvec, extractor()(parameter_vector));}
    double func_T(const vec& tvec, const vec& parameter_vector) const override {return helper_func_T(tvec, extractor()(parameter_vector));}
    
    a_scalar K(const a_vector& tvec, const a_vector& parameter_vector) const override {return helper_K(tvec, extractor()(parameter_vector));}
    a_vector K1(const a_vector& ta_vector, const a_vector& parameter_vector) const override {return helper_K1(ta_vector, extractor()(parameter_vector));}
    a_matrix K2(const a_vector& ta_vector, const a_vector& parameter_vector) const override {return helper_K2(ta_vector, extractor()(parameter_vector));}
    a_scalar K4operatorAABB(const a_vector& ta_vector, const a_matrix& Q1, const a_matrix& Q2, const a_vector& parameter_vector) const override {return helper_K4operatorAABB(ta_vector, extractor()(parameter_vector), Q1, Q2);}
    a_scalar K3K3operatorAABBCC(const a_vector& ta_vector, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {return helper_K3K3operatorAABBCC(ta_vector, extractor()(parameter_vector), Q1, Q2, Q3);}
    a_scalar K3K3operatorABCABC(const a_vector& ta_vector, const a_matrix& Q1, const a_matrix& Q2, const a_matrix& Q3, const a_vector& parameter_vector) const override {return helper_K3K3operatorABCABC(ta_vector, extractor()(parameter_vector), Q1, Q2, Q3);}
    a_scalar tilting_exponent(const a_vector& ta_vector, const a_vector& parameter_vector) const override {return helper_tilting_exponent(ta_vector, extractor()(parameter_vector));}
    a_scalar neg_ll(const a_vector& ta_vector, const a_vector& parameter_vector) const override {return helper_neg_ll(ta_vector, extractor()(parameter_vector));}
    a_scalar func_T(const a_vector& ta_vector, const a_vector& parameter_vector) const override {return helper_func_T(ta_vector, extractor()(parameter_vector));}
    
    // Scratchpad declarations may be overridden if desired
    CGF_base<double, vec, mat>::Scratchpad* scratchpad(const vec& tvec, const vec& parameter_vector) const;
    CGF_base<a_scalar, a_vector, a_matrix>::Scratchpad* scratchpad(const a_vector& tvec, const a_vector& parameter_vector) const;
};

} // namespace CGFs_with_AD
} // namespace saddlepoint

#endif
