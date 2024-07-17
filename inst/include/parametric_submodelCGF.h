#ifndef PARAMETRIC_SUBMODELCGF_H_INCLUDED
#define PARAMETRIC_SUBMODELCGF_H_INCLUDED

#include <utility>
#include <tuple>
#include "BaseWrapper.h"
#include "FunctionObjectWrapper.h"
#include "baseCGF.h"
// #include "saddlepoint_types.h"


// TO DO: document more how the constructors work with AdaptedParametersCallbackCGF

namespace saddlepoint {
namespace CGFs_via_templates {

// A general class to make changes to the parameter(s) of a template-style CGF,
// while leaving the tvec argument (and other arguments such as the vectors passed to the K[x]operator methods) unchanged.
// These changes may include changing the number of parameters.
// As a result FunctionObject must satisfy a fairly specific calling convention.
// For simpler special cases use classes defined below instead, with more natural interfaces,
// e.g., use ReparametrisedCGF for transforming one parameter vector into another parameter vector.
template <class BaseCGF, class CallbackFunctionObject>
class AdaptedParametersCallbackCGF : private BaseWrapper<BaseCGF>, FunctionObjectWrapper<CallbackFunctionObject> {
    // CallbackFunctionObject provides an operator() of the form
        // template <class CallbackObject, class... ParameterTypes>
        // auto operator()(const CallbackObject& co, ParameterTypes... params) const;
    // This operator should process params... in an appropriate way, then execute
        // return co(new_param_1, ..., new_param_k);
    // with whatever number of parameters are appropriate.
    // Do not call the CallbackObject more than once.
    // CallbackObject will be such that, when AdaptParameters::[method_name](tvec, other_arguments..., params...) is called,
    // CallbackFunctionObject's return statement will be equivalent to
        // return static_case<BaseCGF*>(this)->[method_name](tvec, other_arguments..., new_param_1, ..., new_param_k)
    // Here other_arguments... means the method-specific arguments, if any,
    // that are required by [method_name] in addition to the tvec argument and the parameter arguments,
    // such as the three vectors passed to K3operator
private:
    // Abbreviations to access BaseCGF and CallbackFunctionObject provided by internal wrapper classes
    using BaseWrapper<BaseCGF>::base_cgf;
    const CallbackFunctionObject* cfo() const {return FunctionObjectWrapper<CallbackFunctionObject>::function_object();}
public:
    // Constructors can accept any number of arguments
    // First argument (if present) initialises BaseCGF
    // Second and additional arguments (if present) initialise CallbackFunctionObject
    AdaptedParametersCallbackCGF() = default;
    template <class BaseCGFArgType>
    explicit AdaptedParametersCallbackCGF(BaseCGFArgType&& bca)
        : BaseWrapper<BaseCGF>(std::forward<BaseCGFArgType>(bca)) {}
    // Single-argument form: CallbackFunctionObject is default-initialised

    template <class BaseCGFArgType, class CFOArgType, class... ExtraCFOArgTypes>
    AdaptedParametersCallbackCGF(BaseCGFArgType&& bca, CFOArgType&& foa, ExtraCFOArgTypes&&... extra_foa)
        : BaseWrapper<BaseCGF>(std::forward<BaseCGFArgType>(bca)),
          FunctionObjectWrapper<CallbackFunctionObject>(std::forward<CFOArgType>(foa), std::forward<ExtraCFOArgTypes>(extra_foa)...) {}
    // With two or more arguments, first argument initialises BaseCGF, remaining argument(s) initialise CallbackFunctionObject

    template <class t_vector_type, class... ParameterTypes>
    auto K(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        // tvec argument is peeled off and not passed to FunctionObject, and saved to be passed to base_cgf()->K
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    // K1, K2, tilting_exponent, neg_ll, func_T are similar
    template <class t_vector_type, class... ParameterTypes>
    auto K1(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K1(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto K2(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K2(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto tilting_exponent(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->tilting_exponent(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto neg_ll(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->neg_ll(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto func_T(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->func_T(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }

    // Other methods have a specified number of additional arguments which are also peeled off
    template <class t_vector_type, class x_type, class y_type, class... ParameterTypes>
    auto K2operator(t_vector_type&& tvec, x_type&& x, y_type&& y, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K2operator(std::forward<t_vector_type>(tvec), std::forward<x_type>(x), std::forward<y_type>(y),
                                      std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class matrix_type, class... ParameterTypes>
    auto K2operatorAK2AT(t_vector_type&& tvec, matrix_type&& A, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K2operatorAK2AT(std::forward<t_vector_type>(tvec), std::forward<matrix_type>(A),
                                           std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class... ParameterTypes>
    auto K3operator(t_vector_type&& tvec, v1_type&& v1, v2_type&& v2, v3_type&& v3, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K3operator(std::forward<t_vector_type>(tvec),
                                      std::forward<v1_type>(v1), std::forward<v2_type>(v2), std::forward<v3_type>(v3),
                                      std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class v1_type, class v2_type, class v3_type, class v4_type, class... ParameterTypes>
    auto K4operator(t_vector_type&& tvec, v1_type&& v1, v2_type&& v2, v3_type&& v3, v4_type&& v4, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K4operator(std::forward<t_vector_type>(tvec),
                                      std::forward<v1_type>(v1), std::forward<v2_type>(v2),
                                      std::forward<v3_type>(v3), std::forward<v4_type>(v4),
                                      std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class Q1_type, class Q2_type, class... ParameterTypes>
    auto K4operatorAABB(t_vector_type&& tvec, Q1_type&& Q1, Q2_type&& Q2, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K4operatorAABB(std::forward<t_vector_type>(tvec),
                                          std::forward<Q1_type>(Q1), std::forward<Q2_type>(Q2),
                                          std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class... ParameterTypes>
    auto K3K3operatorAABBCC(t_vector_type&& tvec, Q1_type&& Q1, Q2_type&& Q2, Q3_type&& Q3, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K3K3operatorAABBCC(std::forward<t_vector_type>(tvec),
                                              std::forward<Q1_type>(Q1), std::forward<Q2_type>(Q2), std::forward<Q3_type>(Q3),
                                              std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class Q1_type, class Q2_type, class Q3_type, class... ParameterTypes>
    auto K3K3operatorABCABC(t_vector_type&& tvec, Q1_type&& Q1, Q2_type&& Q2, Q3_type&& Q3, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K3K3operatorABCABC(std::forward<t_vector_type>(tvec),
                                              std::forward<Q1_type>(Q1), std::forward<Q2_type>(Q2), std::forward<Q3_type>(Q3),
                                              std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class A1_type, class d1_type, class A2_type, class d2_type, class... ParameterTypes>
    auto K4operatorAABB_factored(t_vector_type&& tvec, A1_type&& A1, d1_type&& d1, A2_type&& A2, d2_type&& d2, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K4operatorAABB_factored(std::forward<t_vector_type>(tvec),
                                                   std::forward<A1_type>(A1), std::forward<d1_type>(d1),
                                                   std::forward<A2_type>(A2), std::forward<d2_type>(d2),
                                                   std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class A1_type, class d1_type, class A2_type, class d2_type, class A3_type, class d3_type, class... ParameterTypes>
    auto K3K3operatorAABBCC_factored(t_vector_type&& tvec, A1_type&& A1, d1_type&& d1, A2_type&& A2, d2_type&& d2, A3_type&& A3, d3_type&& d3,
                                     ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K3K3operatorAABBCC_factored(std::forward<t_vector_type>(tvec),
                                                       std::forward<A1_type>(A1), std::forward<d1_type>(d1),
                                                       std::forward<A2_type>(A2), std::forward<d2_type>(d2),
                                                       std::forward<A3_type>(A3), std::forward<d3_type>(d3),
                                                       std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
        template <class t_vector_type, class A1_type, class d1_type, class A2_type, class d2_type, class A3_type, class d3_type, class... ParameterTypes>
    auto K3K3operatorABCABC_factored(t_vector_type&& tvec, A1_type&& A1, d1_type&& d1, A2_type&& A2, d2_type&& d2, A3_type&& A3, d3_type&& d3,
                                     ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->K3K3operatorABCABC_factored(std::forward<t_vector_type>(tvec),
                                                       std::forward<A1_type>(A1), std::forward<d1_type>(d1),
                                                       std::forward<A2_type>(A2), std::forward<d2_type>(d2),
                                                       std::forward<A3_type>(A3), std::forward<d3_type>(d3),
                                                       std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
    template <class t_vector_type, class... ParameterTypes>
    auto ineq_constraint(t_vector_type&& tvec, ParameterTypes&&... params) const {
        auto callback = [&]<typename... NewParamTypes>(NewParamTypes&&... new_params){
            return base_cgf()->ineq_constraint(std::forward<t_vector_type>(tvec), std::forward<NewParamTypes>(new_params)...);};
        return cfo()->operator()(callback, std::forward<ParameterTypes>(params)...);
    }
};

template <class...> class SPCFO_impl;

template <std::size_t... I, class... FunctionObjects>
class SPCFO_impl<std::index_sequence<I...>, FunctionObjects...> : private FunctionObjectWrapper<FunctionObjects, I>... {
public:
    SPCFO_impl() = default;

    template <class... FOSingleArgTypes>
    explicit SPCFO_impl(FOSingleArgTypes&&... args)
        : FunctionObjectWrapper<FunctionObjects, I>(std::forward<FOSingleArgTypes>(args))... {}
    // Constructor: one argument per function object

    template <class CallbackObject, class... ParameterTypes>
    auto operator()(const CallbackObject& co, const ParameterTypes&... params) const {
        return co(FunctionObjectWrapper<FunctionObjects, I>::function_object()->operator()(params...)...);
    }
};

template <class... FunctionObjects>
using SharedParametersCallbackFunctionObject = SPCFO_impl<std::index_sequence_for<FunctionObjects...>, FunctionObjects...>;

template <class BaseCGF, class... FunctionObjects>
using ReparametrisedCGF = AdaptedParametersCallbackCGF<BaseCGF, SharedParametersCallbackFunctionObject<FunctionObjects...>>;
// BaseCGF is a CGF class expecting k parameters of types T1, ..., Tk
// FunctionObjects are k function object classes, and the i^th class implements a method of the form
    // Ti operator()(ParamTypes... params) const;
// accepting any desired number of input parameters and returning a value of the type T that BaseCGF expects
// ReparametrisedCGF provides methods of the form K(tvec, params...) and K3operator(tvec, params..., v1, v2, v3), etc
// Note that the same input values params... are provided to each of the k function objects

// A common case is that each of the classes in FunctionObjects accepts a single shared vector argument and returns one vector result.
// See for instance the Adaptor class.
// This can be interpreted as a parametric submodel of the original distribution specified by BaseCGF, in which the vector of model parameters
// is transformed by FunctionObjects into one or more distributional parameters suppled to BaseCGF



template <class FunctionObject>
class ExtractPairTupleCallbackFunctionObject : private FunctionObject {
    // Provides a CallbackFunctionObject type to pass to AdaptedParametersCallbackCGF
    // FunctionObject::operator() accepts whatever appropriate number of parameters and returns a std::pair or std::tuple
    // ExtractPairTupleCallbackFunctionObject extracts the entries from the pair or tuple and passes these on as specified by AdaptedParametersCallbackCGF
    // This provides a simple way for a FunctionObject to return more than one parameter to be passed to an underlying CGF
public:
    using FunctionObject::FunctionObject;

    template <class CallbackObject, class ParameterTypes>
    auto operator()(CallbackObject&& co, ParameterTypes&& params) const {
        return std::apply(std::forward<CallbackObject>(co), static_cast<const FunctionObject&>(*this)(std::forward<ParameterTypes>(params)));
    }
};

template <class BaseCGF, class FunctionObject>
using CGFwithExtractor = AdaptedParametersCallbackCGF<BaseCGF, ExtractPairTupleCallbackFunctionObject<FunctionObject>>;
// BaseCGF is a CGF class expecting multiple parameters
// FunctionObject implements a method of the form
    // operator()(ParamTypes... params) const;
// accepting any desired number of input parameters and returning a std::tuple or std::pair
// CGFwithExtractor provides methods of the form K(tvec, params...) and K3operator(tvec, params..., v1, v2, v3), etc
// in which the entries from the pair or tuple FunctionObject(params...) are passed to BaseCGF

// A common case is that BaseCGF expects several distinct parameters (e.g., the shape and rate parameters for a Gamma distribution)
// which are to be extracted as particular entries from a single parameter vector.
// Then FunctionObject() should accept one vector argument, extract the corresponding entries, and store them in a pair or tuple.
// The details (for instance, in which order the parameters should be stored in the vector) can be hidden behind FunctionObject.
// This is useful for using template-based CGFs to implement the virtual-function-based interface of CGF_with_AD,
// in which all parameters must be passed in a single vector argument.

class IdentityCallbackFunctionObject {
public:
    template <class CallbackObject, class... ParameterTypes>
    auto operator()(const CallbackObject& co, ParameterTypes&&... params) const {
        return co(std::forward<ParameterTypes>(params)...);
    }
};
template <class BaseCGF>
using WrapAsCGF = AdaptedParametersCallbackCGF<BaseCGF, IdentityCallbackFunctionObject>;
// For use when BaseCGF is a pointer or reference type, such as CGF_with_AD*
// WrapAsCGF<CGF_as_pointer*> can be used as a template-style CGF object,
// for instance in order to be stored in a container (see SumOfIndependentCGF)
// (This mechanism also allows WrapAsCGF<CGF_as_pointer*> to be used as a base class directly,
// but current code primarily uses the BaseWrapper<> mechanism instead
// since this provides a simpler user interface for specifying template arguments;
// here, for instance, AdaptedParametersCallbackCGF already provides this wrapping.)


} // namespace CGFs_via_templates


namespace CGFs_with_AD {

class Adaptor {
// Abstract base class for an adaptor that accepts a model parameter vector
// and transforms it into distributional parameters to be passed to an underlying CGF
public:
    virtual ~Adaptor() = default;

    virtual vec operator()(const vec& model_parameter_vector) const = 0;
    virtual a_vector operator()(const a_vector& model_parameter_vector) const = 0;
    // Each of these operators returns the distributional parameters as a vector of the same
    // type as the input vector
};

using VectorAdaptor = Adaptor;
// By default, the output from an adaptor is stored in a single vector (possibly of size 1).
// Analogous classes can be used when the output is required to be a scalar or a matrix:
class ScalarAdaptor {
public:
    virtual ~ScalarAdaptor() = default;

    virtual double operator()(const vec& model_parameter_vector) const = 0;
    virtual a_scalar operator()(const a_vector& model_parameter_vector) const = 0;
    // Each of these operators returns the distributional parameters as a scalar of the
    // type matching the input vector
};
class MatrixAdaptor {
public:
    virtual ~MatrixAdaptor() = default;

    virtual mat operator()(const vec& model_parameter_vector) const = 0;
    virtual a_matrix operator()(const a_vector& model_parameter_vector) const = 0;
    // Each of these operators returns the distributional parameters as a matrix of the
    // type matching the input vector
};
// Note: currently MatrixAdaptor is not used
// To date, the only CGF that expects a matrix parameter is linearly_mappedCGF, and this matrix has always been assumed to be fixed

// The parametric_submodelCGF class implements the CGF for models where the parameter
    // vector is transformed by an adaptor, then passed to an underlying CGF.
// The adaptor is stored as a pointer to an abstract
    // base class of type Adaptor
// Functionality is obtained by implementing specific derived classes from Adaptor
// Two such derived classes are implemented below
              // Avoiding a direct alias here to avoid issues with forward declaration of parametric_submodelCGF class.
using parametric_submodelCGF_temp = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<const CGF_with_AD*, const Adaptor*>>;
class parametric_submodelCGF : public parametric_submodelCGF_temp{
public:
  // Forwarding constructor
  parametric_submodelCGF(const CGF_with_AD* base_cgf, const Adaptor* adaptor) : parametric_submodelCGF_temp(base_cgf, adaptor) {}
};

              // // Creator function allocates using new
              // parametric_submodelCGF* make_parametric_submodelCGF_with_AD(const CGF_with_AD* base_cgf,
              //                                                             const Adaptor* a) {
              //     return new parametric_submodelCGF(base_cgf, a);
              // }


// For specific named distributions implemented as template-style CGFs,
// the model parameter vector may be transformed into distributional parameters
// by supplying one adaptor for each parameter supplied to the distribution
// The following class simplifies the syntax for specifying this:
template <class BaseCGF, class... AdaptorTypes>
using SeparateParametersCGF = CGF_with_AD_from_template<saddlepoint::CGFs_via_templates::ReparametrisedCGF<BaseCGF, AdaptorTypes...>>;
// Sample usage:
// SeparateParametersCGF<saddlepoint::CGFs_via_templates::MultinomialCGF, const ScalarAdaptor*, const VectorAdaptor*>
// That is, two adaptors must be provided: one for N (supplying scalar values) and one for prob_vec (supplying vector values)
// Each adaptor is provided as a pointer-to-const of the appropriate type.



class AdaptorFromFunctionPair
: public Adaptor {
    public:
    // Convenience typedefs for pointers to functions
    typedef vec (*d_func_pointer)(const vec&);
    typedef a_vector (*ad_func_pointer)(const a_vector&);

    private:
    d_func_pointer f_d;
    ad_func_pointer f_ad;

    public:
    AdaptorFromFunctionPair(d_func_pointer f_double, ad_func_pointer f_ad_double)
    : f_d(f_double), f_ad(f_ad_double) {}
    // The user supplies two functions that transform the parameter vector
    // One function uses double as its scalar type, the other AD<double>

    vec operator()(const vec& model_parameter_vector) const override {
        return (*f_d)(model_parameter_vector);
    }
    a_vector operator()(const a_vector& model_parameter_vector) const override {
        return (*f_ad)(model_parameter_vector);
    }
};

class SubvectorAdaptor : public Adaptor {
private:
    Eigen::Index pos;
    Eigen::Index len;
public:
    SubvectorAdaptor(Eigen::Index position, Eigen::Index length) : pos(position), len(length) {}

    vec operator()(const vec& model_parameter_vector) const override {
        return model_parameter_vector.segment(pos, len);
    }
    a_vector operator()(const a_vector& model_parameter_vector) const override {
        return model_parameter_vector.segment(pos, len);
    }
};

class SubsetVectorByIndicesAdaptor : public Adaptor {
private:
    Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices;
public:
    explicit SubsetVectorByIndicesAdaptor(Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> indices_) : indices(indices_) {
        if (indices.minCoeff() < 0) throw std::out_of_range("indices cannot be negative - SubsetVectorByIndicesAdaptor");
    }

    vec operator()(const vec& model_parameter_vector) const override {
        // Check that all indices are within range
        if (indices.maxCoeff() >= model_parameter_vector.size()) {
            throw std::out_of_range("indices out of range in SubsetVectorByIndicesAdaptor.");
        }
        vec subset_vector(indices.size());
        for(decltype(indices.size()) i = 0; i < indices.size(); i++) {
                subset_vector(i) = model_parameter_vector(indices(i));
        }
        return subset_vector;
    }
    a_vector operator()(const a_vector& model_parameter_vector) const override {
        // Check that all indices are within range
        if (indices.maxCoeff() >= model_parameter_vector.size()) {
            throw std::out_of_range("indices out of range in SubsetVectorByIndicesAdaptor.");
        }
        a_vector a_subset_vector(indices.size());
        for(decltype(indices.size()) i = 0; i < indices.size(); i++) a_subset_vector(i) = model_parameter_vector(indices(i));
        return a_subset_vector;
    }
};

class SavedVectorAdaptor : public VectorAdaptor {
private:
    vec v;
public:
    explicit SavedVectorAdaptor(const vec& _v) : v(_v) {}

    SavedVectorAdaptor() = delete;
    // No default constructor as otherwise v will be initialised as an empty vector, which is undesirable

    vec operator()(const vec& parameter_vector) const override {return v;}
    // If parameter_vector is provided as type vec, return v as type vec
    a_vector operator()(const a_vector& parameter_vector) const override {return v.cast<a_scalar>();}
    // If parameter_vector is provided as type a_vector, convert v to equivalent type
};
class SavedScalarAdaptor : public ScalarAdaptor {
private:
    double x;
public:
    explicit SavedScalarAdaptor(double _x) : x(_x) {}

    SavedScalarAdaptor() = delete;
    // No default constructor as otherwise v will be initialised as an empty vector, which is undesirable

    double operator()(const vec& parameter_vector) const override {return x;}
    // If parameter_vector is provided as type vec, return v as type vec
    a_scalar operator()(const a_vector& parameter_vector) const override {return static_cast<a_scalar>(x);}
    // If parameter_vector is provided as type a_vector, convert v to equivalent type
};


class ScalarAdaptorFromVectorAdaptor : public ScalarAdaptor {
// This class takes a pointer to a VectorAdaptor as a constructor argument, and it checks that the size of the vector is 1 before returning a scalar.
// It also throws an exception if the size of the vector is not 1.
// Example:
// VectorAdaptor* vec_adaptor_ptr = new SubsetVectorByIndicesAdaptor(...);
// ScalarAdaptor* scalar_adaptor_ptr = new ScalarAdaptorFromVectorAdaptor(vec_adaptor_ptr);
private:
    VectorAdaptor* vec_adaptor_ptr;

public:
    explicit ScalarAdaptorFromVectorAdaptor(VectorAdaptor* _vec_adaptor_ptr) : vec_adaptor_ptr(_vec_adaptor_ptr) {
        if (vec_adaptor_ptr == nullptr)  throw std::invalid_argument("VectorAdaptor is null");
    }

    double operator()(const vec& parameter_vector) const override {
        const auto vector_size = (*vec_adaptor_ptr)(parameter_vector).size();
        // if ((*vec_adaptor_ptr)(parameter_vector).size() != 1)  throw std::invalid_argument("Vector size in the VectorAdaptor is not equal to 1");
         if (vector_size != 1) throw std::invalid_argument("Scalar parameter expected, but received a vector of size " + std::to_string(vector_size));
        return (*vec_adaptor_ptr)(parameter_vector)[0];
    }
    a_scalar operator()(const a_vector& parameter_vector) const override {
        const auto vector_size = (*vec_adaptor_ptr)(parameter_vector).size();
        // if ((*vec_adaptor_ptr)(parameter_vector).size() != 1)  throw std::invalid_argument("Vector size in the VectorAdaptor is not equal to 1");
        if (vector_size != 1) throw std::invalid_argument("Scalar parameter expected, but received a vector of size " + std::to_string(vector_size));
        return (*vec_adaptor_ptr)(parameter_vector)[0];
    }
};





} // namespace CGFs_with_AD
} // namespace saddlepoint


#endif // PARAMETRIC_SUBMODELCGF_H_INCLUDED
