#ifndef SCALARCGF_H_INCLUDED
#define SCALARCGF_H_INCLUDED

#include <cmath>
#include "BaseWrapper.h"
#include "saddlepoint_types.h"

namespace saddlepoint {
namespace CGFs_via_virtual_functions {

// Defines an interface for specifying the CGF of a scalar distribution
// with corresponding functions to handle a vector of i.i.d. scalar entries
// or a vector of independent, but not identically distributed, entries

// This can be done either via templates or via virtual functions
// FOR NOW: via virtual functions, similar to CGF_base
// TO DO: add template version as needed - similar to CGF_with_AD_from_template, that creates instances from a template-style one
// and (possibly) also a version that accepts both plain and AD vectors
template <class vector_type>
// Unlike CGF_base, we do not return any matrices so only one template parameter
// It is assumed that vector_type supports an Eigen-style .array() method that allows for vectorised entrywise operations
class ScalarCGF_baseIID {
public:
    virtual ~ScalarCGF_baseIID() = default;

    // Distribution-specific derivatives, as vectorised functions
    virtual vector_type K_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K1_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K2_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K3_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K4_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    // Returns a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // parameter_vector contains the parameter(s) used across all calculations (as many as needed for the particular distribution)
    // Thus parameter_vector.size() is fixed no matter how large tvec may be


    // Optional functions, in vectorised form, that can be obtained from K, K1 and K2
    // For specific distributions this may be able to be implemented more efficiently, in which case this default function may be overridden
    virtual vector_type tilting_exponent_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns K-t*K1 in each entry
        return K_vectorised_iid(tvec, parameter_vector) - tvec.array() * K1_vectorised_iid(tvec, parameter_vector);
    }
    virtual vector_type neg_ll_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns -(K-t*K1 - 0.5*log(K2)) in each entry
        return 0.5*(2*M_PI*K2_vectorised_iid(tvec, parameter_vector).array()).log() - tilting_exponent_vectorised_iid(tvec, parameter_vector);
    }
    virtual vector_type func_T_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns K4/(8*K2*K2) - 5*K3*K3/(24*K2*K2*K2) in each entry
        vector_type k2val = K2_vectorised_iid(tvec, parameter_vector);
        vector_type k2sq_val = k2val.array() * k2val.array();
        vector_type k3val = K3_vectorised_iid(tvec, parameter_vector);
        vector_type k4val = K4_vectorised_iid(tvec, parameter_vector);
        return k4val.array()/(8*k2sq_val.array()) - 5*k3val.array()*k3val.array()/(24*k2sq_val.array()*k2val.array());
    }
    virtual vector_type ineq_constraint_vectorised_iid(const vector_type& tvec, const vector_type& parameter_vector) const {
        vector_type res;
        return res;
    }
};


template <class vector_type>
class ScalarCGF_baseNonIdentical {
public:
    virtual ~ScalarCGF_baseNonIdentical() = default;

    virtual vector_type K_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K1_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K2_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K3_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    virtual vector_type K4_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const = 0;
    // Returns a vector of the same size as tvec, with K or the appropriate derivative applied to each entry
    // Each entry of tvec has its own parameter entry/entries in parameter_vector
    // Thus parameter_vector.size() will be a multiple of tvec.size() (the multiple is the number of parameters for the particular distribution)
    // TO DO: Document the order of parameter_vector, when there are two or more parameters per entry of tvec
    // e.g. for Gamma: should all shape parameters precede all rate parameters, or should shape/rate pairs be repeated?


    virtual vector_type tilting_exponent_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns K-t*K1 in each entry
        return K_vectorisedNonIdentical(tvec, parameter_vector) - tvec.array() * K1_vectorisedNonIdentical(tvec, parameter_vector);
    }
    virtual vector_type neg_ll_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns -(K-t*K1 - 0.5*log(K2)) in each entry
        return 0.5*(2*M_PI*K2_vectorisedNonIdentical(tvec, parameter_vector).array()).log() - tilting_exponent_vectorisedNonIdentical(tvec, parameter_vector);
    }
    virtual vector_type func_T_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const {
        // Returns K4/(8*K2*K2) - 5*K3*K3/(24*K2*K2*K2) in each entry
        vector_type k2val = K2_vectorisedNonIdentical(tvec, parameter_vector);
        vector_type k2sq_val = (k2val.array() * k2val.array()).matrix();
        vector_type k3val = K3_vectorisedNonIdentical(tvec, parameter_vector);
        return K4_vectorisedNonIdentical(tvec, parameter_vector).array()/(8*k2sq_val.array()) - 5*k3val.array()*k3val.array()/(24*k2sq_val.array()*k2val.array());
    }
    virtual vector_type ineq_constraint_vectorisedNonIdentical(const vector_type& tvec, const vector_type& parameter_vector) const {
        vector_type res;
        return res;
    }
};

} // namespace CGFs_via_virtual_functions


namespace CGFs_via_templates {


template <class ScalarVectorisedIID>
class ScalarVectorisedIID_Defaults {
public:
    // Provides default behaviour for (1) tilting_exponent_vectorised_iid
                                   // (2) neg_ll_vectorised_iid
                                   // (3) func_T_vectorised_iid

    // ScalarVectorisedIID must provide methods K_vectorised_iid, K1_vectorised_iid, K2_vectorised_iid, K3_vectorised_iid, K4_vectorised_iid
    // Any of the three above methods not implemented while creating ScalarVectorisedIID will call the version in ScalarVectorisedIID_Defaults
    template <class t_vector_type, class... ParameterTypes>
    auto tilting_exponent_vectorised_iid(const t_vector_type& tvec, const ParameterTypes&... params) const {
        return ((static_cast<const ScalarVectorisedIID*>(this)->K_vectorised_iid(tvec, params...)).array()
               - ( tvec.array() * (static_cast<const ScalarVectorisedIID*>(this)->K1_vectorised_iid(tvec, params...)).array() )).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto neg_ll_vectorised_iid(const t_vector_type& tvec, const ParameterTypes&... params) const {
        return (0.5*(2*M_PI*static_cast<const ScalarVectorisedIID*>(this)->K2_vectorised_iid(tvec, params...).array()).log() - tilting_exponent_vectorised_iid(tvec, params...)).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto func_T_vectorised_iid(const t_vector_type& tvec, const ParameterTypes&... params) const {
        t_vector_type k2val = static_cast<const ScalarVectorisedIID*>(this)->K2_vectorised_iid(tvec, params...);
        t_vector_type k2sq_val = (k2val.array() * k2val.array()).matrix();
        t_vector_type k3val = static_cast<const ScalarVectorisedIID*>(this)->K3_vectorised_iid(tvec, params...);
        return (static_cast<const ScalarVectorisedIID*>(this)->K4_vectorised_iid(tvec, params...).array()/(8*k2sq_val.array()) - 5*k3val.array()*k3val.array()/(24*k2sq_val.array()*k2val.array())).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto ineq_constraint_vectorised_iid(const t_vector_type& tvec, const ParameterTypes&... params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t< decltype(tvec.eval()) > > vector_type;
        vector_type res;
        return res;
    }
}; // class ScalarVectorisedIID_Defaults
//----------------
//----------------
//----------------
template <class ScalarVectorisedNonIdentical>
class ScalarVectorisedNonIdentical_Defaults {
public:
    template <class t_vector_type, class... ParameterTypes>
    auto tilting_exponent_vectorisedNonIdentical(const t_vector_type& tvec, const ParameterTypes&... params) const {
        return (static_cast<const ScalarVectorisedNonIdentical*>(this)->K_vectorisedNonIdentical(tvec, params...)
               - ( tvec.array() * static_cast<const ScalarVectorisedNonIdentical*>(this)->K1_vectorisedNonIdentical(tvec, params...) )).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto neg_ll_vectorisedNonIdentical(const t_vector_type& tvec, const ParameterTypes&... params) const {
        return (0.5*(2*M_PI*static_cast<const ScalarVectorisedNonIdentical*>(this)->K2_vectorisedNonIdentical(tvec, params...).array()).log() - tilting_exponent_vectorisedNonIdentical(tvec, params...)).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto func_T_vectorisedNonIdentical(const t_vector_type& tvec, const ParameterTypes&... params) const {
        t_vector_type k2val = static_cast<const ScalarVectorisedNonIdentical*>(this)->K2_vectorisedNonIdentical(tvec, params...);
        t_vector_type k2sq_val = (k2val.array() * k2val.array()).matrix();
        t_vector_type k3val = static_cast<const ScalarVectorisedNonIdentical*>(this)->K3_vectorisedNonIdentical(tvec, params...);
        return (static_cast<const ScalarVectorisedNonIdentical*>(this)->K4_vectorisedNonIdentical(tvec, params...).array()/(8*k2sq_val.array()) - 5*k3val.array()*k3val.array()/(24*k2sq_val.array()*k2val.array())).eval();
    }
    template <class t_vector_type, class... ParameterTypes>
    auto ineq_constraint_vectorisedNonIdentical(const t_vector_type& tvec, const ParameterTypes&... params) const {
        typedef typename std::remove_cv_t<typename std::remove_reference_t< decltype(tvec.eval()) > > vector_type;
        vector_type res;
        return res;
    }
}; // class ScalarVectorisedNonIdentical_Defaults


}

} // namespace saddlepoint

namespace saddlepoint {
namespace CGFs_with_AD {

class ScalarCGF_baseIID_with_AD : public CGFs_via_virtual_functions::ScalarCGF_baseIID<vec>,
                                  public CGFs_via_virtual_functions::ScalarCGF_baseIID<a_vector> {
public:
    // To allow overload resolution to operate within ScalarCGF_base_with_AD,
    using ScalarCGF_baseIID<vec>::K_vectorised_iid;
    using ScalarCGF_baseIID<vec>::K1_vectorised_iid;
    using ScalarCGF_baseIID<vec>::K2_vectorised_iid;
    using ScalarCGF_baseIID<vec>::K3_vectorised_iid;
    using ScalarCGF_baseIID<vec>::K4_vectorised_iid;
    using ScalarCGF_baseIID<vec>::tilting_exponent_vectorised_iid;
    using ScalarCGF_baseIID<vec>::neg_ll_vectorised_iid;
    using ScalarCGF_baseIID<vec>::func_T_vectorised_iid;
    using ScalarCGF_baseIID<vec>::ineq_constraint_vectorised_iid;
//----------------------------
    using ScalarCGF_baseIID<a_vector>::K_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::K1_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::K2_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::K3_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::K4_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::tilting_exponent_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::neg_ll_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::func_T_vectorised_iid;
    using ScalarCGF_baseIID<a_vector>::ineq_constraint_vectorised_iid;

};


class ScalarCGF_baseNonIdentical_with_AD : public CGFs_via_virtual_functions::ScalarCGF_baseNonIdentical<vec>,
                                           public CGFs_via_virtual_functions::ScalarCGF_baseNonIdentical<a_vector> {
public:
    using ScalarCGF_baseNonIdentical<vec>::K_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::K1_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::K2_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::K3_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::K4_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::tilting_exponent_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::neg_ll_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::func_T_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<vec>::ineq_constraint_vectorisedNonIdentical;
//----------------------------
    using ScalarCGF_baseNonIdentical<a_vector>::K_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::K1_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::K2_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::K3_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::K4_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::tilting_exponent_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::neg_ll_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::func_T_vectorisedNonIdentical;
    using ScalarCGF_baseNonIdentical<a_vector>::ineq_constraint_vectorisedNonIdentical;

};


template <class templateCGF>
class ScalarCGF_IID_with_AD_from_template : protected templateCGF, public ScalarCGF_baseIID_with_AD {
//private:
//    using BaseWrapper<templateCGF>::base_cgf;
public:
    using templateCGF::templateCGF;

    vec K_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K_vectorised_iid(tvec, parameter_vector);}
    a_vector K_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K_vectorised_iid(tvec, parameter_vector);}

    vec K1_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K1_vectorised_iid(tvec, parameter_vector);}
    a_vector K1_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K1_vectorised_iid(tvec, parameter_vector);}

    vec K2_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K2_vectorised_iid(tvec, parameter_vector);}
    a_vector K2_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K2_vectorised_iid(tvec, parameter_vector);}

    vec K3_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K3_vectorised_iid(tvec, parameter_vector);}
    a_vector K3_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K3_vectorised_iid(tvec, parameter_vector);}

    vec K4_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K4_vectorised_iid(tvec, parameter_vector);}
    a_vector K4_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K4_vectorised_iid(tvec, parameter_vector);}

    vec tilting_exponent_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::tilting_exponent_vectorised_iid(tvec, parameter_vector);}
    a_vector tilting_exponent_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::tilting_exponent_vectorised_iid(tvec, parameter_vector);}

    vec neg_ll_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::neg_ll_vectorised_iid(tvec, parameter_vector);}
    a_vector neg_ll_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::neg_ll_vectorised_iid(tvec, parameter_vector);}

    vec func_T_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::func_T_vectorised_iid(tvec, parameter_vector);}
    a_vector func_T_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::func_T_vectorised_iid(tvec, parameter_vector);}

    vec ineq_constraint_vectorised_iid(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::ineq_constraint_vectorised_iid(tvec, parameter_vector);}
    a_vector ineq_constraint_vectorised_iid(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::ineq_constraint_vectorised_iid(tvec, parameter_vector);}
};



template <class templateCGF>
class ScalarCGF_NonIdenticalwith_AD_from_template : public ScalarCGF_baseNonIdentical_with_AD, protected templateCGF{
//private:
//    using BaseWrapper<templateCGF>::base_cgf;
public:
    using templateCGF::templateCGF;

    vec K_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector K_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K_vectorisedNonIdentical(tvec, parameter_vector);}

    vec K1_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K1_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector K1_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K1_vectorisedNonIdentical(tvec, parameter_vector);}

    vec K2_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K2_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector K2_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K2_vectorisedNonIdentical(tvec, parameter_vector);}

    vec K3_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K3_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector K3_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K3_vectorisedNonIdentical(tvec, parameter_vector);}

    vec K4_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::K4_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector K4_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::K4_vectorisedNonIdentical(tvec, parameter_vector);}

    vec tilting_exponent_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::tilting_exponent_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector tilting_exponent_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::tilting_exponent_vectorisedNonIdentical(tvec, parameter_vector);}

    vec neg_ll_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::neg_ll_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector neg_ll_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::neg_ll_vectorisedNonIdentical(tvec, parameter_vector);}

    vec func_T_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::func_T_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector func_T_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::func_T_vectorisedNonIdentical(tvec, parameter_vector);}

    vec ineq_constraint_vectorisedNonIdentical(const vec& tvec, const vec& parameter_vector) const override {return templateCGF::ineq_constraint_vectorisedNonIdentical(tvec, parameter_vector);}
    a_vector ineq_constraint_vectorisedNonIdentical(const a_vector& tvec, const a_vector& parameter_vector) const override {return templateCGF::ineq_constraint_vectorisedNonIdentical(tvec, parameter_vector);}

};



}
}

#endif // SCALARCGF_H_INCLUDED
