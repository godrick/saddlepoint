#ifndef SCRATCHPAD_H_INCLUDED
#define SCRATCHPAD_H_INCLUDED

namespace saddlepoint {
namespace CGFs_via_virtual_functions {

template <class scalar_type, class vector_type, class matrix_type>
class CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad {
public:
    virtual ~Scratchpad() = default;
    
    virtual Scratchpad_Q* scratchpad_q(const matrix_type& Q) = 0;
    // Creates a new Scratchpad_Q object based on *this, with the square matrix Q added
    
    // Each CGF_base function has a Scratchpad equivalent with the arguments tvec, theta omitted
    virtual const scalar_type& K() = 0;
    virtual const vector_type& K1() = 0;
    virtual const matrix_type& K2() = 0;
    
    virtual const scalar_type& tilting_exponent() = 0;
    virtual const scalar_type& neg_ll() = 0;
    virtual const scalar_type& func_T() = 0;
    // No default implementations are provided; for instance, tilting_exponent() requires the value of tvec, which might not be available
    
    // Note: The above methods are non-const to accommodate the fact that running a method may update the saved internal state.
    // In principle this could be handled using "mutable" but for this application the reuse of previously saved values is close enough
        // to being part of the visible state of the object.
        
    class BasicImplementation;
};

template <class scalar_type, class vector_type, class matrix_type>
typename CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad* CGF_base<scalar_type, vector_type, matrix_type>::scratchpad(const vector_type& tvec, const vector_type& parameter_vector) const {
    return new typename Scratchpad::BasicImplementation(tvec, parameter_vector, this);
}
    
//-----------------------------------------------------------
template <class scalar_type, class vector_type, class matrix_type>
class CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q : public Scratchpad {
    // A similar class where in addition a square matrix Q is fixed
    // This matrix is assigned as the value for all of the arguments Q1, Q2, Q3 in all the K{n}operator{...} methods
public:
    virtual const scalar_type& K4operatorAABB() = 0;
    virtual const scalar_type& K3K3operatorAABBCC() = 0;
    virtual const scalar_type& K3K3operatorABCABC() = 0;
    
    class BasicImplementation;
};

//-----------------------------------------------------------
// Implementation of basic scratchpad classes
//-----------------------------------------------------------
template <class result_type, class FunctionObject>
class SavedResult : private FunctionObject {
    // A convenience object that saves the result of a calculation.
    // When evaluated for the first time, the result is computed using FunctionObject, then saved.
    // On subsequent evaluations, the saved result is returned as a const reference.
    
    // FunctionObject must provide an operator() taking any desired arguments and returning a value that can be assigned to result_type
    // result_type must be DefaultConstructible and Assignable
private:
    result_type val;
    bool computed;
public:
    template <class... ArgTypes>
    const result_type& operator()(ArgTypes&&... args) {
        if (!computed) {
            val = FunctionObject::operator()(std::forward<ArgTypes>(args)...);
        }
        return val;
    }
    
    template <class... ArgTypes>
    SavedResult(ArgTypes&&... args) : FunctionObject(std::forward<ArgTypes>(args)...), computed(false) {} // val is default-constructed
    
    SavedResult() : computed(false) {}
};
//-----------------------------------------------------------
template <class scalar_type, class vector_type, class matrix_type>
class CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad::BasicImplementation : public Scratchpad {
    // The basic implementation stores a copy of the supplied values and a pointer to the CGF object.
    // All methods calculate their values the first time by forwarding these stored values to the CGF object, when needed.
    // Each computed value is saved once computed, and subsequent calls return the saved values.
    
    // Note: this implementation assumes that scalar_type, vector_type, matrix_type are DefaultConstructible and Assignable.
protected:
    const CGF_base* p_cgf;
    vector_type tvec;
    vector_type theta;
    
    struct K_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->K(t, th);}
    };
    SavedResult<scalar_type, K_finder> saved_k;
    
    struct K1_finder {
        vector_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->K1(t, th);}
    };
    SavedResult<vector_type, K1_finder> saved_k1;
    
    struct K2_finder {
        matrix_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->K2(t, th);}
    };
    SavedResult<matrix_type, K2_finder> saved_k2;

    struct tilting_exponent_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->tilting_exponent(t, th);}
    };
    SavedResult<scalar_type, tilting_exponent_finder> saved_te;
    
    struct neg_ll_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->neg_ll(t, th);}
    };
    SavedResult<scalar_type, neg_ll_finder> saved_nll;
    
    struct func_T_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const CGF_base* p) {return p->func_T(t, th);}
    };
    SavedResult<scalar_type, func_T_finder> saved_ft;
    
    friend class CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q::BasicImplementation;
    
public:
    BasicImplementation(const vector_type& t, const vector_type& th, const CGF_base* p) 
        : p_cgf(p), tvec(t), theta(th) {} 
          // saved values are default constructed
    
    Scratchpad_Q* scratchpad_q(const matrix_type& Q) override {
        return new typename CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q::BasicImplementation(Q, this);
    }
    
    // Remaining methods pass values to SavedResult<...> objects, 
        // which check whether the value has been saved already, compute and save it if needed, and return the saved result
    const scalar_type& K() override {return saved_k(tvec, theta, p_cgf);}
    const vector_type& K1() override {return saved_k1(tvec, theta, p_cgf);}
    const matrix_type& K2() override {return saved_k2(tvec, theta, p_cgf);}

    const scalar_type& tilting_exponent() override {return saved_te(tvec, theta, p_cgf);}
    const scalar_type& neg_ll() override {return saved_nll(tvec, theta, p_cgf);}
    const scalar_type& func_T() override {return saved_ft(tvec, theta, p_cgf);}
};
//-----------------------------------------------------------
template <class scalar_type, class vector_type, class matrix_type>
class CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad_Q::BasicImplementation : public Scratchpad_Q {
protected:
    typedef typename CGF_base<scalar_type, vector_type, matrix_type>::Scratchpad::BasicImplementation Basic_SP_type;
    
    matrix_type Q;
    Basic_SP_type* plain_sp;
    
    struct K4operatorAABB_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const matrix_type& q, const CGF_base* p) {
            return p->K4operatorAABB(t, th, q, q);
            }
    };
    SavedResult<scalar_type, K4operatorAABB_finder> saved_k4aabb;
    
    struct K3K3operatorAABBCC_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const matrix_type& q, const CGF_base* p) {
            return p->K3K3operatorAABBCC(t, th, q, q, q);
            }
    };
    SavedResult<scalar_type, K3K3operatorAABBCC_finder> saved_k3k3aabbcc;
    
    struct K3K3operatorABCABC_finder {
        scalar_type operator()(const vector_type& t, const vector_type& th, const matrix_type& q, const CGF_base* p) {
            return p->K3K3operatorABCABC(t, th, q, q, q);
            }
    };
    SavedResult<scalar_type, K3K3operatorABCABC_finder> saved_k3k3abcabc;
    
public:
    BasicImplementation(const matrix_type& q, Basic_SP_type* psp)
        : Q(q), plain_sp(psp) {}
        // saved values are default constructed
    
    Scratchpad_Q* scratchpad_q(const matrix_type& q) override {
        // Note: method creates a new Scratchpad_Q object with a new value q, in addition to the current one
        return plain_sp->scratchpad_q(q);
    }
    
    const scalar_type& K() override {return plain_sp->K();}
    const vector_type& K1() override {return plain_sp->K1();}
    const matrix_type& K2() override {return plain_sp->K2();}
    const scalar_type& tilting_exponent() override {return plain_sp->tilting_exponent();}
    const scalar_type& neg_ll() override {return plain_sp->neg_ll();}
    const scalar_type& func_T() override {return plain_sp->func_T();}
    // Already implemented by plain_sp
    
    const scalar_type& K4operatorAABB() override {
        return saved_k4aabb(plain_sp->tvec, plain_sp->theta, Q, plain_sp->p_cgf);
    }
    const scalar_type& K3K3operatorAABBCC() override {
        return saved_k3k3aabbcc(plain_sp->tvec, plain_sp->theta, Q, plain_sp->p_cgf);
    }
    const scalar_type& K3K3operatorABCABC() override {
        return saved_k3k3abcabc(plain_sp->tvec, plain_sp->theta, Q, plain_sp->p_cgf);
    }
};

} // namespace CGFs_via_virtual_functions
} // namespace saddlepoint

#endif // SCRATCHPAD_H_INCLUDED
