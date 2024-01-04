#ifndef BASEWRAPPER_H_INCLUDED
#define BASEWRAPPER_H_INCLUDED

namespace saddlepoint {
namespace CGFs_via_templates {

// Utility wrapper for "base" CGFs implemented as template classes.
// Typically these CGFs are empty classes (i.e., have no data members)
// and are therefore incorporated as (private or protected) base classes.
// However, pointer or reference types cannot be base classes,
// and this is inconvenient when using CGF_with_AD and other pointer-based CGF classes.
// To provide a consistent interface, inherit template CGF arguments wrapped via BaseWrapper,
// and access via base_cgf() and base_cgf()->base_class_method()
template <class BaseCGF, int n = 0> 
// Optional parameter n is not used but allows a CGF to have multiple base CGFs that may be otherwise the same
class BaseWrapper : protected BaseCGF {
protected:
    using BaseCGF::BaseCGF;
    BaseWrapper(const BaseCGF& bcgf) : BaseCGF(bcgf) {}
    BaseCGF* base_cgf() {return this;}
    const BaseCGF* base_cgf() const {return this;}
};

// If BaseCGF is a pointer or reference (to const), it cannot be a base class and can be stored as a data member.
// Note that BaseWrapper provides no memory management, and the user is responsible for ensuring that
// the objects underlying the pointer or reference remains valid during the whole lifetime of the class using it.
template <class BaseCGFasPointer, int n>
class BaseWrapper<BaseCGFasPointer*, n> {
protected:
    BaseCGFasPointer* p_base;
    explicit BaseWrapper(BaseCGFasPointer* p) : p_base(p) {} 
    BaseCGFasPointer* base_cgf() {return p_base;}
    const BaseCGFasPointer* base_cgf() const {return p_base;}
    BaseWrapper() = delete;
};
template <class BaseCGFasPointer, int n>
class BaseWrapper<const BaseCGFasPointer*, n> {
protected:
    const BaseCGFasPointer* p_base;
    explicit BaseWrapper(const BaseCGFasPointer* p) : p_base(p) {} 
    const BaseCGFasPointer* base_cgf() const {return p_base;} // only const version
    BaseWrapper() = delete;
};

template <class BaseCGFasReference, int n>
class BaseWrapper<BaseCGFasReference&, n> : private BaseWrapper<BaseCGFasReference*> {
protected:
    explicit BaseWrapper(BaseCGFasReference& r) : BaseWrapper<BaseCGFasReference*>(&r) {}
    using BaseWrapper<BaseCGFasReference*>::base_cgf;
};
template <class BaseCGFasReference, int n>
class BaseWrapper<const BaseCGFasReference&, n> : private BaseWrapper<const BaseCGFasReference*> {
protected:
    explicit BaseWrapper(const BaseCGFasReference& r) : BaseWrapper<const BaseCGFasReference*>(&r) {}
    using BaseWrapper<const BaseCGFasReference*>::base_cgf;
};

} // namespace CGFs_via_templates
} // namespace saddlepoint

#endif // BASEWRAPPER_H_INCLUDED
