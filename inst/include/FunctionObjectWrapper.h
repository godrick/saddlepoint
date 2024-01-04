#ifndef FUNCTIONOBJECTWRAPPER_H_INCLUDED
#define FUNCTIONOBJECTWRAPPER_H_INCLUDED

namespace saddlepoint {
namespace CGFs_via_templates {

// Utility wrapper for function objects
// Typically these CGFs are empty classes (i.e., have no data members)
// and can be incorporated as (private or protected) base classes.
// However, pointer or reference types cannot be base classes,
// and this is inconvenient when using the CGF_with_AD-style Adaptor class.
// To provide a consistent interface, wrap all FunctionObjects via FunctionObjectWrapper,
// which provides access to operator() methods, 
// plus a named accessor for the underlying FunctionObject as a pointer.
template <class FunctionObject, int n = 0>
// Optional parameter n is not used but allows creating multiple wrapped FunctionObjects that would be otherwise the same
class FunctionObjectWrapper : protected FunctionObject {
public:
    using FunctionObject::FunctionObject;
    FunctionObjectWrapper(const FunctionObject& fo) : FunctionObject(fo) {}
    using FunctionObject::operator();
    const FunctionObject* function_object() const {return this;}
};

template <class FObyPointer, int n>
class FunctionObjectWrapper<FObyPointer*, n> {
protected:
    FObyPointer* p_fo;
public:
    explicit FunctionObjectWrapper(FObyPointer* p) : p_fo(p) {}
    
    FunctionObjectWrapper() = delete;
    
    FObyPointer* function_object() {return p_fo;}
    const FObyPointer* function_object() const {return p_fo;}
    
    template <class... ArgTypes>
    auto operator()(ArgTypes&&... args) {return function_object()->operator()(std::forward<ArgTypes>(args)...);}
    template <class... ArgTypes>
    auto operator()(ArgTypes&&... args) const {return function_object()->operator()(std::forward<ArgTypes>(args)...);}    
};

template <class FObyPointer, int n>
class FunctionObjectWrapper<const FObyPointer*, n> {
protected:
    const FObyPointer* p_fo;
public:
    explicit FunctionObjectWrapper(const FObyPointer* p) : p_fo(p) {}
    
    FunctionObjectWrapper() = delete;
    
    const FObyPointer* function_object() const {return p_fo;}
    
    template <class... ArgTypes>
    auto operator()(ArgTypes&&... args) const {return function_object()->operator()(std::forward<ArgTypes>(args)...);}    
};

template <class FObyReference, int n>
class FunctionObjectWrapper<FObyReference&, n> : protected FunctionObjectWrapper<FObyReference*> {
public:
    explicit FunctionObjectWrapper(FObyReference& r) : FunctionObjectWrapper<FObyReference*>(&r) {}
    using FunctionObjectWrapper<FObyReference*>::operator();
    using FunctionObjectWrapper<FObyReference*>::function_object;
};

template <class FObyReference, int n>
class FunctionObjectWrapper<const FObyReference&, n> : protected FunctionObjectWrapper<const FObyReference*> {
public:
    explicit FunctionObjectWrapper(const FObyReference& r) : FunctionObjectWrapper<const FObyReference*>(&r) {}
    using FunctionObjectWrapper<const FObyReference*>::operator();
    using FunctionObjectWrapper<const FObyReference*>::function_object;
};

} // namespace CGFs_via_templates
} // namespace saddlepoint

#endif // FUNCTIONOBJECTWRAPPER_H_INCLUDED
