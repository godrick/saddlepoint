#ifndef GRADFUNCT_INCLUDED
#define GRADFUNCT_INCLUDED

#include "saddlepoint_types.h"
#include "atomic_funcs.h"
#include "parametric_submodelCGF.h"

namespace saddlepoint{

Rcpp::List grad_theta_funcT(vec theta, vec tvec,
                            mat negQB, Rcpp::XPtr<parametric_submodelCGF> modelCGF)
{
    a_vector a_theta = theta.cast<CppAD::AD<double>>();
    CppAD::Independent(a_theta);

    atomic_funcs::our_atomic_class_for_supplied_values tvec_afun("tvec_atomic", tvec, negQB);
    a_vector a_tvec(tvec.rows());
    tvec_afun(a_theta, a_tvec);

    a_vector a_funcT(1);
    a_funcT(0) = modelCGF->func_T(a_tvec, a_theta);

    CppAD::ADFun<double>* ADFun_ptr = new CppAD::ADFun<double>;
    ADFun_ptr->Dependent(a_theta, a_funcT);

    return Rcpp::List::create(Rcpp::Named("del.theta.t_hat") = negQB,
                              Rcpp::Named("func.T") = ADFun_ptr->Forward(0, theta),
                              Rcpp::Named("grad.theta.T") = ADFun_ptr->Jacobian(theta));
}

}// namespace saddlepoint
#endif // GRADFUNCT_INCLUDED
