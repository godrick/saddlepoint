#ifndef SUPPLIEDFUNCTIONS_H_INCLUDED
#define SUPPLIEDFUNCTIONS_H_INCLUDED

#include <Rcpp.h>
#include <Eigen/Eigen>

namespace saddlepoint {
namespace CGFs_with_Rcpp {

class SuppliedFunctions
{
private:
    Rcpp::Function fn0;
    Rcpp::Function grad_fn0;
public:
    SuppliedFunctions(Rcpp::Function fn1, Rcpp::Function grad_fn1) : fn0(fn1), grad_fn0(grad_fn1) { }
//----------------------------------------------------------
    Eigen::VectorXd value(Eigen::VectorXd theta) {
        Eigen::VectorXd val;
        // TO DO: Check if this is necessary. I think R should be able to handle this kind of error
        try { val = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(fn0(theta)); }
        catch (const std::exception& e) { throw std::invalid_argument("The provided argument is invalid for the supplied R function."); }
        return val;
    }
    // Eigen::VectorXd value(Eigen::VectorXd theta) {
    //    Eigen::Map<Eigen::VectorXd> val = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(fn0(theta));
    //    return val;
    // }
//----------------------------------------------------------
    Eigen::MatrixXd gradient(Eigen::VectorXd theta) {
        Eigen::VectorXd fn0_val = value(theta);
        Eigen::MatrixXd grad_val;
        // TO DO: Check if this is necessary. I think R should be able to handle this kind of error
        try { grad_val = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(grad_fn0(theta)); }
        catch (const std::exception& e) { throw std::invalid_argument("The provided argument is invalid for the supplied gradient function."); }

        if (grad_val.rows() != fn0_val.size() || grad_val.cols() != theta.size())
            throw std::invalid_argument("Invalid gradient dimension: expected dimension " + std::to_string(fn0_val.size()) + " by " + std::to_string(theta.size()));

        return grad_val;
    }
    // Eigen::MatrixXd gradient(Eigen::VectorXd theta) {
    //    Eigen::Map<Eigen::MatrixXd> grad_val = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(grad_fn0(theta));
    //    return grad_val;
    // }
//----------------------------------------------------------
};

} // namespace CGFs_with_Rcpp
} // namespace saddlepoint

#endif // SUPPLIEDFUNCTIONS_H_INCLUDED
