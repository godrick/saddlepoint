#ifndef ATOMIC_FUNCS_H_INCLUDED
#define ATOMIC_FUNCS_H_INCLUDED

#include <cppad/cppad.hpp>
#include <Eigen/Eigen>

#include "saddlepoint_types.h"


namespace saddlepoint {
namespace atomic_funcs {

//mat2vec function converts a matrix to a AD-vector representation
CppAD::vector<a_scalar> mat2vec(mat x){
  int n=x.size();
  CppAD::vector<a_scalar> res(n);
  for(int i=0;i<n;i++)res[i]=x(i);
  return res;
}

// atomic_logdet class defines an atomic function for the derivative of the log-determinant of a matrix,
    // which can be used in CppAD for automatic differentiation. The function implements the forward mode
    // for up to first order derivatives and the reverse mode for first order derivatives.
// We computing the log-determinant of the matrix by LU decomposition
// This has an adavntage because the LU decomposition of a matrix factorizes the matrix into a lower and upper triangular matrices,
   // then the desired determinant is simply the product of the diagonal elements of the two triangular matrices.
   // Then the logarithm of the determinant is simply the summation the logarithms of the diagonal elements which is computatiionally efficiently.
// Computation of this by LU decomposition is less prone to roundoff errors, making this a stable way in comparison to computation of the determinant directly.
   // This is especially necessary for matrices with large or small elements.

// This class uses some codes from TMB
class atomic_logdet : public CppAD::atomic_three<double> {
public:
     // constructor
     atomic_logdet(const std::string& name) : CppAD::atomic_three<double>(name) { }
private:
    bool for_type(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        CppAD::vector<CppAD::ad_type_enum>&        type_y      ) override
    {
        type_y[0] = CppAD::ad_type_enum::variable_enum ;
        return true;
    }
    bool forward
    (
        const CppAD::vector<double>&              parameter_x  ,
        const CppAD::vector<CppAD::ad_type_enum>& type_x       ,
        size_t                             need_y       ,
        size_t                             order_low    ,
        size_t                             order_up     ,
        const CppAD::vector<double>&              taylor_x     ,
        CppAD::vector<double>&                    taylor_y     ) override
    {
        size_t q1 = order_up + 1;
        size_t n_arg = parameter_x.size();
        size_t n = sqrt(n_arg);

        mat X(n, n);
        for(size_t i = 0; i < n_arg; i++) X.data()[i] = taylor_x[i * q1 + 0];
        mat invX = X.inverse();
        CppAD::vector<a_scalar> invX_vec = mat2vec(invX);

        bool ok = order_up <= 1;
        // Order zero forward mode
        if( order_low <= 0 ){
            //TMB
            mat LU=X.lu().matrixLU();    // Use Eigen LU decomposition
            vec LUdiag = LU.diagonal();
            double res=LUdiag.array().abs().log().sum();
            taylor_y[0] = res;
        }
        if( order_up <= 0 )
            return ok;
        // Order one forward mode.
        if( order_low <= 1 ){
            taylor_y[1]  = 0.0; //initialize sum
            for(size_t i = 0; i < n_arg; i++){
                taylor_y[1] += CppAD::Value(CppAD::Var2Par(invX_vec[i])) * taylor_x[i*q1 + 1];
            }
        }
        if( order_up <= 1 )
            return ok;
        // We only need forward mode upto order 1
        return ok;
    }

    bool reverse(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        size_t                              order_up    ,
        const CppAD::vector<double>&               taylor_x    ,
        const CppAD::vector<double>&               taylor_y    ,
        CppAD::vector<double>&                     partial_x   ,
        const CppAD::vector<double>&               partial_y   ) override
    {
        // TMB
		// ATOMIC_REVERSE  (X^-1*W[0])
		int n=sqrt((double)taylor_x.size());

        mat X(n, n);
        for(size_t i = 0; i < taylor_x.size(); i++) X.data()[i] = taylor_x[i];

        //mat X = vec2mat(taylor_x,n,n);
        mat invX = X.inverse();
        CppAD::vector<a_scalar> vecinvX = mat2vec(invX);

        for(size_t i=0; i<taylor_x.size(); i++) partial_x[i] = CppAD::Value(CppAD::Var2Par(vecinvX[i])) * partial_y[0];
        // ------
        return true;
    }

}; // End of atomic_logdet class


//a_scalar logdet(const a_matrix& arg)
//{
//  size_t nr = size_t( arg.rows() );
//  size_t ny = nr * nr;
//  // -------------------------------------------------------------------
//  // packed version of arg
//  CPPAD_TESTVECTOR(a_scalar) packed_arg(ny);
//    packed_arg[0] = a_scalar( nr );
//    for(size_t i = 0; i < ny; i++) packed_arg[i] = arg.data()[i];
//  // -------------------------------------------------------------------
//
//  static atomic_logdet afun_logdet("atomic_logdet");
//  CPPAD_TESTVECTOR(a_scalar) packed_result(1);
//  afun_logdet(packed_arg, packed_result);
//  // -------------------------------------------------------------------
//  return packed_result[0];
//}
//
//
//double logdet(const mat& X){
//    //TMB
//    mat LU=X.lu().matrixLU();    // Use Eigen LU decomposition
//    mat LUdiag = LU.diagonal();
//    return LUdiag.array().abs().log().sum();
//}

a_scalar logdet_internal(const CPPAD_TESTVECTOR(a_scalar)& packed_arg)
{
  // -------------------------------------------------------------------

  static atomic_logdet afun_logdet("atomic_logdet");
  CPPAD_TESTVECTOR(a_scalar) packed_result(1);
  afun_logdet(packed_arg, packed_result);
  // -------------------------------------------------------------------
  return packed_result[0];
}

template<auto ratc, auto catc, auto opts, auto mratc, auto mcatc>
a_scalar logdet(const Eigen::Matrix<CppAD::AD<double>, ratc, catc, opts, mratc, mcatc>& arg)
{
  size_t nr = size_t( arg.rows() );
  size_t ny = nr * nr;
  // -------------------------------------------------------------------
  // packed version of arg
  CPPAD_TESTVECTOR(a_scalar) packed_arg(ny);
    packed_arg[0] = a_scalar( nr );
    for(size_t i = 0; i < ny; i++) packed_arg[i] = arg.data()[i];
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  return logdet_internal(packed_arg);
}

template<auto ratc, auto catc, auto opts, auto mratc, auto mcatc>
double logdet(const Eigen::Matrix<double, ratc, catc, opts, mratc, mcatc>& X){
    //TMB
    mat LU=X.lu().matrixLU();    // Use Eigen LU decomposition
    mat LUdiag = LU.diagonal();
    return LUdiag.array().abs().log().sum();
}

class our_atomic_class_for_supplied_values : public CppAD::atomic_three<double> {
private:
    const Eigen::VectorXd function_value;
	const Eigen::MatrixXd gradient_value;

	typedef decltype(function_value.rows()) index_type;
	// This type is used instead of size_t as in other CppAD code, to enable correct signed/unsigned comparison
public:
    our_atomic_class_for_supplied_values(const std::string& name, const Eigen::VectorXd& func_val, const Eigen::MatrixXd& grad_val) :
    CppAD::atomic_three<double>(name), function_value(func_val), gradient_value(grad_val)
    { }
private:
// calculate type_y
    bool for_type(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        CppAD::vector<CppAD::ad_type_enum>&        type_y      ) override
    {
		for(index_type i = 0; i < function_value.rows(); i++)
			type_y[i] = CppAD::ad_type_enum::variable_enum ;
        return true;
    }
// forward mode routine called by CppAD
    bool forward(
        const CppAD::vector<double>&              parameter_x  ,
        const CppAD::vector<CppAD::ad_type_enum>& type_x       ,
        size_t                             need_y       ,
        size_t                             order_low    ,
        size_t                             order_up     ,
        const CppAD::vector<double>&              taylor_x     ,
        CppAD::vector<double>&                    taylor_y     ) override
    {
        index_type q1 = order_up + 1;
        //
        // we only implements up to the first order forward mode
        bool ok = order_up <=  1;
        if( ! ok )
            return ok;
        // ------------------------------------------------------------------
        // Zero forward mode.
        if( order_low <= 0 )
        {
			for(index_type i = 0; i < gradient_value.rows(); i++)
				taylor_y[i*q1+0] = function_value[i];
        }
        if( order_up <=  0 )
            return ok;
        // ------------------------------------------------------------------
        // First order forward mode.
        if( order_low <= 1 )
        {   for(index_type i = 0; i < gradient_value.rows(); i++)
                taylor_y[i*q1+1]  = 0.0;
        // ------------------------------------------------------------------
			for(index_type i = 0; i < gradient_value.rows(); i++){//
				for(index_type j = 0; j < gradient_value.cols(); j++){
					taylor_y[i*q1+1]  += gradient_value(i,j) * taylor_x[j*q1+1];
				}
			}
        }
        if( order_up <=  1 )
            return ok;
        // ------------------------------------------------------------------
        return ok;
    }
// reverse mode routine called by CppAD
    bool reverse(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        size_t                              order_up    ,
        const CppAD::vector<double>&               taylor_x    ,
        const CppAD::vector<double>&               taylor_y    ,
        CppAD::vector<double>&                     partial_x   ,
        const CppAD::vector<double>&               partial_y   ) override
    {
        size_t q1 = order_up + 1;
		// ------------------------------------------------------------------
        // we only implement the first order reverse mode
        bool ok = q1 <= 1;
        if( ! ok )
            return ok;
        for(index_type j = 0; j < gradient_value.cols(); j++)
                partial_x[j * q1 + 0] = 0.0;
        // ------------------------------------------------------------------
		for(index_type j = 0; j < gradient_value.cols(); j++){
			for(index_type i = 0; i < gradient_value.rows(); i++){
				partial_x[j*q1+0] += partial_y[i*q1+0] * gradient_value(i,j);
				//partial_x[j*q1+0] += 1.0;
			}
		}
        // ------------------------------------------------------------------
        return ok;
    }
};





template<class FunctionObj>
class atomic_class_for_supplied_functions : public CppAD::atomic_three<double>, protected FunctionObj{
// FunctionObj must provide methods
// VectorXd FunctionObj::value(VectorXd x) and
// MatrixXd FunctionObj::gradient(VectorXd x)
public:
	atomic_class_for_supplied_functions(const std::string& name, const FunctionObj& FO) :
		CppAD::atomic_three<double>(name), FunctionObj(FO) {}
private:
// calculate type_y
    bool for_type(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        CppAD::vector<CppAD::ad_type_enum>&        type_y      ) override
    {
		for(size_t i = 0; i < type_y.size(); i++)
			type_y[i] = CppAD::ad_type_enum::variable_enum ;
        return true;
    }
// forward mode routine called by CppAD
    bool forward(
        const CppAD::vector<double>&              parameter_x  ,
        const CppAD::vector<CppAD::ad_type_enum>& type_x       ,
        size_t                             need_y       ,
        size_t                             order_low    ,
        size_t                             order_up     ,
        const CppAD::vector<double>&              taylor_x     ,
        CppAD::vector<double>&                    taylor_y     ) override
    {
        size_t q1 = order_up + 1;
        //
        // we only implements up to the first order forward mode
        bool ok = order_up <=  1;
        if( ! ok )
            return ok;
		// ------------------------------------------------------------------
		Eigen::VectorXd theta(parameter_x.size());
		for(size_t i = 0; i < parameter_x.size(); i++)
			theta(i) = taylor_x[i*q1 + 0];
		Eigen::VectorXd function_value = FunctionObj::value(theta);
		Eigen::MatrixXd gradient_value = FunctionObj::gradient(theta);

        // ------------------------------------------------------------------
        // Zero forward mode.
        if( order_low <= 0 )
        {
			for(size_t i = 0; i < gradient_value.rows(); i++)
				taylor_y[i*q1+0] = function_value[i];
        }
        if( order_up <=  0 )
            return ok;
        // ------------------------------------------------------------------
                // First order forward mode.
        if( order_low <= 1 )
        {   for(size_t i = 0; i < gradient_value.rows(); i++)
                taylor_y[i*q1+1]  = 0.0;
        // ------------------------------------------------------------------
			for(size_t i = 0; i < gradient_value.rows(); i++){//
				for(size_t j = 0; j < gradient_value.cols(); j++){
					taylor_y[i*q1+1]  += gradient_value(i,j) * taylor_x[j*q1+1];
					//Rcpp::Rcout << "gradient_value " << taylor_y[i*q1+1] << std::endl;
				}
			}
        }
        if( order_up <=  1 )
            return ok;
        // ------------------------------------------------------------------
        return ok;
    }
// reverse mode routine called by CppAD
    bool reverse(
        const CppAD::vector<double>&               parameter_x ,
        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
        size_t                              order_up    ,
        const CppAD::vector<double>&               taylor_x    ,
        const CppAD::vector<double>&               taylor_y    ,
        CppAD::vector<double>&                     partial_x   ,
        const CppAD::vector<double>&               partial_y   ) override
    {
        size_t q1 = order_up + 1;
        Eigen::VectorXd theta(parameter_x.size());
		for(size_t i = 0; i < parameter_x.size(); i++)
			theta(i) = taylor_x[i*q1 + 0];
        Eigen::MatrixXd gradient_value = FunctionObj::gradient(theta);
		// ------------------------------------------------------------------
        // we only implement the first order reverse mode
        bool ok = q1 <= 1;
        if( ! ok )
            return ok;
        for(size_t j = 0; j < gradient_value.cols(); j++)
                partial_x[j * q1 + 0] = 0.0;
        // ------------------------------------------------------------------
		for(size_t j = 0; j < gradient_value.cols(); j++){
			for(size_t i = 0; i < gradient_value.rows(); i++){
				partial_x[j*q1+0] += partial_y[i*q1+0] * gradient_value(i,j);
				//Rcpp::Rcout << "gradient vec: " << partial_x[j*q1+0] << std::endl;
			}
		}
        // ------------------------------------------------------------------
        return ok;
    }
};

} // namespace atomic_funcs
} // namespace saddlepoint

#endif // ATOMIC_FUNCS_H_INCLUDED
