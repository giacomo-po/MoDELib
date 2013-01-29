/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QUADRATURE_H_
#define model_QUADRATURE_H_

#include <assert.h>
#include <Eigen/Dense>
#include <model/Quadrature/GaussLegendre/GaussLegendre.h>
#include <model/Quadrature/UniformOpen.h>



namespace model {	
	
	template<short unsigned int dim, short unsigned int qOrder>
	struct VectorDimTypeSelector{
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		
		static const VectorDim abscissa(const short unsigned int & k, const Eigen::Matrix<double,dim,qOrder>& abscissas){
			return abscissas.col(k);
		}
	};
	
	template<short unsigned int qOrder>
	struct VectorDimTypeSelector<1,qOrder>{
		typedef double VectorDim;
		
		static const VectorDim& abscissa(const short unsigned int & k, const Eigen::Matrix<double,1,qOrder>& abscissas){
			return abscissas(k);
		}
	};
	
	
	/*******************************************************************************************************/
	/*******************************************************************************************************/
	/*! \brief A static class template for generating quadrature abcsissas and weights for arbitrary-order
	 *  and QuadratureRule and dimension.
	 */
	template<short unsigned int dim, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule = GaussLegendre>
	struct Quadrature {
        
		
		typedef typename VectorDimTypeSelector<dim,qOrder>::VectorDim VectorDim;
		
//		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
				
		// a static vector storing abcsissa values
		static const Eigen::Matrix<double,dim,qOrder>	abscissas;
		
		// a static vector storing weight values
		static const Eigen::Matrix<double,1,  qOrder>	weights;		
		
		/* weight ********************************************/
		static const double& weight(const short unsigned int & k) {
			return weights(k);
		}
		
		/* abscissa ******************************************/
		static const VectorDim abscissa(const short unsigned int & k) {
			return VectorDimTypeSelector<dim,qOrder>::abscissa(k,abscissas);
		}
		

		/* integrate_function *******************************/
		template<typename IntegrandType>
		static void integrate(IntegrandType (*fp)(const double&), IntegrandType &intgrl) {
			/*! Integrates function in range [0,1]
			 * @param[in]  fp     the function pointer.
			 * @param[out] intgrl the integral to be returned.
			 *
			 * Sample code: \code
			 * #include <iostream>
			 * #include <Eigen/Dense>
			 * #include <model/Quadrature/Quadrature.h>
			 * 
			 * Eigen::Matrix<double,3,3> integrand(const double& u) {
			 *	return Eigen::Matrix<double,3,3>::Ones()*u*u;
			 * }
			 *	
			 * int main(){
			 *	Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
			 *	Q.integrate(this,&integrand,I);
			 * 	std::cout<< I <<std::endl;	
			 * 	return 0;
			 * }
			 * \endcode
			 */
			for(int k=0;k<qOrder;++k){
				intgrl+=(*fp)(k) * weights(k);
			}			
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const, const Args&...args) {
			/*! Integrates non-static member function of arbitrary class in range [0,1]
			 * @param[out] intgrl The memory area to copy to.
			 * @param[in]  mfp    The member function pointer.
 			 * @param[in]  C    The pointer to ...
			 *
			 * Sample code: \code
			 * #include <iostream>
			 * #include <Eigen/Dense>
			 * #include <model/Quadrature/Quadrature.h>
			 * 
	 		 * class SomeClass{
			 *	
			 *	Eigen::Matrix<double,3,3> integrand(const double& u) const {
			 *		return Eigen::Matrix<double,3,3>::Ones()*u*u;
			 *	}
			 *	
			 *	public:
			 *		
			 *	template <short unsigned int qOrder>
			 *	void compute_integral(){
			 *		Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
			 *		model::Quadrature<1,qOrder,model::GaussLegendre>::integrate(this,&SomeClass::integrand,I);
			 *		std::cout<< I <<std::endl;
			 *	}
			 *		
			 * };
			 *
			 * int main(){
			 * 	SomeClass sc;
			 *  sc.compute_integral<16>();
			 * 	return 0;
			 * }
			 * \endcode
			 */
			for(int k=0;k<qOrder;++k){
				intgrl+=(C->*mfp)(abscissa(k),args...) * weights(k);
			}			
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const, const Args&...args) {
			/*! Integrates non-static member function of arbitrary class in range [0,1]
			 * @param[out] intgrl The memory area to copy to.
			 * @param[in]  mfp    The member function pointer.
 			 * @param[in]  C    The pointer to ...
			 *
			 * Sample code: \code
			 * #include <iostream>
			 * #include <Eigen/Dense>
			 * #include <model/Quadrature/Quadrature.h>
			 * 
	 		 * class SomeClass{
			 *	
			 *	Eigen::Matrix<double,3,3> integrand(const double& u) const {
			 *		return Eigen::Matrix<double,3,3>::Ones()*u*u;
			 *	}
			 *	
			 *	public:
			 *		
			 *	template <short unsigned int qOrder>
			 *	void compute_integral(){
			 *		Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
			 *		model::Quadrature<1,qOrder,model::GaussLegendre>::integrate(this,&SomeClass::integrand,I);
			 *		std::cout<< I <<std::endl;
			 *	}
			 *		
			 * };
			 *
			 * int main(){
			 * 	SomeClass sc;
			 *  sc.compute_integral<16>();
			 * 	return 0;
			 * }
			 * \endcode
			 */
					
			for(int k=0;k<qOrder;++k){
				intgrl+=(C->*mfp)(k,args...) * weights(k);
			}			
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(AnyClass* const C, void (AnyClass::*mfp)(const VectorDim&, const Args&...), const Args&...args) {
			/*! Execute non-static member function of arbitrary class, at quadrature points in range [0,1]
			 
			 */
			for(int k=0;k<qOrder;++k){
				(C->*mfp)(abscissa(k),args...);
			}			
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(AnyClass* const C, void (AnyClass::*mfp)(const int&, const Args&...), const Args&...args) {
			/*! Execute non-static member function of arbitrary class, at quadrature points in range [0,1]
			 
			 */
			for(int k=0;k<qOrder;++k){
				(C->*mfp)(k,args...);
			}			
		}
		
		
	};
	
	
	//! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
	template<short unsigned int dim, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,dim,qOrder> Quadrature<dim,qOrder,QuadratureRule>::abscissas=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<dim,qOrder>(0,0);
		
	//! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
	template<short unsigned int dim, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,1,qOrder> Quadrature<dim,qOrder,QuadratureRule>::weights=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<1,qOrder>(dim,0);

	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

// /solenv/unxmacxp/bin/init-static-template-data
