/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_QUADPOW_H_
#define model_QUADPOW_H_

#include <assert.h>
#include <math.h>
#include <Eigen/Dense>

#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/GaussLegendre.h>


namespace model {
	

	
	template<short unsigned int pOrder, short unsigned int qOrder, template<short unsigned int,short unsigned int> class QuadratureRule=GaussLegendre>
	struct QuadPow {
	
	//	static const Eigen::Matrix<double,1,qOrder> A;
		static const Eigen::Matrix<double,qOrder,pOrder+1>   uPow;
		static const Eigen::Matrix<double,qOrder,pOrder  >  duPow;
		static const Eigen::Matrix<double,qOrder,pOrder-1> dduPow;
		
		
		static Eigen::Matrix<double,qOrder,pOrder+1> uPowFill(){

//			Quadrature<1,qOrder,QuadratureRule> Q;
//			Eigen::Matrix<double,1,qOrder> B = Q.abscissas;
//			std::cout<<Quadrature<1,qOrder,QuadratureRule>::abscissas<<std::endl;

			//std::cout<<Quadrature<1,qOrder,QuadratureRule>::abscissas_(0)<<std::endl;
			if ( Quadrature<1,qOrder>::abscissas.squaredNorm() ){}
			else{
				std::cout<<"THIS IS A BUG WITH GCC IN MAC osX."<<std::endl;
				std::cout<<"To fix it create make a call to the coresponding EXPLICIT Quadrature<1,qOrder>::abscissas before using QuadPow<pOrder,qOrder>. Eg, call Quadrature<1,8>::abscissa. "<<std::endl;
				assert(0);
			}
			

	
			Eigen::Matrix<double,qOrder,pOrder+1> temp;
			
			for (int i=0; i<qOrder; ++i){
				for (int j=0; j<pOrder+1; ++j){
					temp(i,j)=std::pow(Quadrature<1,qOrder,QuadratureRule>::abscissas(i),j);
				}
			}
			
			std::cout<< "Computing powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			
			
			return temp;
		}
		
		static Eigen::Matrix<double,qOrder,pOrder> duPowFill(){
			
			Eigen::Matrix<double,qOrder,pOrder> temp;
			
			for (int i=0; i<qOrder; ++i){
				for (int j=1; j<pOrder+1; ++j){
					temp(i,j-1)=j*std::pow(Quadrature<1,qOrder,QuadratureRule>::abscissas(i),j-1);
				}
			}
			std::cout<< "Computing 1-st derivative of powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			return temp;
		}
		
		static Eigen::Matrix<double,qOrder,pOrder-1> dduPowFill(){
			
			Eigen::Matrix<double,qOrder,pOrder-1> temp;
			
			for (int i=0; i<qOrder; ++i){
				for (int j=2; j<pOrder+1; ++j){
					temp(i,j-2)=j*(j-1)*std::pow(Quadrature<1,qOrder,QuadratureRule>::abscissas(i),j-2);
				}
			}
			std::cout<< "Computing 2-nd derivative of powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			return temp;
		}
	
	};
	
//	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
//	const Eigen::Matrix<double,1,qOrder> QuadPow<pOrder,qOrder,QuadratureRule>::A=Quadrature<1,qOrder,QuadratureRule>::abscissas;
	
	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder+1> QuadPow<pOrder,qOrder,QuadratureRule>::uPow=QuadPow<pOrder,qOrder,QuadratureRule>::uPowFill();
	
	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder  > QuadPow<pOrder,qOrder,QuadratureRule>::duPow=QuadPow<pOrder,qOrder,QuadratureRule>::duPowFill();
	
	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder-1> QuadPow<pOrder,qOrder,QuadratureRule>::dduPow=QuadPow<pOrder,qOrder,QuadratureRule>::dduPowFill();

	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif
