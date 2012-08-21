/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UNIFORMOPEN_H_
#define model_UNIFORMOPEN_H_


#include <assert.h>
//#include <model/Quadrature/QuadratureRecursiveRelation.h>

namespace model {
	
	/**************************************************/
	/* UniformOpen: general case                    */
	/**************************************************/
	/*! \brief Class template defining the GaussLegendre rules for determination 
	 *	of quadrature abscissas and weights.
	 */
	template<short unsigned int dim, short unsigned int qOrder>
	struct UniformOpen{
		
		
		
		UniformOpen(){
			assert(0 && "GaussLegendre: dimensionality not implemented.");
		}
		
	};
	
	
	/**************************************************/
	/* UniformOpen: template specialization dim=1   */
	/**************************************************/
	/*! \brief UniformOpen
	 *	
	 */
	template<short unsigned int qOrder>
	struct UniformOpen<1,qOrder> {
		static Eigen::Matrix<double,2,qOrder> abcsissasAndWeights(){
			Eigen::Matrix<double,2,qOrder> U;
            for(int k=0;k<qOrder;++k){
                U(0,k)= 0.5/qOrder+k*1.0/qOrder;		
                U(1,k)= 1.0/qOrder;
            }
            
           // std::cout<< "Quadrature rule UniformOpen<1,"<<qOrder<<"> (weights in last row):"<<std::endl;
			//std::cout<<std::setprecision(15)<<std::scientific<<U<<std::endl;
            
			return U;
		}
		
	};
		
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

