/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NONUNIFORMOPEN_H_
#define model_NONUNIFORMOPEN_H_

#include <cmath>
#include <assert.h>
#include <Eigen/Dense>

namespace model
{
	
	/**************************************************/
	/* NonUniformOpen: general case                    */
	/**************************************************/
	/*! \brief Class template defining the quadrature rules for determination 
	 *	of quadrature abscissas and weights implemented specifically for discretization of segments with QPs.
	 */
	template<short unsigned int dim, size_t qOrder>
	struct NonUniformOpen
    {
		NonUniformOpen()
        {
			assert(0 && "NonUniformOpen: dimensionality not implemented.");
		}
		
	};
	
	
	/**************************************************/
	/* UniformOpen: template specialization dim=1   */
	/**************************************************/
	/*!\brief NonUniformOpen<1,N> places quadrature points
     * at the segments such that their density is high near the end points 
     * low away from the end points
	 */
	template<size_t qOrder>
	struct NonUniformOpen<1,qOrder>
    {
		static Eigen::Matrix<double,2,qOrder> abcsissasAndWeights()
        {
			Eigen::Matrix<double,2,qOrder> U;
            for(size_t k=0;k<qOrder;++k)
            {
                // U(0,k)= 0.5/qOrder+k*1.0/qOrder;	//1.0 forces double conversion
                U(0,k)= 0.5*(cos((2.0*double (k)+1.0)*M_PI/2/qOrder)+1.0);
                // std::cout<<"Absicca"<<U(0,k)<<std::endl;
                U(1,k)= 1.0/qOrder;
            }
			return U;
		}
		
	};
} // namespace model
#endif

