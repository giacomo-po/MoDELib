/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GAUSSLEGENDRE_H_
#define model_GAUSSLEGENDRE_H_


#include <assert.h>





namespace model {
	
	/**************************************************/
	/* GaussLegendre: general case                    */
	/**************************************************/
	/*! \brief Class template defining the GaussLegendre rules for determination 
	 *	of quadrature abscissas and weights.
	 */
	template<short unsigned int dim, short unsigned int qOrder>
	struct GaussLegendre{
		
//        static_assert(false,"GaussLegendre: CASE NOT IMPLEMENTED.");
		
		
		GaussLegendre(){
			assert(0 && "GaussLegendre: dimensionality not implemented.");
		}
		
	};
    
    
	
	
	//////////////////////////////////////////////////////////////
} // namespace model


#include <model/Quadrature/GaussLegendre/includeGaussLegendre1D.h>
#include <model/Quadrature/GaussLegendre/includeGaussLegendre2D.h>
#include <model/Quadrature/GaussLegendre/includeGaussLegendre3D.h>

#endif

