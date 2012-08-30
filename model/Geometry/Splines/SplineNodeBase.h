/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINENODEBASE_H_
#define model_SPLINENODEBASE_H_

#include <vector>
#include <assert.h>
#include <iterator>
#include <boost/tuple/tuple.hpp>
#include <Eigen/Dense>
#include <model/Network/NetworkNode.h>
//#include <model/Utilities/InitializationData.h>
//#include <model/Network/EdgeExpansion.h>
#include <model/Network/Operations/EdgeExpansion.h>
#include <model/Math/CompileTimeMath/Pow.h>


namespace model {
	
	/*******************************************************************************/
	/* class template SplineNodeBase: general case                                 */
	/*******************************************************************************/
	template <typename Derived,	short unsigned int dim, short unsigned int corder,  typename InterpolationType>
	class SplineNodeBase {
		
	public:
		SplineNodeBase(){
			std::cout<<"Sorry, this Node type has not been implemented"<<std::endl;
			assert(0);
		};
	}; 
	
	#include "model/Geometry/Splines/SplineNodeBase_Hermite.h"
	#include "model/Geometry/Splines/SplineNodeBase_CatmullRom.h"
	
}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#endif