/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINENODE_H_
#define model_SPLINENODE_H_

//#include <vector>
//#include <assert.h>
//#include <iterator>
//#include <Eigen/Dense>
//#include <model/Network/NetworkNode.h>
//#include <model/Network/Operations/EdgeExpansion.h>
//#include <model/Math/CompileTimeMath/Pow.h>

namespace model
{
	
	/**************************************************************************/
	/* class template SplineNodeBase: general case                            */
	/**************************************************************************/
	template <typename Derived,	short unsigned int dim, short unsigned int corder,  typename InterpolationType>
	class SplineNode
    {
		
	public:
//		SplineNode()
//        {
//			std::cout<<"This Node type has not been implemented"<<std::endl;
//			assert(0);
//		};
	};
		
}

//#include "model/Geometry/Splines/SplineNode_Hermite.h"
#include "model/Geometry/Splines/SplineNode_CatmullRom.h"

#endif
