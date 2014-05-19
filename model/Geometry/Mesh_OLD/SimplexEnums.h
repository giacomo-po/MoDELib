/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_SIMPLEXENUMS_H_
#define  model_SIMPLEXENUMS_H_

namespace model {
	
	template<short unsigned int order>
	struct SimplexEnums{
		enum {nVertices=SimplexEnums<order-1>::nVertices+1};
		enum {nEdges=SimplexEnums<order-1>::nEdges+SimplexEnums<order-1>::nVertices};
	};
	
	template<>
	struct SimplexEnums<0>{
		enum {nVertices=1};
		enum {nEdges=0};
	};	
	
}
#endif
