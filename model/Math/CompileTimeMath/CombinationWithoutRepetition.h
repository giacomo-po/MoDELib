/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_COMBINATIONWITHOUTREPETITION_H_
#define model_COMBINATIONWITHOUTREPETITION_H_


#include <model/CompileTimeMath/Factorial.h>
#include <model/CompileTimeMath/PermutationWithoutRepetition.h>

namespace model {
	
	
	template <int N, int k>
	struct CwoR{
		enum{value=PwoR<N,k>::value / Factorial<k>::value};
	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace ctmath

#endif
