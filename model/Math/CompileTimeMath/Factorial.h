/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_FACTORIAL_H_
#define model_FACTORIAL_H_

namespace model {
	
	//////////////////////////////////////////////////////////////////////////
	template< int N >
	struct Factorial {
		enum { value = Factorial<(N>0)?(N-1):0>::value * ((N>0)? N:1)};
	};
	
	template<>
	struct Factorial<0> {
		enum { value = 1 };
	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace model

#endif
