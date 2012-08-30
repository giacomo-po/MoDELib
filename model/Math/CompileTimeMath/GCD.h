/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GCD_H_
#define model_GCD_H_

namespace model {
	
	//////////////////////////////////////////////////////////////////////////
	template< size_t A, size_t B >
	struct GCD {
		enum { value = GCD< B, A % B >::value };
	};
	
	template< size_t A >
	struct GCD<A,0> {
		enum { value = A };
	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace model

#endif