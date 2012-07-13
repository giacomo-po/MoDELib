/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_POW_H_
#define model_POW_H_

//#include <assert.h>
#include <boost/static_assert.hpp> // to be repaced with static_assert in C++0x

namespace model {
	
	//////////////////////////////////////////////////////////////////////////
	template<const int x, const int n>
	struct Pow {
		BOOST_STATIC_ASSERT(n>0);
		enum { value = x * Pow<x,n-1>::value};
		
		
//		/////////////////////////////////////////////////////////
//		// cast operator: allows a Pow<n> instance to be used as a regular int
//		operator const int & () const {
//			return static_cast<const int>(value);
//		}
		
	};
	
	template<const int x>
	struct Pow<x,0> {
		enum { value = 1};
	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace model

#endif
