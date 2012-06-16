/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BINOMIAL_H_
#define model_BINOMIAL_H_

#include <model/Utilities/Factorial.h>

namespace model {
	
	//////////////////////////////////////////////////////////////////////////
	template< int N, int k>
	struct Binomial {
		enum { value = Factorial<N>::value / (Factorial<N-k>::value * Factorial<k>::value)};
		
//		template<typename T>
//		T operator+(const T & rhs){return value+rhs;}
//		
//		template<typename T>
//		T operator-(const T & rhs){return value-rhs;}
//		
//		template<typename T>
//		T operator*(const T & rhs){return value*rhs;}
//		
//		template<typename T>
//		T operator/(const T & rhs){return value/rhs;}
	};
	
//	template<>
//	struct Factorial<0> {
//		enum { value = 1 };
//	};
	
	//////////////////////////////////////////////////////////////////////////
} // end namespace model

#endif
