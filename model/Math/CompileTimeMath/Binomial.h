/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BINOMIAL_H_
#define model_BINOMIAL_H_

#include <model/Math/CompileTimeMath/Factorial.h>

namespace model
{
	/*\brief A compile-time class template for computation of binomial coefficients
     */
	template< int N, int k>
	struct Binomial
    {
//		enum { value = Factorial<N>::value / (Factorial<N-k>::value * Factorial<k>::value)};
        static constexpr int value = factorial(N) / (factorial(N-k) * factorial(k));
		
	};

} // end namespace model

#endif
