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
	/*\brief A compile-time implementation of binomial coefficients
     */
    
    constexpr int binomial(int N, int k)
    {
        return factorial(N) / (factorial(N-k) * factorial(k));
    }
    
}

#endif
