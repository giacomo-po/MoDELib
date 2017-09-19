/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_FACTORIAL_H_
#define model_FACTORIAL_H_

namespace model
{
	
    /*!\brief A compile-time implementation of the factorial function         */
    constexpr int factorial(int n)
    {
        return n <= 1? 1 : (n * factorial(n - 1));
    }
    
}
#endif
