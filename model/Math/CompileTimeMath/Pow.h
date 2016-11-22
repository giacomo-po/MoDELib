/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_POW_H_
#define model_POW_H_

namespace model
{
    
    /*!\brief A class template for compile-time calculation of powers.
     */
    template<int x, int n>
    struct Pow
    {
        static_assert (n>0, "n MUST BE > 0.");
        static constexpr int value = x * Pow<x,n-1>::value;
    };
    
    template<const int x>
    struct Pow<x,0>
    {
        static constexpr int value = 1;
    };
    
}
#endif
