/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NumericPrecision_H_
#define model_NumericPrecision_H_

#include <cfloat>

namespace model
{
	
    template<typename T>
    struct NumericPrecision
    {
        static constexpr T epsilon=DBL_EPSILON;
    };

    template<>
    struct NumericPrecision<double>
    {
        static constexpr double epsilon=DBL_EPSILON;
    };

    
    template<>
    struct NumericPrecision<float>
    {
        static constexpr double epsilon=FLT_EPSILON;
    };

} // namespace model
#endif
