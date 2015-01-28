/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PERMUTATIONWITHOUTREPETITION_H_
#define model_PERMUTATIONWITHOUTREPETITION_H_

#include <model/Math/CompileTimeMath/Factorial.h>

namespace model
{

    template <int N, int k>
	struct PermutationWithoutRepetition
    {
//		enum{value=Factorial<N>::value / Factorial<(k>=0 && k<=N && N>=0)?(N-k):0>::value * (k>=0 && k<=N && N>=0)};
        static constexpr int value=Factorial<N>::value / Factorial<(k>=0 && k<=N && N>=0)?(N-k):0>::value * (k>=0 && k<=N && N>=0);

    };

}
#endif
