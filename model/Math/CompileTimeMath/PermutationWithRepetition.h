/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_PERMUTATIONWITHREPETITION_H_
#define model_PERMUTATIONWITHREPETITION_H_

#include <model/CompileTimeMath/Pow.h>

namespace model
{
	
	template<int N, int K>
	struct PermutationWithRepetition
    {
	/* The permutation P (n, k) =n Pk is the number of ways to form ordered sets of k objects taken from a pool of n objects if repetition is allowed. */
	
	/* This should be computed as N^K rows of K columns */
		
		Eigen::Matrix<int,Pow<N,K>::value,K> permutations;
		
		PermutationWithRepetition(const Eigen::Matrix<int,N,1>& pool)
        {
		
		}
	};

}
#endif
