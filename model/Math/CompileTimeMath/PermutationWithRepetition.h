/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
 
#ifndef model_PERMUTATIONWITHREPETITION_H_
#define model_PERMUTATIONWITHREPETITION_H_

#include <Eigen/Dense>
#include <CTM.h>

namespace model
{
    template <int _k>
    struct PermutationWithRepetition
    {
	/* The permutation P (n, k) =n Pk is the number of ways to form ordered sets of k objects taken from a pool of n objects if repetition is allowed. */
	
	/* This should be computed as N^K rows of K columns */
        
        static constexpr int k=_k;
//        static_assert(N>0,"N MUST BE >0.");
        static_assert(k>0,"k MUST BE >0.");
		
		
        template<typename T,int N>
		static Eigen::Matrix<T,k,CTM::pow(N,k)> permute(const Eigen::Matrix<T,1,N>& pool)
        {

            Eigen::Matrix<T,k-1,CTM::pow(N,k-1)> temp=PermutationWithRepetition<k-1>::permute(pool);
            Eigen::Matrix<T,k,CTM::pow(N,k)> temp1;
            for(size_t p=0;p<N;++p)
            {
                temp1.template block<  1,CTM::pow(N,k-1)>(0,p*CTM::pow(N,k-1))=Eigen::Matrix<T,1,CTM::pow(N,k-1)>::Constant(pool(p));
                temp1.template block<k-1,CTM::pow(N,k-1)>(1,p*CTM::pow(N,k-1))=temp;
            }
            
            return temp1;
		}
	};
    
    template <>
    struct PermutationWithRepetition<1>
    {
        /* The permutation P (n, k) =n Pk is the number of ways to form ordered sets of k objects taken from a pool of n objects if repetition is allowed. */
        
        /* This should be computed as N^K rows of K columns */
        
        static constexpr int k=1;
        
        
        template<typename T,int N>
        static const Eigen::Matrix<T,k,CTM::pow(N,k)>& permute(const Eigen::Matrix<T,1,N>& pool)
        {
            return pool;
        }
    };

}
#endif
