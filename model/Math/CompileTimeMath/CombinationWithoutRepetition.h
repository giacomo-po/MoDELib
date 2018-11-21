/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CombinationWithoutRepetition_H_
#define model_CombinationWithoutRepetition_H_

#include <Eigen/Dense>
#include <Binomial.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
	/*!Class template that computes C(N,k), the combination without repetition 
     * of k values taken from a pool of N values. 
     *
     * C(N,k) represents the number of ways to form un-ordered sets of k objects 
     * taken from a pool of N objects if repetition is not allowed.
     */
	template <int N, int k>
	struct CombinationWithoutRepetition
    {
        static_assert(N>0,"N MUST BE >0.");
        static_assert(k>0,"k MUST BE >0.");
        static_assert(k<=N,"k MUST BE <=N.");
        static constexpr int value=Binomial<N,k>::value;
//		enum{value=Binomial<N,k>::value};
        
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          *\returns a matrix having in each row the 
          */
            Eigen::Matrix<T,value,k> temp;
            // Use Pascal rule to fill the combinations: bin(N,k) = bin(n-1,k-1) + bin(n-1,k)
            temp<<(Eigen::Matrix<T,k,Binomial<N-1,k-1>::value>() <<	Eigen::Matrix<T,1,Binomial<N-1,k-1>::value>().Constant(pool(0)),
                   CombinationWithoutRepetition<N-1,k-1>::template combine<T>(pool.template segment<N-1>(1)).transpose()).finished().transpose(),
            CombinationWithoutRepetition<N-1,k>::template combine<T>(pool.template segment<N-1>(1));

            
            return temp;
            
        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // SPECIALIZATION FOR THE TERM BIN(N-1,k). The stop condition is reached at BIN(k,k)
    template <int k>
	struct CombinationWithoutRepetition<k,k>
    {
        static_assert(k>0,"k MUST BE >0.");
        enum{N=k};
		enum{value=Binomial<N,k>::value};
        
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          */
            return pool;
        }
        
	};
	
    /**************************************************************************/
    /**************************************************************************/
    // SPECIALIZATION FOR THE TERM BIN(N-1,k-1). The stop condition is reached at BIN(N,1)
    template <int N>
	struct CombinationWithoutRepetition<N,1>
    {
        static_assert(N>0,"N MUST BE >0.");
        enum{k=1};
		enum{value=Binomial<N,k>::value};
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          */
            return pool;            
        }
        
	};

    /**************************************************************************/
} // end namespace ctmath

#endif
