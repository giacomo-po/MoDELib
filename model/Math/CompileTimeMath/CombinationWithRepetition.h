/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CombinationWithRepetition_H_
#define model_CombinationWithRepetition_H_

#include <Eigen/Dense>
#include <CTM.h>

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
	template <int _N, int _k>
	struct CombinationWithRepetition
    {
        static constexpr int N=_N;
        static constexpr int k=_k;
        static_assert(N>0,"N MUST BE >0.");
        static_assert(k>0,"k MUST BE >0.");
        static constexpr int value=CTM::binomial(N+k-1,k);
      
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          *\returns a matrix having in each row the 
          */
            Eigen::Matrix<T,value,k> temp;
//            // Use Pascal rule to fill the combinations: bin(N,k) = bin(n-1,k-1) + bin(n-1,k)
            temp<<(Eigen::Matrix<T,k,CTM::binomial(N+k-2,k-1)>() <<	Eigen::Matrix<T,1,CTM::binomial(N+k-2,k-1)>().Constant(pool(0)),
                    CombinationWithRepetition<N,k-1>::template combine<T>(pool).transpose()).finished().transpose(),
                    CombinationWithRepetition<N-1,k>::template combine<T>(pool.template segment<N-1>(1));
            
            return temp;
            
        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // End recursion in n at n=1
    template <int k>
	struct CombinationWithRepetition<1,k>
    {        
        
        static_assert(k>0,"k MUST BE >0.");
        static constexpr int N=1;
		static constexpr int value=CTM::binomial(N+k-1,k);
        
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          */
            return Eigen::Matrix<T,value,k>::Constant(pool(N-1));
        }
        
	};

    /**************************************************************************/
    /**************************************************************************/
    // End recursion in k at k=1
    template <int _N>
	struct CombinationWithRepetition<_N,1>
    {
                
        static constexpr int N=_N;
        static_assert(N>0,"N MUST BE >0.");
        static constexpr int k=1;
		static constexpr int value=CTM::binomial(N+k-1,k);
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          */
            return pool;            
        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // Case N=1, k=1
    template <>
	struct CombinationWithRepetition<1,1>
    {
        static constexpr int N=1;
        static constexpr int k=1;
		static constexpr int value=CTM::binomial(N+k-1,k);
        
        template <typename T>
        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        {/*!@param[in] pool a row vector of N values
          */
            return pool;
        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // Case k=0
    template <int _N>
	struct CombinationWithRepetition<_N,0>
    {
        static constexpr int N=_N;
        static_assert(N>0,"N MUST BE >0.");
        static constexpr int k=0;
		static constexpr int value=CTM::binomial(N+k-1,k);
        
//        template <typename T>
//        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
//        {/*!@param[in] pool a row vector of N values
//          */
//            return pool;
//        }
        
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // Case N=1, k=0
    template <>
	struct CombinationWithRepetition<1,0>
    {
        
        static constexpr int N=1;
//        static_assert(N>0,"N MUST BE >0.");
        static constexpr int k=0;
		static constexpr int value=CTM::binomial(N+k-1,k);
        
        //        template <typename T>
        //        static Eigen::Matrix<T,value,k> combine(const Eigen::Matrix<T,1,N>& pool)
        //        {/*!@param[in] pool a row vector of N values
        //          */
        //            return pool;
        //        }
        
	};

}
#endif
