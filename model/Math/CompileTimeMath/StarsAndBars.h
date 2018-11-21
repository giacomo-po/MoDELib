/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StarsAndBars_H_
#define model_StarsAndBars_H_

#include <Eigen/Dense>
#include <CombinationWithRepetition.h>

namespace model
{
    /**************************************************************************/
    /**************************************************************************/
	/*!Class template that computes the possible sets of N integers whose sum
     * is a given integer k. This is the solution of the Diophantine equation
     *  \f$x_1+\ldots + x_N=k\f$.
     */
	template <int N, int k,int j=0>
	struct StarsAndBars
    {
        static_assert(N>0,"N MUST BE >0.");
        static_assert(k>0,"k MUST BE >0.");
        static_assert(j>=0,"j MUST BE >0.");
        static_assert(j<k,"j MUST BE <k.");
        
//        static constexpr int N=_N; // number of positive integers
//        static constexpr int k=_k; // value of the sum of the N integers
//        static constexpr int j=_j; // value of the sum of the N integers
        static constexpr int r1=CombinationWithRepetition<N,k-j-1>::value;
        static constexpr int r2=CombinationWithRepetition<N-1,k-j>::value;
        
        static Eigen::Matrix<int,r1+r2,N> sAb()
        {/*!@param[in] pool a row vector of N values
          *\returns a matrix having in each row one of the combinatins of 
          * integers whose sum is k
          */
            return (Eigen::Matrix<int,r1+r2,N>() << (StarsAndBars<N,k,j+1>::sAb()),
            (Eigen::Matrix<int,N,r2>() <<	Eigen::Matrix<int,1,r2>().Constant(j),
             StarsAndBars<N-1,k-j,0>::sAb().transpose()).finished().transpose()).finished();
        }
        
	};
    
    
    /**************************************************************************/
    /**************************************************************************/
    // End recursion in j at j=k
    template <int N, int k>
	struct StarsAndBars<N,k,k>
    {
        static_assert(N>0,"N MUST BE >0.");
//        static constexpr int N=_N;
//        static constexpr int k=_k;
//        static constexpr int j=_k;
        static constexpr int r=CombinationWithRepetition<N-1,k-k>::value; // this is = 1, since C(n,0)=1
        
        static Eigen::Matrix<int,r,N> sAb()
        {
            Eigen::Matrix<int,r,N> temp;
            temp << k, StarsAndBars<N-1,0,0>::sAb();
            return temp;
        }
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // End recursion in k at k=0, (which implies j=0)
    template <int N>
	struct StarsAndBars<N,0,0>
    {
        static_assert(N>0,"N MUST BE >0.");
//        static constexpr int N=_N;
//        static constexpr int k=0;
//        static constexpr int j=k;
		static constexpr int r=CombinationWithRepetition<N,0>::value; // this is = 1, since C(n,0)=1
        
        static Eigen::Matrix<int,r,N> sAb() // this is an 1xN matrix
        {
            return Eigen::Matrix<int,r,N>::Zero();
        }
	};
    
    /**************************************************************************/
    /**************************************************************************/
    // End recursion in n at n=1
    template <int k,int j>
	struct StarsAndBars<1,k,j>
    {
        static_assert(k>0,"k MUST BE >0.");
//        static constexpr int N=1;
		static constexpr int r=CombinationWithRepetition<1,k-j>::value; // this is = 1, since C(1,k-j)=1
        
        static Eigen::Matrix<int,r,1> sAb() // this is a 1x1 matrix
        {
            return Eigen::Matrix<int,r,1>::Constant(k-j);
        }
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Template Specialization Case <1,k,k>
    template <int k>
    struct StarsAndBars<1,k,k>
    {
        static_assert(k>0,"k MUST BE >0.");
//        static constexpr int N=1;
//        static constexpr int j=k;
        static constexpr int r=CombinationWithRepetition<1,0>::value; // this is = 1, since C(1,k-j)=1
        
        static Eigen::Matrix<int,r,1> sAb() // this is a 1x1 matrix
        {
            return Eigen::Matrix<int,r,1>::Constant(0);
        }
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Template Specialization Case <1,0,0>
    template <>
    struct StarsAndBars<1,0,0>
    {        
//        static constexpr int N=1;
//        static constexpr int k=0;
//        static constexpr int j=k;
        static constexpr int r=CombinationWithRepetition<1,0>::value; // this is = 1, since C(1,0)=1
        
        static Eigen::Matrix<int,r,1> sAb() // this is a 1x1 matrix
        {
            return Eigen::Matrix<int,r,1>::Constant(0);
        }
        
	};
}
#endif
