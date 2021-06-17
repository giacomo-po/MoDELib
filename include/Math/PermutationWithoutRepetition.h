/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PERMUTATIONWITHOUTREPETITION_H_
#define model_PERMUTATIONWITHOUTREPETITION_H_

//#include <Factorial.h>
//#include <Pow.h>
#include <CTM.h>

namespace model
{
    
    
    
    /**************************************************************************/
    /**************************************************************************/
    template <int k>
    struct PermutationWithoutRepetition
    {
        //		enum{value=Factorial<N>::value / Factorial<(k>=0 && k<=N && N>=0)?(N-k):0>::value * (k>=0 && k<=N && N>=0)};
        static constexpr int value(int N)
        {
            return CTM::factorial(N) / CTM::factorial((k>=0 && k<=N && N>=0)?(N-k):0) * (k>=0 && k<=N && N>=0);
        }
        
        template <typename T,int M>
        static Eigen::Matrix<T,1,M-1> reducePool(const Eigen::Matrix<T,1,M>& pool,const int& p)
        {
            Eigen::Matrix<T,1,M-1> temp;
            int n=0;
            for (int i=0;i<M;++i)
            {
                if(i!=p)
                {
                    temp(n)=pool(i);
                    n++;
                }
            }
            return temp;
        }
        
        template<typename T,int N>
        static Eigen::Matrix<T,k,value(N)> permute(const Eigen::Matrix<T,N,1>& pool)
        {
            return permute(pool.transpose().eval());
        }
        
        
        template<typename T,int N>
        static Eigen::Matrix<T,k,value(N)> permute(const Eigen::Matrix<T,1,N>& pool)
        {
            assert(N>=k);
            
            
            Eigen::Matrix<T,k,value(N)> temp1;
            const int c(PermutationWithoutRepetition<k-1>::value(N-1));
            for(size_t p=0;p<N;++p)
            {
                temp1.template block<  1,c>(0,p*c)=Eigen::Matrix<T,1,c>::Constant(pool(p));
                temp1.template block<k-1,c>(1,p*c)=PermutationWithoutRepetition<k-1>::permute(reducePool(pool,p));
            }
            
            return temp1;
        }
        
        
        template<typename T,int N>
        static Eigen::Matrix<T,k,value(N)*CTM::pow(2,k)> permuteWithPlusMinusSign(const Eigen::Matrix<T,N,1>& pool)
        {
            return permuteWithPlusMinusSign(pool.transpose().eval());
        }
        
        template<typename T,int N>
        static Eigen::Matrix<T,k,value(N)*CTM::pow(2,k)> permuteWithPlusMinusSign(const Eigen::Matrix<T,1,N>& pool)
        {
            assert(N>=k);
            
            
            Eigen::Matrix<T,k,value(N)*CTM::pow(2,k)> temp1;
            const int c(PermutationWithoutRepetition<k-1>::value(N-1)*CTM::pow(2,k-1));
            for(size_t p=0;p<N;++p)
            {
                temp1.template block<  1,c>(0,(2*p+0)*c)=Eigen::Matrix<T,1,c>::Constant( pool(p));
                temp1.template block<  1,c>(0,(2*p+1)*c)=Eigen::Matrix<T,1,c>::Constant(-pool(p));
                
                Eigen::Matrix<T,k-1,c> temp=PermutationWithoutRepetition<k-1>::permuteWithPlusMinusSign(reducePool(pool,p));
                temp1.template block<k-1,c>(1,(2*p+0)*c)=temp;
                temp1.template block<k-1,c>(1,(2*p+1)*c)=temp;
                //                std::cout<<<<std::endl;
            }
            
            return temp1;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PermutationWithoutRepetition<1>
    {
        static constexpr int k=1;
        //		enum{value=Factorial<N>::value / Factorial<(k>=0 && k<=N && N>=0)?(N-k):0>::value * (k>=0 && k<=N && N>=0)};
        static constexpr int value(int N)
        {
            return CTM::factorial(N) / CTM::factorial((k>=0 && k<=N && N>=0)?(N-k):0) * (k>=0 && k<=N && N>=0);
        }
        
        
        
        template<typename T,int N>
        static const Eigen::Matrix<T,k,value(N)>& permute(const Eigen::Matrix<T,1,N>& pool)
        {
            assert(N>=k);
            return pool;
        }
        
        template<typename T,int N>
        static const Eigen::Matrix<T,k,2*N> permuteWithPlusMinusSign(const Eigen::Matrix<T,1,N>& pool)
        {
            return (Eigen::Matrix<T,k,2*N>()<<pool,-pool).finished();
        }
        
    };
    
}
#endif
