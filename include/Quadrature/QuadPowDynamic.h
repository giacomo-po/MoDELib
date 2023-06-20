/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_QUADPOWDYNAMIC_H_
#define model_QUADPOWDYNAMIC_H_

#include <assert.h>
#include <math.h>
#include <float.h>
#include <Eigen/Dense>
#include <QuadPow.h>
#include <QuadratureDynamic.h>



namespace model
{
    /**************************************************************************/
    template<
    short unsigned int pOrder,
    template <short unsigned int, size_t> class QuadratureRule,
    size_t...Orders>
    struct QuadPowDynamic;

    /**************************************************************************/
    template<short unsigned int pOrder,
    template <short unsigned int, size_t> class QuadratureRule,
    size_t qOrder,
    size_t... otherOrders>
    struct QuadPowDynamic<pOrder,QuadratureRule,qOrder,otherOrders...> : public QuadPowDynamic<pOrder,QuadratureRule,otherOrders...>
    {
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder+1>& uPow(const int& N)
        {
            return (N==qOrder)? QuadPow<pOrder,qOrder,QuadratureRule>::uPow : QuadPowDynamic<pOrder,QuadratureRule,otherOrders...>::uPow(N);
        }
        
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder>& duPow(const int& N)
        {
            return (N==qOrder)? QuadPow<pOrder,qOrder,QuadratureRule>::duPow : QuadPowDynamic<pOrder,QuadratureRule,otherOrders...>::duPow(N);
        }
        
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder-1>& dduPow(const int& N)
        {
            return (N==qOrder)? QuadPow<pOrder,qOrder,QuadratureRule>::dduPow : QuadPowDynamic<pOrder,QuadratureRule,otherOrders...>::dduPow(N);
        }
        
        static int lowerOrder(const int& k)
        {
            return QuadratureDynamic<1,QuadratureRule,qOrder,otherOrders...>::lowerOrder(k);
        }
        
        static int lowerOrder(const double& d)
        {
            return QuadratureDynamic<1,QuadratureRule,qOrder,otherOrders...>::lowerOrder(d);
        }
    };
    
    /**************************************************************************/
    template<short unsigned int pOrder,
    template <short unsigned int, size_t> class QuadratureRule,
    size_t qOrder>
    struct QuadPowDynamic<pOrder,QuadratureRule,qOrder>
    {
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder+1>& uPow(const int& N)
        {
            if(N!=qOrder)
            {
                throw std::runtime_error("quadrature order N not found in QuadPowDynamic");
            }
            return  QuadPow<pOrder,qOrder,QuadratureRule>::uPow;
        }
        
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder>& duPow(const int& N)
        {
            if(N!=qOrder)
            {
                throw std::runtime_error("quadrature order N not found in QuadPowDynamic");
            }
            return  QuadPow<pOrder,qOrder,QuadratureRule>::duPow;
        }
        
        static const Eigen::Matrix<double,Eigen::Dynamic,pOrder-1>& dduPow(const int& N)
        {
            if(N!=qOrder)
            {
                throw std::runtime_error("quadrature order N not found in QuadPowDynamic");
            }
            return  QuadPow<pOrder,qOrder,QuadratureRule>::dduPow;
        }
    };

    
} // namespace model
#endif
