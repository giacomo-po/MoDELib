/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakDiff_H_
#define model_LinearWeakDiff_H_

//#include <AreSameType.h>
#include <type_traits> // std::is_same
#include <LinearWeakExpression.h>
#include <ExpressionRef.h>

namespace model
{
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief A class template that performs the sum of two TrialExpressionBase.
     */
    template <typename T1,typename T2>
    struct  LinearWeakDiff : public LinearWeakExpression<LinearWeakDiff<T1,T2> >
    {
        
        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING TrialExpressionBase OF DIFFERENT TRIALFUNCTIONS.");
        typedef typename T1::TrialFunctionType TrialFunctionType;
        
        //! first operand
        ExpressionRef<T1> op1;
        //! second operand
        ExpressionRef<T2> op2;
        
        /**********************************************************************/
        LinearWeakDiff(const T1& x, const T2& y) :
        /* init list           */ op1(x),
        /* init list           */ op2(y)
        {/*!
          */
            
        }
        
        /**********************************************************************/
        LinearWeakDiff(T1&& x, const T2& y) :
        /* init list           */ op1(std::move(x)),
        /* init list           */ op2(y)
        {/*!
          */
            
        }
        
        /**********************************************************************/
        LinearWeakDiff(const T1& x, T2&& y) :
        /* init list           */ op1(x),
        /* init list           */ op2(std::move(y))
        {/*!
          */
            
        }
        
        /**********************************************************************/
        LinearWeakDiff(T1&& x, T2&& y) :
        /* init list           */ op1(std::move(x)),
        /* init list           */ op2(std::move(y))
        {/*!
          */
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            return op1().globalVector()-op2().globalVector();
        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1, typename T2>
    LinearWeakDiff<T1,T2> operator-(const LinearWeakExpression<T1>& op1,const LinearWeakExpression<T2>& op2)
    {
        return  LinearWeakDiff<T1,T2>(op1.derived(),op2.derived());
    }
    
    /**************************************************************************/
    template <typename T1, typename T2>
    LinearWeakDiff<T1,T2> operator-(LinearWeakExpression<T1>&& op1,const LinearWeakExpression<T2>& op2)
    {
        return  LinearWeakDiff<T1,T2>(std::move(op1.derived()),op2.derived());
    }
    
    /**************************************************************************/
    template <typename T1, typename T2>
    LinearWeakDiff<T1,T2> operator-(const LinearWeakExpression<T1>& op1,LinearWeakExpression<T2>&& op2)
    {
        return  LinearWeakDiff<T1,T2>(op1.derived(),std::move(op2.derived()));
    }
    
    /**************************************************************************/
    template <typename T1, typename T2>
    LinearWeakDiff<T1,T2> operator-(LinearWeakExpression<T1>&& op1,LinearWeakExpression<T2>&& op2)
    {
        return  LinearWeakDiff<T1,T2>(std::move(op1.derived()),std::move(op2.derived()));
    }
    
    
    
}	// close namespace
#endif
