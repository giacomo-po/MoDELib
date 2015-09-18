/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakSum_H_
#define model_LinearWeakSum_H_

#include <model/Utilities/AreSameType.h>
#include <model/FEM/WeakForms/LinearWeakExpression.h>

namespace model
{
 	/**************************************************************************/
	/**************************************************************************/
	/*! \brief A class template that performs the sum of two TrialExpressionBase.
	 */
    template <typename T1,typename T2>
	struct  LinearWeakSum : public LinearWeakExpression<LinearWeakSum<T1,T2> >
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING TrialExpressionBase OF DIFFERENT TRIALFUNCTIONS.");
        typedef typename T1::TrialFunctionType TrialFunctionType;
        
        //! first operand
        const T1 op1; // expression could be a temporary, so copy by value
        //! second operand
        const T2 op2; // expression could be a temporary, so copy by value
        
        /**********************************************************************/
        LinearWeakSum(const T1& x, const T2& y) :
        /* init list           */ op1(x),
        /* init list           */ op2(y)
        {/*!
          */
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            return op1.globalVector()+op2.globalVector();
        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1, typename T2>
    LinearWeakSum<T1,T2> operator+(const LinearWeakExpression<T1>& op1,const LinearWeakExpression<T2>& op2)
    {
        return  LinearWeakSum<T1,T2>(op1,op2);
    }
    
}	// close namespace
#endif
