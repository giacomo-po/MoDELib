/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TestExpression_H_
#define model_TestExpression_H_

//#include <model/FEM/WeakForms/LinearWeakForm.h>
//#include <model/FEM/WeakForms/BilinearWeakForm.h>



namespace model
{
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<typename TrialExpressionType>
	struct TestExpression : public TrialExpressionType
    {
        
        
        typedef TestExpression<TrialExpressionType> TestExpressionType ;
//        typedef typename TrialExpressionType::TrialFunctionType TrialFunctionType;
//        //        constexpr static int rows=TrialExpressionType::rows;
//        
//        
//        const TrialExpressionType trialExp;
//        
        /**********************************************************************/
        TestExpression(const TrialExpressionType& trialExp) :
        /* init list */ TrialExpressionType(trialExp)
        {
            
        }
//
//        /**********************************************************************/
//        const TrialExpressionType& trial() const
//        {
//            return trialExp;
//        }
        
//        template <typename T>
//        BilinearWeakForm<TrialExpressionType,T> operator, (const TrialExpressionBase<T>& trialE) const
//        {
//            return BilinearWeakForm<TrialExpressionType,T>(*this,trialE);
//        }
        
        
    };
    
//    template <typename T1,typename T2>
//    LinearWeakForm<T1,T2> operator, (const TestExpression<T1>& testE, const EvalExpression<T2>& evalE)
//    {
//        return LinearWeakForm<T1,T2>(testE,evalE);
//    }
//    
//    template <typename T1>
//    LinearWeakForm<T1,Constant<double,1,1> > operator, (const TestExpression<T1>& testE, const double& c)
//    {
//        return LinearWeakForm<T1,Constant<double,1,1> >(testE,make_constant(c));
//    }
//    
//    template <typename T1, int rows, int cols>
//    LinearWeakForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testE, const Eigen::Matrix<double,rows,cols>& c)
//    {
//        return LinearWeakForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testE,make_constant(c));
//    }
    
    
//    template <typename T1,typename T2>
//    BilinearWeakForm<T1,T2> operator, (const TestExpression<T2>& testE, const TrialExpressionBase<T1>& trialE)
//    {
//        return BilinearWeakForm<T1,T2>(testE,trialE);
//    }
	
    

    
}	// close namespace
#endif



//    /**************************************************************************/
//	/**************************************************************************/
//	template<typename _TrialExpressionType>
//	struct TestExpression
//    {
//
//
//        typedef _TrialExpressionType TrialExpressionType;
//        typedef typename TrialExpressionType::TrialFunctionType TrialFunctionType;
////        constexpr static int rows=TrialExpressionType::rows;
//
//
//        const TrialExpressionType trialExp;
//
//        /**********************************************************************/
//        TestExpression(const TrialExpressionType& trialE_in) :
//        /* init list */ trialExp(trialE_in)
//        {
//
//        }
//
//        /**********************************************************************/
//        const TrialExpressionType& trial() const
//        {
//            return trialExp;
//        }
//
//
//    };
