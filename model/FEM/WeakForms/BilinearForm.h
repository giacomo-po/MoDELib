/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearForm_H_
#define model_BilinearForm_H_

#include <model/FEM/TrialOperators/TrialExpressionBase.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TerminalColors.h>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    /*! A class template representing the product between a TestExpression and
     * a TrialExpression.
     */
    template <typename T1,typename T2>
	struct BilinearForm
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE CREATING A BilinearForm OF DIFFERENT TrialFunction TYPES.");
        static_assert(T1::rows==T2::rows,"YOU ARE CREATING A BilinearForm BETWEEN A TrialExpression AND A TestExpression WITH DIFFERENT NUMBER OF ROWS");

        typedef typename TestExpression<T1>::TestExpressionType TestExpressionType;
        typedef T2 TrialExpressionType;
        typedef typename T2::TrialFunctionType TrialFunctionType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::FiniteElementType FiniteElementType;
        
        
        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        
        const TestExpressionType testExpr;
        const T2 trialExpr;
        const size_t gSize;
        

        /**********************************************************************/
        BilinearForm(const TestExpression<T1>& testE, const TrialExpressionBase<T2>& trialE) :
        /* init list */ testExpr(testE), // cast testE to its base T2 type
        /* init list */ trialExpr(trialE.derived()), // cast trialE to its derived T1 type
        /* init list */ gSize(trialExpr.nodeSize()*dofPerNode)
        {
            
//            std::cout<<greenColor<<"Creating BilinearForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
    };
    
    
    template <typename T1,typename T2>
    BilinearForm<T1,T2> operator, (const TestExpression<T1>& testE, const TrialExpressionBase<T2>& trialE)
    {
        return BilinearForm<T1,T2>(testE,trialE);
    }
    
}	// close namespace
#endif

