/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearForm_H_
#define model_BilinearForm_H_

#include <type_traits> // std::is_same
#include <TrialExpressionBase.h>
//#include <AreSameType.h>
#include <TerminalColors.h>
#include <ExpressionRef.h>

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
        
        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE CREATING A BilinearForm OF DIFFERENT TrialFunction TYPES.");
        static_assert(T1::rows==T2::rows,"YOU ARE CREATING A BilinearForm BETWEEN A TrialExpression AND A TestExpression WITH DIFFERENT NUMBER OF ROWS");

        typedef typename TestExpression<T1>::TestExpressionType TestExpressionType;
        typedef T2 TrialExpressionType;
        
        typedef typename T2::TrialFunctionType TrialFunctionType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::FiniteElementType FiniteElementType;
        
        
//        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
//        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        
        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;

        
//        const TestExpressionType&   testExpr;
//        const TrialExpressionType& trialExpr;
//        const size_t gSize;
        
        ExpressionRef<TestExpression<T1>> testExpr;
        ExpressionRef<TrialExpressionType> trialExpr;
        
        /**********************************************************************/
        BilinearForm(const TestExpression<T1>& testE, const T2& trialE) :
        /* init list */ testExpr(testE), // cast testE to its base T2 type
        /* init list */ trialExpr(trialE) // cast trialE to its derived T1 type
//        /* init list */ gSize(trialExpr.nodeSize()*dofPerNode)
        {
            
//            std::cout<<greenColor<<"Creating BilinearForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        BilinearForm(TestExpression<T1>&& testE, const T2& trialE) :
        /* init list */ testExpr(std::move(testE)), // cast testE to its base T2 type
        /* init list */ trialExpr(trialE) // cast trialE to its derived T1 type
        //        /* init list */ gSize(trialExpr.nodeSize()*dofPerNode)
        {
            
            //            std::cout<<greenColor<<"Creating BilinearForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        BilinearForm(const TestExpression<T1>& testE, T2&& trialE) :
        /* init list */ testExpr(testE), // cast testE to its base T2 type
        /* init list */ trialExpr(std::move(trialE)) // cast trialE to its derived T1 type
        //        /* init list */ gSize(trialExpr.nodeSize()*dofPerNode)
        {
            
            //            std::cout<<greenColor<<"Creating BilinearForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        BilinearForm(TestExpression<T1>&& testE, T2&& trialE) :
        /* init list */ testExpr(std::move(testE)), // cast testE to its base T2 type
        /* init list */ trialExpr(std::move(trialE)) // cast trialE to its derived T1 type
        //        /* init list */ gSize(trialExpr.nodeSize()*dofPerNode)
        {
            
            //            std::cout<<greenColor<<"Creating BilinearForm: gSize="<<gSize<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
//        template <typename AbscissaType>
        ElementMatrixType operator()(const ElementType& ele, const BaryType& bary) const
        {
//            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
            return testExpr().sfm(ele,bary).transpose()*trialExpr().sfm(ele,bary);
        }
        
    };
    
    
    template <typename T1,typename T2>
    BilinearForm<T1,T2> operator, (const TestExpression<T1>& testExp, const TrialExpressionBase<T2>& trialExp)
    {
        return BilinearForm<T1,T2>(testExp,trialExp.derived());
    }

    template <typename T1,typename T2>
    BilinearForm<T1,T2> operator, (TestExpression<T1>&& testExp, const TrialExpressionBase<T2>& trialExp)
    {
        return BilinearForm<T1,T2>(std::move(testExp),trialExp.derived());
    }
    
    template <typename T1,typename T2>
    BilinearForm<T1,T2> operator, (const TestExpression<T1>& testExp, TrialExpressionBase<T2>&& trialExp)
    {
        return BilinearForm<T1,T2>(testExp,std::move(trialExp.derived()));
    }
    
    template <typename T1,typename T2>
    BilinearForm<T1,T2> operator, (TestExpression<T1>&& testExp, TrialExpressionBase<T2>&& trialExp)
    {
        return BilinearForm<T1,T2>(std::move(testExp),std::move(trialExp.derived()));
    }

    
}	// close namespace
#endif

