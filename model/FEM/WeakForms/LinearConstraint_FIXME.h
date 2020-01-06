/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearConstraint_H_
#define model_LinearConstraint_H_

#include <Eigen/Dense>
#include <EvalExpression.h>
#include <EvalFunction.h>
#include <TrialExpression.h>
#include <ExpressionRef.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename EvalType,typename TrialType>
    struct LinearConstraint
    {
        
//        static_assert(TestType::rows==EvalType::rows,"YOU ARE CREATING A LinearForm BETWEEN EXPRESSIONS WITH DIFFERENT NUMBER OF ROWS");
        
        //    public:
        
        typedef typename TrialType::TrialFunctionType  TrialFunctionType;
        constexpr static int evalCols=EvalType::cols;
        
        
        //    private:
//        typedef LinearForm<TestType,EvalType> LinearFormType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        
        constexpr static int dim=TypeTraits< TrialFunctionType>::dim;
        constexpr static int dofPerNode=TypeTraits< TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits< TrialFunctionType>::dofPerElement;
        
        ExpressionRef<EvalType>   evalExp;
        ExpressionRef<TrialType> trialExp;
        
        /**********************************************************************/
        LinearConstraint(const EvalFunction<EvalType>& eval, const EvalType& eval) :
        /* init list */ testExp(test),
        /* init list */ evalExp(eval)
        {
            
        }
        
        /**********************************************************************/
        LinearForm(TestExpression<TestType>&& test, const EvalType& eval) :
        /* init list */ testExp(std::move(test)),
        /* init list */ evalExp(eval)
        {
            
        }
        
        /**********************************************************************/
        LinearForm(const TestExpression<TestType>& test, EvalType&& eval) :
        /* init list */ testExp(test),
        /* init list */ evalExp(std::move(eval))
        {
            
        }
        
        /**********************************************************************/
        LinearForm(TestExpression<TestType>&& test, EvalType&& eval) :
        /* init list */ testExp(std::move(test)),
        /* init list */ evalExp(std::move(eval))
        {
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,evalCols> operator()(const ElementType& ele,
                                                                const BaryType& bary) const
        {
            return testExp().sfm(ele,bary).transpose()*evalExp()(ele,bary);
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operators
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (const TestExpression<T1>& testExp,
                                 const   EvalFunction<T2>& evalFunc)
    {
        return LinearForm<T1,T2>(testExp,evalFunc.derived());
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (/**/ TestExpression<T1>&& testExp,
                                 const  EvalFunction<T2>& evalFunc)
    {
        return LinearForm<T1,T2>(std::move(testExp),evalFunc.derived());
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (const TestExpression<T1>&  testExp,
                                 /*   */ EvalFunction<T2>&& evalFunc)
    {
        return LinearForm<T1,T2>(testExp,std::move(evalFunc.derived()));
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (TestExpression<T1>&& testExp,
                                 EvalFunction<T2>&& evalFunc)
    {
        return LinearForm<T1,T2>(std::move(testExp),std::move(evalFunc.derived()));
    }
    
}
#endif


