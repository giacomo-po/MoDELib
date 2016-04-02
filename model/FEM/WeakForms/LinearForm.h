/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearForm_H_
#define model_LinearForm_H_

#include <Eigen/Dense>
#include <model/FEM/TrialOperators/EvalExpression.h>
#include <model/FEM/TrialOperators/TestExpression.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DerivedTest,typename DerivedEval>
    class LinearForm : public LinearWeakExpression<LinearForm<DerivedTest,DerivedEval>>
    {
        
        static_assert(DerivedTest::rows==DerivedEval::rows,"YOU ARE CREATING A LinearForm BETWEEN EXPRESSIONS WITH DIFFERENT NUMBER OF ROWS");
        
    public:
        
        typedef typename DerivedTest::TrialFunctionType  TrialFunctionType;
        constexpr static int evalCols=DerivedEval::cols;
        
        
    private:
        typedef LinearForm<DerivedTest,DerivedEval> LinearFormType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        
        constexpr static int dim=TypeTraits< TrialFunctionType>::dim;
        constexpr static int dofPerNode=TypeTraits< TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits< TrialFunctionType>::dofPerElement;
        
        
        
        
        
    public:
        
        const DerivedTest  testExp;
        const DerivedEval  evalExp;
        
        //        /**********************************************************************/
        //        LinearForm(const TestExpression<DerivedTest>& testE, const EvalExpression<DerivedEval>& evalE) :
        //        /* init list */ testExp(testE.trial()),
        //        /* init list */ evalExp(evalE.derived())
        //        {
        //
        //        }
        
        /**********************************************************************/
        LinearForm(const TestExpression<DerivedTest>& test, const EvalExpression<DerivedEval>& eval) :
        //        /* init list */ testExp(testE.trial()),
        /* init list */ testExp(test),
        /* init list */ evalExp(eval.wrappedExp)
        {
            
        }
        
        /**********************************************************************/
        // Do not allow construction from temporary evalE
        //        LinearForm(const TestExpression<DerivedTest>& testE, EvalExpression<DerivedEval>&& evalE) = delete ;
        
        //        LinearForm(const TestExpression<DerivedTest>& testE, const DerivedEval&& evalE) = delete ;
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,evalCols> operator()(const ElementType& ele, const BaryType& bary) const
        {
            return testExp.sfm(ele,bary).transpose()*evalExp(ele,bary);
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operators
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (const TestExpression<T1>& testExp, const EvalExpression<T2>& evalExp)
    {
        return LinearForm<T1,T2>(testExp,evalExp);
    }
    
    template <typename T1>
    LinearForm<T1,Constant<double,1,1> > operator, (const TestExpression<T1>& testE, const double& c)
    {
        return LinearForm<T1,Constant<double,1,1> >(testE,make_constant(c).eval());
    }
    
//    template <typename T1, int rows, int cols>
//    LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testE,
//                                                                                   const Eigen::Matrix<double,rows,cols>& c)
//    {
//        return LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testE,make_constant(c).eval());
//    }
    
    template <typename T1, int rows, int cols>
    LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testE,
                                                                                   const Constant<Eigen::Matrix<double,rows,cols>,rows,cols>& c)
    {
        return LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testE,c.eval());
    }
    
}	// close namespace
#endif


