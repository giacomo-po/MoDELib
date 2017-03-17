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
    template <typename TestType,typename EvalType>
    class LinearForm : public LinearWeakExpression<LinearForm<TestType,EvalType>>
    {
        
        static_assert(TestType::rows==EvalType::rows,"YOU ARE CREATING A LinearForm BETWEEN EXPRESSIONS WITH DIFFERENT NUMBER OF ROWS");
        
    public:
        
        typedef typename TestType::TrialFunctionType  TrialFunctionType;
        constexpr static int evalCols=EvalType::cols;
        
        
    private:
        typedef LinearForm<TestType,EvalType> LinearFormType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        
        constexpr static int dim=TypeTraits< TrialFunctionType>::dim;
        constexpr static int dofPerNode=TypeTraits< TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits< TrialFunctionType>::dofPerElement;
        
        
        
        
        
    public:
        
//        const TestType&  testExp;
//        const TestType&  testExp;
//        const EvalType&  evalExp;
        
        internal::ExpressionRef<TestType> testExp;
        internal::ExpressionRef<EvalType> evalExp;

        
        //        /**********************************************************************/
        //        LinearForm(const TestExpression<TestType>& testE, const EvalExpression<EvalType>& evalE) :
        //        /* init list */ testExp(testE.trial()),
        //        /* init list */ evalExp(evalE.derived())
        //        {
        //
        //        }
        
        /**********************************************************************/
        LinearForm(const TestExpression<TestType>& test, const EvalExpression<EvalType>& eval) :
        //        /* init list */ testExp(testE.trial()),
        /* init list */ testExp(test.trialExp()),
        /* init list */ evalExp(eval.trialExp())
        {
            
        }
        
        /**********************************************************************/
        LinearForm(TestExpression<TestType>&& test, const EvalExpression<EvalType>& eval) :
        //        /* init list */ testExp(testE.trial()),
        /* init list */ testExp(std::move(test.trialExp())),
        /* init list */ evalExp(eval.trialExp())
        {
            
        }
        
        /**********************************************************************/
        LinearForm(const TestExpression<TestType>& test, EvalExpression<EvalType>&& eval) :
        //        /* init list */ testExp(testE.trial()),
        /* init list */ testExp(test.trialExp()),
        /* init list */ evalExp(std::move(eval.trialExp()))
        {
            
        }
        
        /**********************************************************************/
        LinearForm(TestExpression<TestType>&& test, EvalExpression<EvalType>&& eval) :
        //        /* init list */ testExp(testE.trial()),
        /* init list */ testExp(std::move(test.trialExp())),
        /* init list */ evalExp(std::move(eval.trialExp()))
        {
            
        }
        
        /**********************************************************************/
        // Do not allow construction from temporary evalE
        //        LinearForm(const TestExpression<TestType>& testE, EvalExpression<EvalType>&& evalE) = delete ;
        
        //        LinearForm(const TestExpression<TestType>& testE, const EvalType&& evalE) = delete ;
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dofPerElement,evalCols> operator()(const ElementType& ele, const BaryType& bary) const
        {
            return testExp().sfm(ele,bary).transpose()*evalExp()(ele,bary);
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operators
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (const TestExpression<T1>& testExp,
                                 const EvalExpression<T2>& evalExp)
    {
                std::cout<<"LinearForm operator, 1"<<std::endl;
        return LinearForm<T1,T2>(testExp,evalExp);
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (TestExpression<T1>&& testExp,
                                 const EvalExpression<T2>& evalExp)
    {
                std::cout<<"LinearForm operator, 2"<<std::endl;
        return LinearForm<T1,T2>(std::move(testExp),evalExp);
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (const TestExpression<T1>& testExp,
                                EvalExpression<T2>&& evalExp)
    {
                std::cout<<"LinearForm operator, 3"<<std::endl;
        return LinearForm<T1,T2>(testExp,std::move(evalExp));
    }
    
    template <typename T1,typename T2>
    LinearForm<T1,T2> operator, (TestExpression<T1>&& testExp,
                                 EvalExpression<T2>&& evalExp)
    {
        std::cout<<"LinearForm operator, 4"<<std::endl;
        return LinearForm<T1,T2>(std::move(testExp),std::move(evalExp));
    }
    
//    template <typename T1>
//    LinearForm<T1,Constant<double,1,1> > operator, (const TestExpression<T1>& testExp,
//                                                    const double& c)
//    {
////        return LinearForm<T1,Constant<double,1,1> >(testExp,make_constant(c).eval());
//        return LinearForm<T1,Constant<double,1,1> >(testExp,make_constant(c));
//
//    }
    
//    template <typename T1, int rows, int cols>
//    LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testE,
//                                                                                   const Eigen::Matrix<double,rows,cols>& c)
//    {
//        return LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testE,make_constant(c).eval());
//    }
    
//    template <typename T1, int rows, int cols>
//    LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testExp,
//                                                                                   const Constant<Eigen::Matrix<double,rows,cols>,rows,cols>& c)
//    {
////        return LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testExp,c.eval());
//        return LinearForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testExp,eval(c));
//
//    }
    
}	// close namespace
#endif


