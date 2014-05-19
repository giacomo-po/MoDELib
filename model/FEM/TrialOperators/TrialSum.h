/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialSum_H_
#define model_TrialSum_H_

#include <stdlib.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TypeTraits.h>
#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>

namespace model
{
 	/**************************************************************************/
	/**************************************************************************/
	/*! \brief A class template that performs the sum of two TrialExpressionBase.
	 */
    template <typename T1,typename T2>
	class TrialSum : public TrialBase<typename T1::TrialFunctionType>,
    /*            */ public TrialExpressionBase<TrialSum<T1,T2> >
    {
        
        static_assert(AreSameType<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING TrialExpressionBase OF DIFFERENT TRIALFUNCTIONS.");
        static_assert(T1::rows==T2::rows,"YOU ARE SUMMING TrialExpressionBaseS WITH DIFFERENT NUMBER OF ROWS");
        typedef typename T1::TrialFunctionType _TrialFunctionType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;

        //! first operand
        const T1 op1; // expression could be a temporary, so copy by value
        //! second operand
        const T2 op2; // expression could be a temporary, so copy by value
        
    public:
        
        constexpr static int rows=T1::rows;
        
        
//        /**********************************************************************/
//        TrialSum(const TrialExpressionBase<T1>& x, const TrialExpressionBase<T2>& y) :
//        /* base initialization */ TrialBaseType(x.derived().trial()),
//        /* init list           */ op1(x.derived()),
//        /* init list           */ op2(y.derived())
//        {/*!
//          */
//            
//            if(&op1.trial()!=&op2.trial())
//            {
//                std::cout<<"&op1.trial()="<<&op1.trial()<<std::endl;
//                std::cout<<"&op2.trial()="<<&op2.trial()<<std::endl;
//                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
//                //std::exit(EXIT_FAILURE);
//            }
//        }

        /**********************************************************************/
        TrialSum(const T1& x, const T2& y) :
        /* base initialization */ TrialBaseType(x.trial()),
        /* init list           */ op1(x),
        /* init list           */ op2(y)
        {/*!
          */
            
            if(&op1.trial()!=&op2.trial())
            {
                std::cout<<"&op1.trial()="<<&op1.trial()<<std::endl;
                std::cout<<"&op2.trial()="<<&op2.trial()<<std::endl;
                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
                //std::exit(EXIT_FAILURE);
            }
        }

        
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,TypeTraits<_TrialFunctionType>::dofPerElement> sfm(const typename TypeTraits<_TrialFunctionType>::ElementType& ele,
                                                                                     const typename TypeTraits<_TrialFunctionType>::BaryType& bary) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the sum of the shape-function-matrices of the operands,
          * evaluated at bary on ele.
          */
            return op1.sfm(ele,bary)+op2.sfm(ele,bary);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,rows*TypeTraits<_TrialFunctionType>::dim,TypeTraits<_TrialFunctionType>::dofPerElement> sfmGrad(const typename TypeTraits<_TrialFunctionType>::ElementType& ele,
                                                                                                                             const typename TypeTraits<_TrialFunctionType>::BaryType& bary) const
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the sum of gradients of the shape-function-matrices of the
          * operands, evaluated at bary on ele.
          */
            return op1.sfmGrad(ele,bary)+op2.sfmGrad(ele,bary);
        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1, typename T2>
    TrialSum<T1,T2> operator+(const TrialExpressionBase<T1>& op1,const TrialExpressionBase<T2>& op2)
    {
        return TrialSum<T1,T2>(op1,op2);
    }
    
}	// close namespace
#endif


//    /**************************************************************************/
//    template <int N, typename FE, typename T2>
//    TrialSum<TrialExpression<TrialFunction<N,FE> >,T2> operator+(const TrialFunction<N,FE>& op1,const TrialExpressionBase<T2>& op2)
//    {
//        return TrialSum<TrialExpression<TrialFunction<N,FE> >,T2>(TrialExpression<TrialFunction<N,FE> >(op1),op2);
//    }

//    /**************************************************************************/
//    template <int N, typename FE, typename T2>
//    TrialSum<TrialExpression<TrialFunction<N,FE> >,T2> operator+(const TrialExpressionBase<T2>& op2,const TrialFunction<N,FE>& op1)
//    {
//        return TrialSum<TrialExpression<TrialFunction<N,FE> >,T2>(TrialExpression<TrialFunction<N,FE> >(op1),op2);
//    }

//    /**************************************************************************/
//    template <int N1, typename FE1, int N2, typename FE2>
//    TrialSum<TrialExpression<TrialFunction<N1,FE1> >,TrialExpression<TrialFunction<N2,FE2> > > operator+(const TrialFunction<N1,FE1>& op1, const TrialFunction<N2,FE2>& op2)
//    {
//        return TrialSum<TrialExpression<TrialFunction<N1,FE1> >,TrialExpression<TrialFunction<N2,FE2> > >(TrialExpression<TrialFunction<N1,FE1> >(op1),TrialExpression<TrialFunction<N2,FE2> >(op2));
//    }
