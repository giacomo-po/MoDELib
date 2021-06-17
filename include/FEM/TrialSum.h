/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialSum_H_
#define model_TrialSum_H_

#include <type_traits> // std::is_same
#include <utility> // for std::move
#include <stdlib.h>
//#include <AreSameType.h>
//#include <TypeTraits.h>
#include <TrialBase.h>
#include <TrialExpressionBase.h>

namespace model
{
 	/**************************************************************************/
	/**************************************************************************/
	/*! \brief A class template that performs the sum of two TrialExpressionBase.
	 */
    template <typename T1,typename T2>
	class TrialSum : //public TrialBase<typename T1::TrialFunctionType>,
    /*            */ public TrialExpressionBase<TrialSum<T1,T2> >
    {
        
        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING TrialExpressionBase OF DIFFERENT TRIALFUNCTIONS.");
        static_assert(T1::rows==T2::rows,"YOU ARE SUMMING TrialExpressionBaseS WITH DIFFERENT NUMBER OF ROWS");
        typedef typename T1::TrialFunctionType _TrialFunctionType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType _ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::BaryType _BaryType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;
        constexpr static int dim=TypeTraits<_TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<_TrialFunctionType>::dofPerElement;
        
//        //! first operand
//        const T1& op1;
//        //! second operand
//        const T2& op2;
        
//        internal::ExpressionRef<T1> op1;
//        internal::ExpressionRef<T2> op2;

        
    public:
        
        constexpr static int rows=T1::rows;
        
//        /**********************************************************************/
//        TrialSum(const T1& x, const T2& y) :
//        /* base initialization */ TrialBaseType(x.trial()),
//        /* init list           */ op1(x),
//        /* init list           */ op2(y)
//        {/*!
//          */
//            
//            if(&op1().trial()!=&op2().trial())
//            {
//                std::cout<<"&op1.trial()="<<&op1().trial()<<std::endl;
//                std::cout<<"&op2.trial()="<<&op2().trial()<<std::endl;
//                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
//                //std::exit(EXIT_FAILURE);
//            }
//        }
        
//        /**********************************************************************/
//        TrialSum(T1&& x, const T2& y) :
//        /* base initialization */ TrialBaseType(x.trial()),
//        /* init list           */ op1(std::move(x)),
//        /* init list           */ op2(y)
//        {/*!
//          */
//            
//            if(&op1().trial()!=&op2().trial())
//            {
//                std::cout<<"&op1.trial()="<<&op1().trial()<<std::endl;
//                std::cout<<"&op2.trial()="<<&op2().trial()<<std::endl;
//                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
//                //std::exit(EXIT_FAILURE);
//            }
//        }
//        
//        /**********************************************************************/
//        TrialSum(const T1& x,T2&& y) :
//        /* base initialization */ TrialBaseType(x.trial()),
//        /* init list           */ op1(x),
//        /* init list           */ op2(std::move(y))
//        {/*!
//          */
//            
//            if(&op1().trial()!=&op2().trial())
//            {
//                std::cout<<"&op1.trial()="<<&op1().trial()<<std::endl;
//                std::cout<<"&op2.trial()="<<&op2().trial()<<std::endl;
//                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
//                //std::exit(EXIT_FAILURE);
//            }
//        }
//        
//        /**********************************************************************/
//        TrialSum(T1&& x, T2&& y) :
//        /* base initialization */ TrialBaseType(x.trial()),
//        /* init list           */ op1(std::move(x)),
//        /* init list           */ op2(std::move(y))
//        {/*!
//          */
//            
//            if(&op1().trial()!=&op2().trial())
//            {
//                std::cout<<"&op1.trial()="<<&op1().trial()<<std::endl;
//                std::cout<<"&op2.trial()="<<&op2().trial()<<std::endl;
//                std::cout<<"SUMMING EXPRESSIONS OF DIFFERENT TRIAL FUNCTIONS! Exiting."<<std::endl;
//                //std::exit(EXIT_FAILURE);
//            }
//        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,rows,TypeTraits<_TrialFunctionType>::dofPerElement> sfm(const _ElementType& ele,
                                                                                     const  _BaryType& bary)
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the sum of the shape-function-matrices of the operands,
          * evaluated at bary on ele.
          */
            return T1::sfm(ele,bary)+T2::sfm(ele,bary);
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,rows*dim,dofPerElement> sfmGrad(const _ElementType& ele,
                                                             const  _BaryType& bary)
        {/*!@param[in] elem the element
          * @param[in] bary the barycentric cooridinate
          *\returns the sum of gradients of the shape-function-matrices of the
          * operands, evaluated at bary on ele.
          */
            return T1::sfmGrad(ele,bary)+T2::sfmGrad(ele,bary);
        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const _ElementType& ele, const _BaryType& bary) const
//        {/*!@param[in] ele the element
//          * @param[in] bary the vector of barycentric coordinates
//          * \returns the value of the Derived expression at bary.
//          */
//            return sfm(ele,bary)*this->trial().dofs(ele);
//        }
//        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P,
//                                                const Simplex<dim,dim>* guess) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            const std::pair<bool,const _ElementType*> temp=this->trial().fe.searchWithGuess(P,guess);
//            //            return (temp.first? sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second)) : Eigen::Matrix<double,rows,1>::Zero());
//            Eigen::Matrix<double,rows,1> val(Eigen::Matrix<double,rows,1>::Zero());
//            if(temp.first)
//            {
//                val=sfm(*(temp.second),temp.second->simplex.pos2bary(P))*this->trial().dofs(*(temp.second));
//            }
//            return val;
//        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(op1.trial.fe.mesh.begin()->second));
//        }
        
    };
    
    
    /**************************************************************************/
    template <typename T1, typename T2>
    TrialSum<T1,T2> operator+(const TrialExpressionBase<T1>& op1,
                              const TrialExpressionBase<T2>& op2)
    {
        return TrialSum<T1,T2>(op1,op2);
    }
    
    template <typename T1, typename T2>
    TrialSum<T1,T2> operator+(TrialExpressionBase<T1>&& op1,
                              const TrialExpressionBase<T2>& op2)
    {
        return TrialSum<T1,T2>(std::move(op1),op2);
    }
    
    template <typename T1, typename T2>
    TrialSum<T1,T2> operator+(const TrialExpressionBase<T1>& op1,
                            TrialExpressionBase<T2>&& op2)
    {
        return TrialSum<T1,T2>(op1,std::move(op2));
    }
    
    template <typename T1, typename T2>
    TrialSum<T1,T2> operator+(TrialExpressionBase<T1>&& op1,
                            TrialExpressionBase<T2>&& op2)
    {
        return TrialSum<T1,T2>(std::move(op1),std::move(op2));
    }
    
}	// close namespace
#endif
