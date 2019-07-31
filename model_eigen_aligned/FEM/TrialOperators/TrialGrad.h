/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialGrad_H_
#define model_TrialGrad_H_

#include <utility> // for std::move
#include <type_traits> //static_assert
//#include <memory> // shared_ptr
#include <TypeTraits.h>
#include <Simplex.h>
#include <TrialBase.h>
#include <TrialExpressionBase.h>
#include <ExpressionRef.h>

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T>
    class TrialGrad : //public TrialBase<typename T::TrialFunctionType>,
    /*             */ public TrialExpressionBase<TrialGrad<T> >
    {
        
        
        
        //const T& op;
//        internal::ExpressionRef<T> op;

    public:
        
        typedef typename T::TrialFunctionType TrialFunctionType;
//        typedef TrialBase<TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=TrialFunctionType::dim*T::rows;
        constexpr static int cols=T::cols;
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;
//        const std::shared_ptr<T> op;

//        const T& op;
        ExpressionRef<T> op;

        /**********************************************************************/
        TrialGrad(const T& exp) :
        /* init */ op(exp)
        {
//            std::cout<<"TrialGrad constructor 1"<<std::endl;
        }
        
        /**********************************************************************/
        TrialGrad(T&& exp) :
        /* init */ op(std::move(exp))
        {
            //            std::cout<<"TrialGrad constructor 1"<<std::endl;
        }
        
//        /**********************************************************************/
//        TrialGrad(const TrialExpressionBase<T>& exp)
////        /* base init */ TrialBaseType(x.derived().trial()),
////        /* init list */ op(x.derived())
//        {
//            std::cout<<"TrialGrad constructor 1"<<std::endl;
//        }
        
        //        /**********************************************************************/
//        TrialGrad(const TrialExpressionBase<T>& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
//        /* init list */ op(x.derived())
//        {
//            std::cout<<"TrialGrad constructor 1"<<std::endl;
//        }
        
//        /**********************************************************************/
//        TrialGrad(TrialExpressionBase<T>&& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
//        /* init list */ op(std::move(x.derived()))
//        {
//            std::cout<<"TrialGrad constructor 2"<<std::endl;
//        }
        
        /**********************************************************************/
        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                    const BaryType& bary) const
        {
            return op().sfmGrad(ele,bary);
        }
        
        //        /**********************************************************************/
        //        void sfmGrad(const typename TF::ElementType& ele, const Eigen::Matrix<double,TF::dim+1,1>& bary) const
        //        {
        //            static_assert(0,"SECOND GRADIENTS ARE NOT SUPPORTED YET.");
        //        }
        
//        /**********************************************************************/
//        EvalMatrixType operator()(const ElementType& ele,
//                                  const BaryType& bary) const
//        {
//            return TrialBase<TrialFunctionType>::eval(*this,ele,bary);
//        }
//
//        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P,
//                                                const Simplex<dim,dim>* guess) const
//        {
//            return TrialBase<TrialFunctionType>::eval(*this,P,guess);
//        }
//        
//        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {
//            return TrialBase<TrialFunctionType>::eval(*this,P);
//        }

    };
    

    
    
    /**************************************************************************/
    template <typename Derived>
    TrialGrad<Derived> grad(const TrialExpressionBase<Derived>& x)
    {
        return TrialGrad<Derived>(x.derived());
    }
    
    /**************************************************************************/
    template <typename Derived>
    TrialGrad<Derived> grad(TrialExpressionBase<Derived>&& x)
    {
        return TrialGrad<Derived>(std::move(x.derived()));
    }

    
    
    
//    /**************************************************************************/
//    template <typename T>
//    TrialGrad<T> grad(TrialExpressionBase<T>&& x)
//    {
//        return TrialGrad<T>(std::move(x));
//    }
//    
}
#endif


//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(op.trial.fe.mesh.begin()->second));
//        }
