/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDef_H_
#define model_TrialDef_H_

#include <utility> // for std::move
#include <type_traits> //static_assert
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
    struct TrialDef : public TrialExpressionBase<TrialDef<T> >
    {
        
        typedef typename T::TrialFunctionType TrialFunctionType;
//        typedef TrialBase<TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=(TrialFunctionType::dim*(TrialFunctionType::dim+1))/2;
        constexpr static int cols=T::cols;
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;


        static_assert(TypeTraits<TrialFunctionType>::dim==TypeTraits<TrialFunctionType>::nComponents,"SYMMETRIC GRADIENT (DEF) CAN ONLY BE COMPUTED IF nComponents==dim.");
        
//        const T& op;
        ExpressionRef<T> op;

        
    //public:
        

        /**********************************************************************/
        TrialDef(const T& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(x)
        {
            
        }
        
        /**********************************************************************/
        TrialDef(T&& x) :
//        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(std::move(x))
        {
            
        }
        
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,dofPerElement> sfm(const ElementType& ele,
                                                     const BaryType& bary) const
        {
            return op().sfmDef(ele,bary);
        }
        
        //        /**********************************************************************/
        //        void sfmGrad(const typename TF::ElementType& ele, const Eigen::Matrix<double,TF::dim+1,1>& bary) const
        //        {
        //            static_assert(0,"SECOND GRADIENTS ARE NOT SUPPORTED YET.");
        //        }

    };
    
    
    /**************************************************************************/
    template <typename T>
    TrialDef<T> def(const TrialExpressionBase<T>& x)
    {
        return TrialDef<T>(x.derived());
    }
    
    /**************************************************************************/
    template <typename T>
    TrialDef<T> def(TrialExpressionBase<T>&& x)
    {
        return TrialDef<T>(std::move(x.derived()));
    }
    
}	// close namespace
#endif


//        /**********************************************************************/
//        Eigen::Matrix<double,rows,1> operator()(const Eigen::Matrix<double,dim,1>& P) const
//        {/*!@param[in] P the position vector
//          * @param[in] guess the Simplex where the search starts
//          * \returns the value of the Derived expression at P.
//          */
//            return operator()(P,&(op.trial.fe.mesh.begin()->second));
//
//        }
