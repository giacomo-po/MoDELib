/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDef_H_
#define model_TrialDef_H_

#include <type_traits> //static_assert
#include <model/Utilities/TypeTraits.h>
#include <model/FEM/TrialOperators/TrialBase.h>
#include <model/FEM/TrialOperators/TrialExpressionBase.h>

namespace model
{
    
	
    /**************************************************************************/
	/**************************************************************************/
    template <typename T>
    class TrialDef : public TrialBase<typename T::TrialFunctionType>,
    /*            */ public TrialExpressionBase<TrialDef<T> >
    {
        
        typedef typename T::TrialFunctionType _TrialFunctionType;
        typedef TrialBase<_TrialFunctionType> TrialBaseType;

        static_assert(TypeTraits<_TrialFunctionType>::dim==TypeTraits<_TrialFunctionType>::nComponents,"SYMMETRIC GRADIENT (DEF) CAN ONLY BE COMPUTED IF nComponents==dim.");
        
        const T op; // Input expression could be a temporary, so copy by value
        
    public:
        
        constexpr static int rows=(_TrialFunctionType::dim*(_TrialFunctionType::dim+1))/2;

        /**********************************************************************/
        TrialDef(const TrialExpressionBase<T>& x) :
        /* base init */ TrialBaseType(x.derived().trial()),
        /* init list */ op(x.derived())
        {
            
        }
        
        
        /**********************************************************************/
        Eigen::Matrix<double,rows,TypeTraits<_TrialFunctionType>::dofPerElement> sfm(const typename TypeTraits<_TrialFunctionType>::ElementType& ele,
                                                                                     const typename TypeTraits<_TrialFunctionType>::BaryType& bary) const
        {
            return op.sfmDef(ele,bary);
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
        return TrialDef<T>(x);
    }
    
}	// close namespace
#endif
