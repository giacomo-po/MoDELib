/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TestExpression_H_
#define model_TestExpression_H_

#include <utility> // for std::move
#include <model/FEM/TrialOperators/ExpressionRef.h>


namespace model
{
    
//    template<typename Derived>
//    struct TestExpression
//    {
//        /**********************************************************************/
//        Derived& derived()
//        {/*!\returns A reference to the Derived object
//          */
//            return *static_cast<Derived*>(this);
//        }
//        
//        /**********************************************************************/
//        Derived* p_derived()
//        {/*!\returns A pointer to the Derived object
//          */
//            return  static_cast<Derived*>(this);
//        }
//        
//        /**********************************************************************/
//        const Derived& derived() const
//        {/*! A const reference to the Derived object
//          */
//            return *static_cast<const Derived*>(this);
//        }
//        
//        /**********************************************************************/
//        const Derived* p_derived() const
//        {/*!\returns A pointer to the const Derived object
//          */
//            return  static_cast<const Derived*>(this);
//        }
//        
//        /**********************************************************************/
//        operator const Derived&() const
//        {/*!\returns the const Derived object (cast operator)
//          */
//            return derived();
//        }
//    };
    

    /**************************************************************************/
    /**************************************************************************/
    template<typename TrialExpressionType>
    struct TestExpression : public internal::ExpressionRef<TrialExpressionType>
    {
        
//        typedef TestExpression<TrialExpressionType> TestExpressionType;
//        typedef typename TrialExpressionType::ShapeFunctionMatrixType     ShapeFunctionMatrixType;
//        typedef typename TrialExpressionType::BaryType BaryType;
//        typedef typename TrialExpressionType::ElementType ElementType;
        
        //const TrialExpressionType& trialExp;
//        internal::ExpressionRef<TrialExpressionType> trialExp;

        typedef internal::ExpressionRef<TrialExpressionType> BaseType;
        
        /**********************************************************************/
        TestExpression(const TrialExpressionType& trialExp_in) :
        /* init list */ BaseType(trialExp_in)
        {
            
        }
        
        /**********************************************************************/
        TestExpression(TrialExpressionType&& trialExp_in) :
        /* init list */ BaseType(std::move(trialExp_in))
        {
            
        }
        
//        /**********************************************************************/
//        ShapeFunctionMatrixType sfm(const ElementType& ele, const BaryType& bary) const
//        {/*!\param[in] ele a finite element
//          * \param[in] bary the vector of barycentric coordinates
//          * \returns The matrix of shape functions for the element, evaluated at
//          * bary.
//          */
//            return trialExp().sfm(ele,bary);
//        }
//        
//        /**********************************************************************/
//        size_t gSize() const
//        {
//            return trialExp().gSize();
//        }
        
    };
    
}	// close namespace
#endif


/**************************************************************************/
/**************************************************************************/
//	template<typename TrialExpressionType>
//	struct TestExpression : public TrialExpressionType
//    {
//
//        typedef TestExpression<TrialExpressionType> TestExpressionType ;
//
//        /**********************************************************************/
//        TestExpression(const TrialExpressionType& trialExp) :
//        /* init list */ TrialExpressionType(trialExp)
//        {
//
//        }
//
//    };
