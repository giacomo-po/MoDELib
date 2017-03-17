/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EvalExpression_H_
#define model_EvalExpression_H_

#include <utility> // for std::move
#include <model/FEM/TrialOperators/ExpressionRef.h>


namespace model
{
    
//    /**************************************************************************/
//    /**************************************************************************/
//    template<typename Derived>
//    struct EvalExpression
//    {
//        
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
//        
//    };
    
    //    /**************************************************************************/
	/**************************************************************************/
	template<typename TrialExpressionType>
    struct EvalExpression
    {
        
//        const T& wrappedExp;
        internal::ExpressionRef<TrialExpressionType> trialExp;
        
//        EvalExpression(const EvalExpression<T>& other) : wrappedExp(other.wrappedExp)
//        {
//            std::cout<<"EvalExpression CopyConstructor"<<std::endl;
//        }
        
        EvalExpression(const TrialExpressionType& exp) : trialExp(exp)
        {
            std::cout<<"EvalExpression constructor 1"<<std::endl;
        }
        
        EvalExpression(TrialExpressionType&& exp) : trialExp(std::move(exp))
        {
            std::cout<<"EvalExpression constructor 2"<<std::endl;
        }

        const TrialExpressionType& operator()() const
        {
            return trialExp();
        }
        
    };
    
}	// close namespace
#endif

//        /**********************************************************************/
//        const Derived& derived() const
//        {/*! A const reference to the Derived object
//          */
//            return *static_cast<const Derived*>(this);
//        }
