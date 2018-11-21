/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialExpressionBase_H_
#define model_TrialExpressionBase_H_

#include <iostream>
#include <iomanip>
#include <TypeTraits.h>
#include <TestExpression.h>
#include <EvalExpression.h>
#include <TrialDomainView.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    /*!\brief A class template that provides the base for all the expressions
     * involving a TrialFunction. The expression-template mechanism is based
     * on the CRTP pattern.
     */
    template<typename Derived>
    struct TrialExpressionBase
    {
        TrialExpressionBase()= default;
        TrialExpressionBase(const TrialExpressionBase<Derived>& other) = delete; // do not copy TrialExpressionBase
        TrialExpressionBase(TrialExpressionBase<Derived>&& other) = default;
        
        /**********************************************************************/
        Derived& derived()
        {/*!\returns A reference to the Derived object
          */
            return *static_cast<Derived*>(this);
        }
        
        /**********************************************************************/
        Derived* p_derived()
        {/*!\returns A pointer to the Derived object
          */
            return  static_cast<Derived*>(this);
        }
        
        /**********************************************************************/
        const Derived& derived() const
        {/*! A const reference to the Derived object
          */
            return *static_cast<const Derived*>(this);
        }
        
        /**********************************************************************/
        const Derived* p_derived() const
        {/*!\returns A pointer to the const Derived object
          */
            return  static_cast<const Derived*>(this);
        }
        
        /**********************************************************************/
        operator const Derived&() const
        {/*!\returns the const Derived object (cast operator)
          */
            return derived();
        }
        
        /**********************************************************************/
        TrialDomainView<Derived,0> onDomain() const
        {
            return TrialDomainView<Derived,0>(eval(*this));
        }
        
        /**********************************************************************/
        TrialDomainView<Derived,1> onBoundary() const
        {
            return TrialDomainView<Derived,1>(eval(*this));
        }
        
    };
    
    /**************************************************************************/
    template<typename Derived>
    TestExpression<Derived> test(const TrialExpressionBase<Derived>& trial)
    {
        return TestExpression<Derived>(trial.derived());
    }
    
    template<typename Derived>
    TestExpression<Derived> test(TrialExpressionBase<Derived>&& trial)
    {
        return TestExpression<Derived>(std::move(trial.derived()));
    }
    
    /**************************************************************************/
    template<typename Derived>
    EvalExpression<Derived> eval(const TrialExpressionBase<Derived>& trial)
    {
        return EvalExpression<Derived>(trial.derived());
    }
    
    template<typename Derived>
    EvalExpression<Derived> eval(TrialExpressionBase<Derived>&& trial)
    {
        return EvalExpression<Derived>(std::move(trial.derived()));
    }

}	// close namespace
#endif


//    template<typename Derived>
//    TestExpression<typename std::remove_reference<Derived>::type> test(Derived&& trial)
//    {
//    EXPERIMENTING WITH STD::FORWARD TO REMOVE ALL OPERATOR AND CONSTRUCTOR OVERLOADS
//        std::cout<<"test function 2"<<std::endl;
//
//        return TestExpression<typename std::remove_reference<Derived>::type>(std::forward<Derived>(trial));
//    }
//
//


//        /**********************************************************************/
//        EvalExpression<Derived> eval() const
//        {
//            return EvalExpression<Derived>(derived());
//        }
