/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakExpression_H_
#define model_BilinearWeakExpression_H_


#include <BilinearWeakExpression.h>


namespace model
{
    
    template <typename LHS, typename RHS>
    class WeakProblem;

    
    /**************************************************************************/
	/**************************************************************************/
    template <typename Derived>
	struct BilinearWeakExpression
    {
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
                
    };
    
    
    /**********************************************************************/
    template< typename B,typename L>
    WeakProblem<B,L> operator==(const BilinearWeakExpression<B>& bwe,
                               const LinearWeakExpression<L>& lwe)
    {
        return WeakProblem<B,L>(bwe.derived(),lwe.derived());
    }
    
    /**********************************************************************/
    template< typename B,typename L>
    WeakProblem<B,L> operator==(BilinearWeakExpression<B>&& bwe,
                                const LinearWeakExpression<L>& lwe)
    {
        return WeakProblem<B,L>(std::move(bwe.derived()),lwe.derived());
    }
    
    /**********************************************************************/
    template< typename B,typename L>
    WeakProblem<B,L> operator==(const BilinearWeakExpression<B>& bwe,
                                LinearWeakExpression<L>&& lwe)
    {
        return WeakProblem<B,L>(bwe.derived(),std::move(lwe.derived()));
    }
    /**********************************************************************/
    template< typename B,typename L>
    WeakProblem<B,L> operator==(BilinearWeakExpression<B>&& bwe,
                                LinearWeakExpression<L>&& lwe)
    {
        return WeakProblem<B,L>(std::move(bwe.derived()),std::move(lwe.derived()));
    }
    
}	// close namespace
#endif

