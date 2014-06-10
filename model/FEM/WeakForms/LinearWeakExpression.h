/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakExpression_H_
#define model_LinearWeakExpression_H_




namespace model
{
    


    
    /**************************************************************************/
	/**************************************************************************/
    template <typename Derived>
	struct LinearWeakExpression
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
    
    
}	// close namespace
#endif

