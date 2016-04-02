/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EvalExpression_H_
#define model_EvalExpression_H_

namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
	template<typename T>
    struct EvalExpression
    {
        
        const T wrappedExp;
        
        EvalExpression(const T& exp) : wrappedExp(exp)
        {
        
        }
        
//        /**********************************************************************/
//        const Derived& derived() const
//        {/*! A const reference to the Derived object
//          */
//            return *static_cast<const Derived*>(this);
//        }
        
    };
    
}	// close namespace
#endif