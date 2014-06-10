/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DirichletCondition_H_
#define model_DirichletCondition_H_

#include <deque>


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename T>
	struct DirichletCondition
    {
        
//        const TrialFunctionType& tf;
        
        std::deque<size_t> gIDcontainer;
        
        /**********************************************************************/
        template <typename FiniteElementType>
        DirichletCondition(const FiniteElementType& fe)
//        DirichletCondition(const TrialFunctionType& tf_in) : tf(tf_in)
        {
//            typedef typename TrialFunctionType::NodeType NodeType;
            
            for (int n=0;n<fe.nodeSize();++n)
            {
//                if(T::template)
//                {
//                    
//                }
            }
            
        }
        
    };
    
    
    
    
    
}	// close namespace
#endif