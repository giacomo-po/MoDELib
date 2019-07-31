/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SymInteractionDeque_H_
#define model_SymInteractionDeque_H_

#include <deque>

#include <utility>  // std::pair

namespace model
{
	
    template <typename IterType>
    struct SymInteractionDeque :
    /* inheritance  */ public std::deque<std::pair<IterType,IterType> >
    {
        
        
        /**********************************************************************/
        SymInteractionDeque(const IterType& first, const IterType& last)
        {
            for (IterType iterA=first;iterA!=last;++iterA)
            {
                for (IterType iterB=iterA;iterB!=last;++iterB)
                {
                    this->emplace_back(iterA,iterB);
                }
            }
            
        }
        
    };
    
    
} // namespace model
#endif

