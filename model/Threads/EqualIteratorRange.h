/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EqualIteratorRange_H_
#define model_EqualIteratorRange_H_

#include <vector>

#include <utility>  // std::pair
#include <iterator> // std::distance, std::advance
#include <assert.h> // std::distance, std::advance

namespace model
{
	
    template <typename ContainerType>
    struct EqualIteratorRange :
    /* inheritance  */ public std::vector<std::pair<typename ContainerType::const_iterator,
    /*                                           */ typename ContainerType::const_iterator>>
    {
        
        typedef typename ContainerType::const_iterator IterType;
        
        /**********************************************************************/
        EqualIteratorRange(const IterType& first, const IterType& last, const int& nRanges)
        {
            assert(nRanges>0 && "NUMBER OF RANGES MUST BE >0");
            
            const auto nElments(std::distance(first,last));
            assert( nElments > 0 && "RANGE MUST HAVE AT LEAST ONE ELEMENT");
            
            const size_t quotient = nElments / nRanges;
            const size_t remainder= nElments % nRanges;
            
            this->reserve(nRanges);
            
            IterType currentBegin(first);
            IterType currentEnd(first);
            
            for (int i=0;i<nRanges;++i)
            {
                std::advance(currentEnd,(i<remainder)? quotient+1 : quotient);
                this->emplace_back(currentBegin,currentEnd);
                currentBegin=currentEnd;
            }
            
        }
        
    };
    
    
	
	/************************************************************/
	/************************************************************/
} // namespace model
#endif
