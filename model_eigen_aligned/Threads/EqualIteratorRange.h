/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EqualIteratorRange_H_
#define model_EqualIteratorRange_H_

#include <deque>

#include <utility>  // std::pair
#include <iterator> // std::distance, std::advance
#include <assert.h> // std::distance, std::advance

namespace model
{
    
    template <typename IterType>
    struct EqualIteratorRange : public std::deque<std::pair<IterType,IterType>>
    {
        
        /**********************************************************************/
        EqualIteratorRange(const IterType& first,
                           const IterType& last,
                           const size_t& nRanges,
                           const bool& checkRanges=false)
        {
            
            assert(nRanges>0 && "NUMBER OF RANGES MUST BE >0");
            
            const auto nElments(std::distance(first,last));
            //assert( nElments > 0 && "RANGE MUST HAVE AT LEAST ONE ELEMENT");
            
            const size_t quotient = nElments / nRanges;
            const size_t remainder= nElments % nRanges;
            
            // Fill ranges
            IterType currentBegin(first);
            IterType currentEnd(first);
            
            for (size_t i=0;i<nRanges;++i)
            {
                std::advance(currentEnd,(i<remainder)? quotient+1 : quotient);
                this->emplace_back(currentBegin,currentEnd);
                currentBegin=currentEnd;
            }
            
            //check
            if(checkRanges)
            {
                auto temp(nElments*0);
                for (size_t i=0;i<nRanges;++i)
                {
                    temp+=std::distance(this->operator[](i).first,this->operator[](i).second);
                }
                assert(temp==nElments);
            }
            
        }
        
    };
    
} // namespace model
#endif
