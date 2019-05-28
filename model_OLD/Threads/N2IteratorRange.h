

#ifndef model_N2IteratorRange_h
#define model_N2IteratorRange_h

#include <deque>
#include <utility>
#include <iterator>

namespace model
{
    
    template<typename IteratorType>
    struct N2IteratorRange : public std::deque<std::pair<IteratorType,IteratorType>>
    {
        
        N2IteratorRange(const IteratorType& first,const IteratorType& last,const size_t& nDiv)
        {
            size_t distance(std::distance(first,last));
            
            const size_t interactions(distance*(distance+1)/2);
            const size_t quotient(interactions/nDiv);
            
            IteratorType currentBegin(first);
            IteratorType currentEnd(first);
            
            size_t weight=0;
            while(currentEnd!=last)
            {
                weight+=distance;
                currentEnd++;

                if(weight>=quotient)
                {
                    this->emplace_back(currentBegin,currentEnd);
                    currentBegin=currentEnd;
                    weight=0;
                }
                distance--;
//                currentEnd++;
            }
            if(currentBegin!=last)
            {
                this->emplace_back(currentBegin,last);
            }
            
            assert(this->size()<=nDiv);

            
        }
        
    };
}


#endif
