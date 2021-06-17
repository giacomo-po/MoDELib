/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ParallelFor_H_
#define model_ParallelFor_H_

#include <thread>
#include <vector>
#include <assert.h>
//#include <iostream>

#include <Model/Threads/EqualIteratorRange.h>

namespace model {
	
    struct ParallelFor :
    /* inheritance  */ private std::vector<std::thread>
    {
        
        const long int nThreads;
        
        
        ParallelFor(const long int& n = std::thread::hardware_concurrency()) :
        /* init list */ nThreads(n)
        {
//            std::cout<<"Creating ParallelFor with "<<nThreads<<" threads."<<std::endl;

            assert(nThreads>0 && "NUMBER OF THREADS MUST BE >0");
            this->reserve(nThreads);
        }
        
        template< typename AnyClass, class Function, class... Args >
        void runRange(const size_t& _begin, const size_t& _end, AnyClass* const C, Function&& foo, Args&&... args )
        {
            
            assert( _end > _begin && "RANGE MUST HAVE AT LEAST ONE ELEMENT");
            //std::vector<std::pair<size_t,size_t> > chunks;
            //chunks.reserve(max_threads);
            size_t range_size= (_end-_begin) / nThreads;
//            ParticleContainerIteratorType cur_iter = this->begin();
//            for(int i = 0; i < max_threads - 1; ++i)
//            {
//                ParticleContainerIteratorType last_iter = cur_iter;
//                std::advance(cur_iter, chunk_size);
//                chunks.push_back(std::make_pair(last_iter, cur_iter));
//            }
//            chunks.push_back(std::make_pair(cur_iter, this->end()));
            
            size_t currentBegin(_begin);
            for (long int i=0;i<nThreads;++i)
            {
                size_t currentEnd( (i<nThreads)? currentBegin+range_size : _end );
                this->emplace_back(foo,C,currentBegin,currentEnd,args...);
                currentBegin=currentEnd;
            }
            
            for (long int i=0;i<this->size();++i)
            {
                this->operator[](i).join();
            }
            
            this->clear();
            
        }
        
    };
	
	/************************************************************/	
	/************************************************************/
} // namespace model
#endif
