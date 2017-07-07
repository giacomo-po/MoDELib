/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PARALLELEXECUTE_H_
#define model_PARALLELEXECUTE_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iterator> // std::advance
#include <utility>  // std::pair
#include <map>
#include <model/Threads/EqualIteratorRange.h>


namespace model
{
	
	template <typename VertexType, typename EdgeType>
	class ParallelExecute
    {
		

        typedef std::map<size_t,VertexType> NetworkVertexMapType;
//		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
//		typedef boost::ptr_map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
        typedef std::map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
        //! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
		
	public:
		
		/* Constructor ********************************************************/
		ParallelExecute(NetworkVertexMapType& networkVertexMapRef_in,
        /*            */ NetworkEdgeMapType&     networkEdgeMapRef_in) : networkVertexMapRef(networkVertexMapRef_in),
		/*                                                            */ networkEdgeMapRef(networkEdgeMapRef_in){}
		
        
		/* execute ************************************************************/
		void execute(void (EdgeType::*Lfptr)(void))
        {
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
            EqualIteratorRange<typename NetworkEdgeMapType::iterator> eir(networkEdgeMapRef.begin(),networkEdgeMapRef.end(),nThreads);

#pragma omp parallel for
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (typename NetworkEdgeMapType::iterator linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                {
                    (linkIter->second.*Lfptr)();
                }
            }

//            for (unsigned int k=0;k<networkEdgeMapRef.size();++k)
//            {
//                typename NetworkEdgeMapType::iterator linkIter(networkEdgeMapRef.begin()); //  the data within a parallel region is private to each thread
//                std::advance(linkIter,k);
//                (linkIter->second.*Lfptr)();
//            }
#else
            for (typename NetworkEdgeMapType::iterator linkIter=networkEdgeMapRef.begin();linkIter!=networkEdgeMapRef.end();++linkIter)
            {
                (linkIter->second.*Lfptr)();
            }
#endif
		}
        
		/* execute ************************************************************/
		template <typename T>
		void execute(void (EdgeType::*Lfptr)(const T&), const T & input)
        {
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
            EqualIteratorRange<typename NetworkEdgeMapType::iterator> eir(networkEdgeMapRef.begin(),networkEdgeMapRef.end(),nThreads);
            
#pragma omp parallel for
            for (int thread=0;thread<eir.size();thread++)
            {
                for (typename NetworkEdgeMapType::iterator linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                {
                    (linkIter->second.*Lfptr)(input);
                }
            }
            //#pragma omp parallel for
//            for (unsigned int k=0;k<networkEdgeMapRef.size();++k){
//                typename NetworkEdgeMapType::iterator linkIter(networkEdgeMapRef.begin()); //  the data within a parallel region is private to each thread
//                std::advance(linkIter,k);
//                (linkIter->second.*Lfptr)(input);
//            }
#else
            for (typename NetworkEdgeMapType::iterator linkIter=networkEdgeMapRef.begin();linkIter!=networkEdgeMapRef.end();++linkIter)
            {
                (linkIter->second.*Lfptr)(input);
            }
#endif
		}

	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

