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
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#endif

#include <iostream>
#include <memory>   // for auto_ptr
#include <utility>  // for std::pair
//#include <assert.h>
#include <boost/ptr_container/ptr_map.hpp>

namespace model {
	
	template <typename VertexType, typename EdgeType>
	class ParallelExecute{
		
		
		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
		typedef boost::ptr_map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
		
	public:
		
		/* Constructor **********************************************/
		ParallelExecute(NetworkVertexMapType& networkVertexMapRef_in,
		/*            */ NetworkEdgeMapType&     networkEdgeMapRef_in) : networkVertexMapRef(networkVertexMapRef_in),
		/*                                                            */ networkEdgeMapRef(networkEdgeMapRef_in){}
		

		/* execute **********************************************/
		void execute(void (EdgeType::*Lfptr)(void)){
#ifdef _OPENMP
			#pragma omp parallel 
			#pragma omp single
{
std::cout<<"[executing with "<<omp_get_num_threads()<<" threads]"<<std::flush;
#endif
			
			for (typename NetworkEdgeMapType::iterator linkIter=networkEdgeMapRef.begin();linkIter!=networkEdgeMapRef.end();++linkIter){
#ifdef _OPENMP
 		    #pragma omp task firstprivate(linkIter)	
#endif			
			(linkIter->second->*Lfptr)();
			}
#ifdef _OPENMP
		  	#pragma omp taskwait
}
#endif			
		}

		/* execute **********************************************/
		template <typename T>
		void execute(void (EdgeType::*Lfptr)(const T&), const T & input){
#ifdef _OPENMP
			#pragma omp parallel 
			#pragma omp single
{
std::cout<<"[executing with "<<omp_get_num_threads()<<" threads]"<<std::flush;
#endif
			
			for (typename NetworkEdgeMapType::iterator linkIter=networkEdgeMapRef.begin();linkIter!=networkEdgeMapRef.end();++linkIter){
#ifdef _OPENMP
 		    #pragma omp task firstprivate(linkIter)	
#endif			
			(linkIter->second->*Lfptr)(input);
			}
#ifdef _OPENMP
		  	#pragma omp taskwait
}
#endif			
		}


		/* vertexExecute **********************************************/
		void execute(void (VertexType::*Vfptr)(void)){
#ifdef _OPENMP
			#pragma omp parallel 
			#pragma omp single
{
std::cout<<"[executing with "<<omp_get_num_threads()<<" threads]"<<std::flush;
#endif
			
			for (typename NetworkVertexMapType::iterator vertexIter=networkVertexMapRef.begin();vertexIter!=networkVertexMapRef.end();++vertexIter){
#ifdef _OPENMP
 		    #pragma omp task firstprivate(vertexIter)	
#endif			
			(vertexIter->second->*Vfptr)();
			}
#ifdef _OPENMP
		  	#pragma omp taskwait
}
#endif			
		}

		/* execute **********************************************/
		template <typename T>
		void execute(void (VertexType::*Vfptr)(const T&), const T & input){
#ifdef _OPENMP
			#pragma omp parallel 
			#pragma omp single
{
std::cout<<"[executing with "<<omp_get_num_threads()<<" threads]"<<std::flush;
#endif
			
			for (typename NetworkVertexMapType::iterator vertexIter=networkVertexMapRef.begin();vertexIter!=networkVertexMapRef.end();++vertexIter){
#ifdef _OPENMP
 		    #pragma omp task firstprivate(vertexIter)	
#endif			
			(vertexIter->second->*Vfptr)(input);
			}
#ifdef _OPENMP
		  	#pragma omp taskwait
}
#endif			
		}


		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

