/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXCONNECTION_H_
#define model_VERTEXCONNECTION_H_

//#include <iostream>
#include <memory>   // for auto_ptr
#include <utility>  // for std::pair
#include <assert.h>
#include <boost/ptr_container/ptr_map.hpp>

namespace model {
	
	template <typename VertexType, typename EdgeType>
	class VertexConnection{
		
		
		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
		typedef boost::ptr_map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
		
	public:
		
		/* Constructor **********************************************/
		VertexConnection(NetworkVertexMapType& networkVertexMapRef_in,
		/*            */ NetworkEdgeMapType&     networkEdgeMapRef_in) : networkVertexMapRef(networkVertexMapRef_in),
		/*                                                            */ networkEdgeMapRef(networkEdgeMapRef_in){}
		
		/* connect **************************************************/
//		template <typename FlowType>
		template <typename ...EdgeArgTypes>
		bool connect(const size_t& i, const size_t& j, const EdgeArgTypes&... edgeArgs){
			assert(i!=j && "TRYING TO CONNECT A VERTEX TO ITSELF.");
			typename NetworkVertexMapType::iterator Ni(networkVertexMapRef.find(i));
			assert(Ni!=networkVertexMapRef.end() && "TRYING TO CONNECT (i->j): VERTEX i DOES NOT EXIST.");
			typename NetworkVertexMapType::iterator Nj(networkVertexMapRef.find(j));
			assert(Nj!=networkVertexMapRef.end() && "TRYING TO CONNECT (i->j): VERTEX j DOES NOT EXIST.");			
			assert(networkEdgeMapRef.find(std::make_pair(i,j))==networkEdgeMapRef.end() && "EDGE ALREADY EXISTS. CANNOT CONNECT");
			assert(networkEdgeMapRef.find(std::make_pair(j,i))==networkEdgeMapRef.end() && "OPPOSITE EDGE ALREADY EXISTS. CANNOT CONNECT");
			std::auto_ptr<EdgeType> pL (new EdgeType(std::make_pair(Ni->second,Nj->second), edgeArgs...) );
			assert(networkEdgeMapRef.insert(std::make_pair(i,j), pL ).second && "CANNOT INSERT EDGE IN networkEdgeMapRef.");
			return true;
		}
		
		/* disconnect ***********************************************/
		template<bool removeIsolatedNodes>
		bool disconnect(const size_t& i, const size_t& j){
			typename NetworkEdgeMapType::iterator iterIJ(networkEdgeMapRef.find(std::make_pair(i,j))); // look for edge (i->j)
			if (iterIJ!=networkEdgeMapRef.end()){   // edge (i->j) is in the edgeMap
				networkEdgeMapRef.erase(iterIJ);	// remove (i->j)
			}
			else{ // edge (i->j) is not in the edgeMap
				typename NetworkEdgeMapType::iterator iterJI(networkEdgeMapRef.find(std::make_pair(j,i))); // look for edge (j->i)
				assert(iterJI!=networkEdgeMapRef.end() && "NEITHER (i->j) nor (j->i) ARE IN networkEdgeMapRef.");
				networkEdgeMapRef.erase(iterJI);	// remove (j->i)
			}
			if (removeIsolatedNodes){				
				typename NetworkVertexMapType::iterator Ni(networkVertexMapRef.find(i));
				assert(Ni!=networkVertexMapRef.end() && "NODE i DOES NOT EXIST");
				if (Ni->second->is_isolated()){
					networkVertexMapRef.erase(Ni); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(i), gives a bug
					//	WHY IS THIS TRIGGERED????				assert(networkVertexMapRef.find(i)==networkVertexMapRef.end() && "NODE i IS STILL IN networkVertexMapRef AFTER ERASE.");
				}
				typename NetworkVertexMapType::iterator Nj(networkVertexMapRef.find(j));
				assert(Nj!=networkVertexMapRef.end() && "NODE j DOES NOT EXIST");
				if (Nj->second->is_isolated()){
					networkVertexMapRef.erase(Nj); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(j), gives a bug
					// 	WHY IS THIS TRIGGERED????					assert(networkVertexMapRef.find(j)==networkVertexMapRef.end() && "NODE j IS STILL IN networkVertexMapRef AFTER ERASE.");
				}
			}
			return true;
		}
		
		/************************************************************/
		// remove (a node)
		template<bool removeIsolatedNodes>
		bool remove(const size_t& k){
//						std::cout<<"I'm here 2"<<std::endl;
			
			typename NetworkVertexMapType::iterator Nk(networkVertexMapRef.find(k));
			assert(Nk!=networkVertexMapRef.end() && "REMOVING NON-EXISTING VERTEX.");
			
			
			
			//bool success=0;
			//			isNetworkNodeType Nk=node(k);
			
			//			if (Nk.first){
			
			// Define BoolLinkMFPsize_t as a member-function-pointer of class EdgeType with size_t input and bool output
			typedef bool (EdgeType::*BoolLinkMFPsize_t)(const size_t &) const;
			BoolLinkMFPsize_t MFP = &EdgeType::isIncident;
			disconnect_if<removeIsolatedNodes>(MFP,k);
			//	disconnect_if<removeIsolatedNodes>(&EdgeType::isIncident,k);	// WHY IS THIS NOT WORKING????????
			
			typename NetworkVertexMapType::iterator Nkk(networkVertexMapRef.find(k));	
			if (Nkk!=networkVertexMapRef.end()){
				networkVertexMapRef.erase(Nkk); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(j), gives a bug
			}
			assert(networkVertexMapRef.find(k)==networkVertexMapRef.end() && "NODE k IS STILL IN networkVertexMapRef AFTER ERASE.");
			
			
			//				success=NetworkNodeContainerType::erase(k) == 1;
			//			}			
			
			return true;
		}
		
		
		/************************************************************/
		// disconnect_if
		template<bool removeIsolatedNodes>
		size_t disconnect_if(bool (EdgeType::*Lfptr)(void) const){
			// LOOPING AND ERASING IS DANGEROUS!!! YOU MUST INCREMENT THE ITERATOR BEFORE ERASING
			
//			std::cout<<"I'm here 0"<<std::endl;
			size_t count = 0;
			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();){
				if ( (edgeIter->second->*Lfptr)() ) {
					typename NetworkEdgeMapType::iterator toBeDisconnected(edgeIter);
					++edgeIter;
					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second->source->sID,toBeDisconnected->second->sink->sID);					
				}
				else{
					++edgeIter;
				}
			}
			return count;
		}
		
		/************************************************************/
		// disconnect_if
		template<bool removeIsolatedNodes, typename T>
		size_t disconnect_if(bool (EdgeType::*Lfptr)(const T &) const, const T & input){
			size_t count = 0;
			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();){
				if ( (edgeIter->second->*Lfptr)(input) ) {
					// LOOPING AND ERASING IS DANGEROUS!!! YOU MUST INCREMENT THE ITERATOR BEFORE ERASING
					typename NetworkEdgeMapType::iterator toBeDisconnected(edgeIter);
					++edgeIter;
					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second->source->sID,toBeDisconnected->second->sink->sID);					
				}
				else{
					++edgeIter;
				}
			}
			return count;
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif


//#include <model/Network/Operations/VertexFinder.h>
//#include <model/Network/Operations/EdgeFinder.h>
//		typedef typename VertexFinder<VertexType>::isNetworkVertexType isNetworkVertexType;
//		typedef typename EdgeFinder<EdgeType>::isNetworkEdgeType isNetworkEdgeType;



//		
//		
//		
//		/************************************************************/
//		// remove (a node)
//		template<bool removeIsolatedNodes>
//		bool remove(const size_t& k){
//			
//			//			std::cout<<"About to remove"<<std::endl;
//			
//			isNetworkVertexType Nk(VertexFinder<VertexType>(networkVertexMapRef).node(k));
//			assert(Nk.first && "TRYING TO REMOVE NON-EXISTING VERTEX");
//			
//			//			bool success=0;
//			
//			//			if (Nk.first){
//			//						}			
//			// Define BoolLinkMFPsize_t as a member-function-pointer of class EdgeType with size_t input and bool output
//			typedef bool (EdgeType::*BoolEdgeMFPsize_t)(const size_t &) const;
//			BoolEdgeMFPsize_t MFP = &EdgeType::isIncident;
//			disconnect_if<removeIsolatedNodes>(MFP,k);
//			//	disconnect_if<removeIsolatedNodes>(&EdgeType::isIncident,k);	// WHY IS THIS NOT WORKING????????
//			//success=NetworkNodeContainerType::erase(k) == 1;
//			
//			bool success(networkVertexMapRef.erase(k) == 1); // THIS IS UNSAFE, CHANGE IN THE FOLLOWING
//			//		assert(networkVertexMapRef.erase(k) == 1 && "UNABLE TO REMOVE VERTEX FROM networkVertexMapRef");
//			
//			return success;
//		}
//		
//		/************************************************************/
//		// disconnect_if
//		template<bool removeIsolatedNodes>
//		size_t disconnect_if(bool (EdgeType::*Lfptr)(void) const){
////			std::cout<<"About to disconnect_if"<<std::endl;
//
//			size_t count = 0;
//			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();){
//				if ( (edgeIter->second->*Lfptr)() ) {
//					// LOOPING AND ERASING IS DANGEROUS!!! YOU MUST INCREMENT THE ITERATOR BEFORE ERASING
//					typename NetworkEdgeMapType::iterator toBeDisconnected = edgeIter;
//					++edgeIter;
//					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second->source->sID,toBeDisconnected->second->sink->sID);					
//				}
//				else{
//					++edgeIter;
//				}
//			}
//			return count;
//		}
//		
//		/************************************************************/
//		// disconnect_if
//		template<bool removeIsolatedNodes, typename T>
//		size_t disconnect_if(bool (EdgeType::*Lfptr)(const T &) const, const T & input){
//			
////			std::cout<<"About to disconnect_if_T"<<std::endl;
//
//			
//			size_t count = 0;
//			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();){
//				if ( (edgeIter->second->*Lfptr)(input) ) {
////					std::cout<<"I'm here 1"<<std::endl;
//					// LOOPING AND ERASING IS DANGEROUS!!! YOU MUST INCREMENT THE ITERATOR BEFORE ERASING
//					typename NetworkEdgeMapType::iterator toBeDisconnected = edgeIter;
////										std::cout<<"I'm here 2"<<std::endl;
//					++edgeIter;
////										std::cout<<"I'm here 3"<<std::endl;
//					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second->source->sID,toBeDisconnected->second->sink->sID);					
////									std::cout<<"I'm here 4"<<std::endl;
//				}
//				else{
//					++edgeIter;
//				}
//			}
//			
////			std::cout<<"I'm here 5"<<std::endl;
//
//			return count;
//		}



//isNetworkVertexType Nj=VertexFinder<VertexType>(networkVertexMapRef).node(j);				
//				//if (Nj.first ) {
//					if (Nj.second->is_isolated()){
//						//remove<0>(Nj.second->sID);
//						//networkVertexMapRef.erase(j);
//						networkVertexMapRef.erase(networkVertexMapRef.find(j));
//						assert(networkVertexMapRef.find(j)==networkVertexMapRef.end() && "NODE j IS IN networkVertexMapRef AFTER ERASE.");
//					}
//				//}




//			bool success=( networkEdgeMapRef.erase( std::make_pair(i,j) ) + networkEdgeMapRef.erase( std::make_pair(j,i) ) ) == 1;


//			NetworkVertexMapType::iterator iterI(networkVertexMapRef.find(i));
//			assert()

// remove isoltaed nodes

//isNetworkVertexType Ni=VertexFinder<VertexType>(networkVertexMapRef).node(i);

//if (Ni.first) {

//}


//		/************************************************************/
//		// disconnect
//		template<bool removeIsolatedNodes>
//		bool disconnect(const size_t& i, const size_t& j){
//			
//			isNetworkVertexType Ni(VertexFinder<VertexType>(networkVertexMapRef).node(i));
//			assert(Ni.first && "COULD NOT FIND NODE i IN DISCONNECTING");
//			
//			isNetworkVertexType Nj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
//			assert(Nj.first && "COULD NOT FIND NODE i IN DISCONNECTING");
//			
//			isNetworkEdgeType Eij(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,j));
//			if (Eij.first){
//				assert(networkEdgeMapRef.erase(std::make_pair(i,j))==1 && "UNABLE TO ERASE EDGE I->J FROM networkEdgeMapRef.");
//			}
//			else{
//				isNetworkEdgeType Eji(EdgeFinder<EdgeType>(networkEdgeMapRef).link(j,i));
//				assert(Eji.first && "IN DISCONNECTING (I,J), NEITHER I->J NOR J->I EXIST.");
//				assert(networkEdgeMapRef.erase(std::make_pair(j,i))==1 && "UNABLE TO ERASE EDGE J->I FROM networkEdgeMapRef.");
//			}
//			
//			
////			unsigned int dIJ();
//			
////						std::cout<<"About to disconnect"<<std::endl;
//			
////			bool success=( networkEdgeMapRef.erase( std::make_pair(i,j) ) + networkEdgeMapRef.erase( std::make_pair(j,i) ) ) == 1;
////			assert(( networkEdgeMapRef.erase(std::make_pair(i,j)) + networkEdgeMapRef.erase(std::make_pair(j,i)) ) == 1 && "UNABLE TO REMOVE EDGE FROM networkEdgeMapRef");
//
//			
//			// remove isoltaed nodes
////			if (success){
////				isNetworkVertexType Ni(VertexFinder<VertexType>(networkVertexMapRef).node(i));
////				isNetworkVertexType Nj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
//				
//				if (Ni.first && removeIsolatedNodes) {
//					if (Ni.second->is_isolated()){
//						unsigned int remv=networkVertexMapRef.erase(i);
//						std::cout<<"Number of Erases of i is:"<<remv<<std::endl;
//						assert(remv==1 && "UNABLE TO REMOVE ISOLATED VERTEX from networkVertexMapRef");
//
////						assert(networkVertexMapRef.erase(i) == 1 && "UNABLE TO REMOVE ISOLATED VERTEX from networkVertexMapRef");
//						//remove<0>(Ni.second->sID);
//					}
//				}
//				
//				if (Nj.first && removeIsolatedNodes) {
//					if (Nj.second->is_isolated()){
//						assert(networkVertexMapRef.erase(j) == 1 && "UNABLE TO REMOVE ISOLATED VERTEX from networkVertexMapRef");
////						remove<0>(Nj.second->sID);
//					}
//				}
////			}
//			
////			return success;
////			std::cout<<"Finish Disconnect"<<std::endl;
//			return true;
//		}
