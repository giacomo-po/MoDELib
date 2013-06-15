/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DEPTHFIRST_H_
#define model_DEPTHFIRST_H_

#include <assert.h>
#include <tuple>

namespace model {
	
	template <typename NodeType, typename LinkType>
	class DepthFirst{		
		
		NodeType& nodeRef;
		
		
		/*****************************************************************************************/
		/* search ********************************************************************************/
		bool search (std::set<size_t>& searchedNodes, const size_t& ID, const size_t& N = ULONG_MAX) const{
			bool found(nodeRef.sID==ID);
			if (N!=0 && !found){
				assert(searchedNodes.insert(nodeRef.sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::const_iterator NeighborIter=nodeRef.neighborhood().begin();NeighborIter!=nodeRef.neighborhood().end();++NeighborIter){
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						found=std::get<0>(NeighborIter->second)->search(searchedNodes,ID,N-1); // ask if neighbor can reach
						if (found){
							break;
						}
					}
				}	
			}
			return found;
		}
		
		
		/*****************************************************************************************/
		/* execute *******************************************************************************/
		template <typename T>
		void execute(std::set<size_t>& searchedNodes, std::set<std::pair<size_t,size_t> >& searchedLinks, void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			(nodeRef.p_derived()->*Nfptr)(input); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(nodeRef.sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=nodeRef.neighborhood().begin();NeighborIter!=nodeRef.neighborhood().end();++NeighborIter){
					if (!std::get<2>(NeighborIter->second)==0){
						if (searchedLinks.find(std::get<1>(NeighborIter->second)->nodeIDPair)==searchedLinks.end()){  // neighbor not searched
							(std::get<1>(NeighborIter->second)->*Lfptr)(input); // execute Lfptr on connecting link
							assert(searchedLinks.insert(std::get<1>(NeighborIter->second)->nodeIDPair).second && "CANNOT INSERT CURRENT LINK IN SEARCHED LINKS"); // this node has been searched
						}
					}
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->execute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		
		/* nodeExecute *****************************************************************/
		template <typename T>
		void nodeExecute(std::set<size_t>& searchedNodes, void (Derived::*Nfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			(nodeRef.p_derived()->*Nfptr)(input); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(nodeRef.sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=nodeRef.neighborhood().begin();NeighborIter!=nodeRef.neighborhood().end();++NeighborIter){
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->nodeExecute(searchedNodes,Nfptr,input, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		
		/* nodeExecute *****************************************************************/
		void nodeExecute(std::set<size_t>& searchedNodes, void (Derived::*Nfptr)(void), const size_t& N = ULONG_MAX){
			(nodeRef.p_derived()->*Nfptr)(); // execute Nfptr on this node
			if (N!=0){
				assert(searchedNodes.insert(nodeRef.sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
				for (typename NeighborContainerType::iterator NeighborIter=nodeRef.neighborhood().begin();NeighborIter!=nodeRef.neighborhood().end();++NeighborIter){
					if (searchedNodes.find(std::get<0>(NeighborIter->second)->sID)==searchedNodes.end()){  // neighbor not searched
						std::get<0>(NeighborIter->second)->nodeExecute(searchedNodes,Nfptr, N-1); // continue executing on neighbor
					}
				}	
			}
		}
		
	public:
		
		/* Costructor with node ******************************************************************/
		DepthFirst(const NodeType& nodeRef_in) : nodeRef(nodeRef_in){ }

		/* search ********************************************************************************/
		bool search (const size_t& ID, const size_t& N = ULONG_MAX) const {
			std::set<size_t> searchedNodes;
			return search(searchedNodes,ID,N);
		}
		
		/* execute *******************************************************************************/
		template <typename T>
		void execute(void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			std::set<std::pair<size_t,size_t> > searchedLinks;
			execute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N);
		}
		
		/* nodeExecute ***************************************************************************/
		template <typename T>
		void nodeExecute(void (Derived::*Nfptr)(const T&), const T & input, const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			nodeExecute(searchedNodes,Nfptr,input, N);
		}
		
		void nodeExecute(void (Derived::*Nfptr)(void), const size_t& N = ULONG_MAX){
			std::set<size_t> searchedNodes;
			nodeExecute(searchedNodes,Nfptr, N);
		}
		
	};
		
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
} // namespace model
#endif



