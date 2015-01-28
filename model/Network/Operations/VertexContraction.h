/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_VERTEXCONTRACTION_H_
#define model_VERTEXCONTRACTION_H_

#include <float.h> // for FLT_EPSILON, remove this once FlowComparison Class is created

#include <model/Network/Operations/VertexFinder.h>
#include <model/Network/Operations/EdgeFinder.h>
#include <model/Network/Operations/VertexInsertion.h>
#include <model/Network/Operations/VertexConnection.h>

namespace model {
	
	template <typename VertexType, typename EdgeType>
	class VertexContraction{
		
		typedef typename EdgeType::FlowType FlowType;
		
		typedef typename VertexType::NeighborContainerType NeighborContainerType; 
		
		typedef typename VertexFinder<VertexType>::isConstNetworkVertexType isConstNetworkVertexType;
		
		typedef typename EdgeFinder<EdgeType>::isConstNetworkEdgeType isConstNetworkEdgeType;
		
		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
		typedef boost::ptr_map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
		
		
		
		/* contractHelper ***************************************************/
		template <bool removeIsolatedNodes>
		void contractHelper(const size_t& i, const size_t& j){
			
			const isConstNetworkVertexType Vi(VertexFinder<VertexType>(networkVertexMapRef).node(i));
			assert(Vi.first && "CONTRACTING NON EXISTING VERTEX i");
			const isConstNetworkVertexType Vj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
			assert(Vj.first && "CONTRACTING NON EXISTING VERTEX j");
            assert(i!=j && "IN CONTRACTING (i,j), i AND j MUST BE DISTINCT");
			
			for (typename NeighborContainerType::const_iterator nIter=Vj.second->outNeighborhood().begin(); nIter!=Vj.second->outNeighborhood().end();++nIter){
				const size_t k(std::get<0>(nIter->second)->sID);
				if(k!=i){
					const FlowType fjk(std::get<1>(nIter->second)->flow);
					const isConstNetworkEdgeType eik(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,k));
					if (eik.first){
						const FlowType fik(eik.second->flow);
						VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(i,k);
						if((fik+fjk).norm()>FLT_EPSILON){
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fik+fjk);
						}
					}
					else{
						const isConstNetworkEdgeType eki(EdgeFinder<EdgeType>(networkEdgeMapRef).link(k,i));
						if(eki.first){
							const FlowType fki(eki.second->flow);
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(k,i);
							if((fki-fjk).norm()>FLT_EPSILON){
								VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fki-fjk);
							}
						}
						else {
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fjk);
						}					
					}
				}
			}
			
			for (typename NeighborContainerType::const_iterator nIter=Vj.second->inNeighborhood().begin(); nIter!=Vj.second->inNeighborhood().end();++nIter){
				const size_t k(std::get<0>(nIter->second)->sID);
				if (k!=i){
					const FlowType fkj(std::get<1>(nIter->second)->flow);
					const isConstNetworkEdgeType eik(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,k));
					if (eik.first){
						const FlowType fik(eik.second->flow);
						VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(i,k);
						if((fik-fkj).norm()>FLT_EPSILON){
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fik-fkj);
						}
					}
					else{
						const isConstNetworkEdgeType eki(EdgeFinder<EdgeType>(networkEdgeMapRef).link(k,i));
						if(eki.first){
							const FlowType fki(eki.second->flow);
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(k,i);
							if((fki+fkj).norm()>FLT_EPSILON){
								VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fki+fkj);
							}
						}
						else {
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fkj);
						}					
					}
				}
			}
			
			VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template remove<removeIsolatedNodes>(j); // removeIsolatedNodes is used here because annihilation leaves isolated nodes attached to j
			const isConstNetworkVertexType newVj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
			assert(!newVj.first && "CONTRACTING LEFT j");		
		}
		
		
		
		
		
	public:
		//		VertexContraction(const VertexType& V1in, const VertexType& V2in) : V1(V1in), V2(V2in) {}
		/* Constructor **********************************************/
		VertexContraction(NetworkVertexMapType& networkVertexMapRef_in,
		/*             */   NetworkEdgeMapType& networkEdgeMapRef_in) : networkVertexMapRef(networkVertexMapRef_in),
		/*                                                           */   networkEdgeMapRef(networkEdgeMapRef_in  ){}
		
		
		
		/* contract *************************************************/
		template <typename ...NodeArgTypes>
		void contract(const size_t& i, const size_t& j, const NodeArgTypes&... NodeInput){
			
			const size_t newID(VertexInsertion<VertexType>(networkVertexMapRef).insert(NodeInput...)); // CHANGE THIS LIKE EXPAND
			
			/* - Call contractHelper with removeIsolatedNodes=0 to make sure that j survives
			 * - This will not create isolated nodes
			 * - if i was connected to j now newID is connected to j
			 */
			contractHelper<0>(newID,i); 
			
			/* Now call contractHelper with removeIsolatedNodes=1 to make sure that isolated nodes are removed
			 * - if newID->j (or reverse) is the only connection, then newID will also be removed
			 * If newID->j (and reverse) do not exists, then newID will survive because remove<1> j does not affect newID
			 */
			contractHelper<1>(newID,j);
			
			// if j was isolated then newID will now also be isolated. Remove it.
			const isConstNetworkVertexType Vn(VertexFinder<VertexType>(networkVertexMapRef).node(newID));
			
			if(Vn.first)
            {
				if(Vn.second->is_isolated())
                { // WHAT IF THE ONLY NEIGHBOR OF I IS J? THEN NEW_ID IS ISOLATED
					VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template remove<0>(newID); 
				}
			}
			
			
		}
		
		/* contractSecond *******************************************/
		//		template <typename ...NodeArgTypes>
		void contractSecond(const size_t& i, const size_t& j){
			
			// If i->j (or reverse) exists and is the only connection this will remove both i and j
			// If i->j (and reverse) do not exists, then i will survive
			contractHelper<1>(i,j);
			// if j was isolated then newID will now also be isolated. Remove it.
			const isConstNetworkVertexType Vi(VertexFinder<VertexType>(networkVertexMapRef).node(i));
			
			if(Vi.first){
				if(Vi.second->is_isolated()){ // WHAT IF THE ONLY NEIGHBOR OF I IS J? THEN NEW_ID IS ISOLATED
					VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template remove<0>(i); 
				}
			}
			
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif
