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
#include <map> 

#include <model/Network/Operations/VertexFinder.h>
#include <model/Network/Operations/EdgeFinder.h>
#include <model/Network/Operations/VertexInsertion.h>
#include <model/Network/Operations/VertexConnection.h>
#include <model/MPI/MPIcout.h>

#define VerboseContract(N,x) if(verboseContract>=N){model::cout<<x;}


namespace model
{
    
    /****************************************************************/
    /****************************************************************/
    template <typename VertexType, typename EdgeType>
    struct ContractingVertices
    {
        
        typedef typename EdgeFinder<EdgeType>::isConstNetworkEdgeType isConstNetworkEdgeType;
        
        const VertexType& v0;
        const VertexType& v1;
        const isConstNetworkEdgeType E;
        
        ContractingVertices(const VertexType& v0in,const VertexType& v1in, const isConstNetworkEdgeType& E_in) :
        /* init */v0(v0in),
        /* init */v1(v1in),
        /* init */E(E_in)
        {}
        
    };
	
    /****************************************************************/
    /****************************************************************/
	template <typename VertexType, typename EdgeType>
	class VertexContraction
    {
		
		typedef typename EdgeType::FlowType FlowType;
		typedef typename VertexType::NeighborContainerType NeighborContainerType;
		typedef typename VertexFinder<VertexType>::isConstNetworkVertexType isConstNetworkVertexType;
		typedef typename EdgeFinder<EdgeType>::isConstNetworkEdgeType isConstNetworkEdgeType;
//        typedef typename EdgeFinder<EdgeType>::isNetworkEdgeType isNetworkEdgeType;
        typedef std::map<size_t,VertexType> NetworkVertexMapType;
        typedef std::map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;

        //! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;

		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
        /**********************************************************************/
		template <bool removeIsolatedNodes>
		void contractHelper(const size_t& i, const size_t& j)
        {/* In contractHelper, j is destroyed
          */
			
			const isConstNetworkVertexType Vi(VertexFinder<VertexType>(networkVertexMapRef).node(i));
			assert(Vi.first && "CONTRACTING NON EXISTING VERTEX i");
			const isConstNetworkVertexType Vj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
			assert(Vj.first && "CONTRACTING NON EXISTING VERTEX j");
            assert(i!=j && "IN CONTRACTING (i,j), i AND j MUST BE DISTINCT");
			
            // Loop over out-neighbors of j and connect to i
			for (typename NeighborContainerType::const_iterator nIter=Vj.second->outNeighborhood().begin(); nIter!=Vj.second->outNeighborhood().end();++nIter)
            {
				const size_t k(std::get<0>(nIter->second)->sID);
				if(k!=i)
                {
					const FlowType fjk(std::get<1>(nIter->second)->flow);
					const isConstNetworkEdgeType eik(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,k));
					if (eik.first)
                    {
						const FlowType fik(eik.second->flow);
						VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(i,k);
						if((fik+fjk).squaredNorm()>FLT_EPSILON)
                        {
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fik+fjk);
						}
					}
					else{
						const isConstNetworkEdgeType eki(EdgeFinder<EdgeType>(networkEdgeMapRef).link(k,i));
						if(eki.first)
                        {
							const FlowType fki(eki.second->flow);
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(k,i);
							if((fki-fjk).squaredNorm()>FLT_EPSILON)
                            {
								VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fki-fjk);
							}
						}
						else
                        {
//							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fjk);
                            VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,EdgeRef<EdgeType>(*std::get<1>(nIter->second)));
                        }
					}
				}
			}
			
            // Loop over in-neighbors of j and connect to i
			for (typename NeighborContainerType::const_iterator nIter=Vj.second->inNeighborhood().begin(); nIter!=Vj.second->inNeighborhood().end();++nIter)
            {
				const size_t k(std::get<0>(nIter->second)->sID);
				if (k!=i)
                {
					const FlowType fkj(std::get<1>(nIter->second)->flow);
					const isConstNetworkEdgeType eik(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,k));
					if (eik.first)
                    {
						const FlowType fik(eik.second->flow);
						VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(i,k);
						if((fik-fkj).squaredNorm()>FLT_EPSILON)
                        {
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,k,fik-fkj);
						}
					}
					else{
						const isConstNetworkEdgeType eki(EdgeFinder<EdgeType>(networkEdgeMapRef).link(k,i));
						if(eki.first)
                        {
							const FlowType fki(eki.second->flow);
							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(k,i);
							if((fki+fkj).squaredNorm()>FLT_EPSILON)
                            {
								VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fki+fkj);
							}
						}
						else
                        {
//							VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,fkj);
                            VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(k,i,EdgeRef<EdgeType>(*std::get<1>(nIter->second)));
                        }
					}
				}
			}
			
			VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template remove<removeIsolatedNodes>(j); // removeIsolatedNodes is used here because annihilation leaves isolated nodes attached to j
			const isConstNetworkVertexType newVj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
			assert(!newVj.first && "CONTRACTING LEFT j");		
		}
		
	public:

        static int verboseContract;
        
        /**********************************************************************/
		VertexContraction(NetworkVertexMapType& networkVertexMapRef_in,
		/*             */   NetworkEdgeMapType& networkEdgeMapRef_in) : networkVertexMapRef(networkVertexMapRef_in),
		/*                                                           */   networkEdgeMapRef(networkEdgeMapRef_in  )
        {}
		
        /**********************************************************************/
		template <typename ...NodeArgTypes>
		void contract(const size_t& i, const size_t& j, const NodeArgTypes&... NodeInput)
        {
            VerboseContract(1,"contracting "<<i<<" "<<j<<std::flush;);
            const isConstNetworkVertexType Vi(VertexFinder<VertexType>(networkVertexMapRef).node(i));
            assert(Vi.first && "CONTRACTING NON EXISTING VERTEX i");
            const isConstNetworkVertexType Vj(VertexFinder<VertexType>(networkVertexMapRef).node(j));
            assert(Vj.first && "CONTRACTING NON EXISTING VERTEX j");


            const isConstNetworkEdgeType tempE1(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,j));
            const isConstNetworkEdgeType tempE2(tempE1.first? tempE1 : EdgeFinder<EdgeType>(networkEdgeMapRef).link(j,i));
            
			const size_t newID(VertexInsertion<VertexType>(networkVertexMapRef).insert(ContractingVertices<VertexType,EdgeType>(*Vi.second,*Vj.second,tempE2),NodeInput...).first->first); // CHANGE THIS LIKE EXPAND
			VerboseContract(1," into "<<newID<<std::endl;);
            
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
		
        /**********************************************************************/
		void contractSecond(const size_t& i, const size_t& j)
        {
            VerboseContract(1,"contracting second of "<<i<<" "<<j<<std::endl;);

			// If i->j (or reverse) exists and is the only connection this will remove both i and j
			// If i->j (and reverse) do not exists, then i will survive
			contractHelper<1>(i,j);
            
			// if j was isolated then newID will now also be isolated. Remove it.
			const isConstNetworkVertexType Vi(VertexFinder<VertexType>(networkVertexMapRef).node(i));
			if(Vi.first)
            {
				if(Vi.second->is_isolated())
                { // WHAT IF THE ONLY NEIGHBOR OF I IS J? THEN NEW_ID IS ISOLATED
					VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template remove<0>(i); 
				}
			}
			
		}
		
	};
    
    template <typename VertexType, typename EdgeType>
    int VertexContraction<VertexType,EdgeType>::verboseContract=0;
	
} // namespace model
#endif
