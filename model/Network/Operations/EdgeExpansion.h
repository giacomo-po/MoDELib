/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EDGEEXPANSION_H_
#define model_EDGEEXPANSION_H_

#include <utility>  // for std::pair
#include <model/Network/Operations/EdgeFinder.h>
#include <model/Network/Operations/VertexInsertion.h>
#include <model/Network/Operations/VertexConnection.h>
#include <model/MPI/MPIcout.h>

#define VerboseExpand(N,x) if(verboseExpand>=N){model::cout<<x;}


namespace model
{
	
	/****************************************************************/
	/****************************************************************/
	template <typename EdgeType>
	struct EdgeRef
    {
	
		const EdgeType& E;
		EdgeRef(const EdgeType& Ein) : E(Ein) {}

	};

	
	/****************************************************************/
	/****************************************************************/
	template <typename VertexType, typename EdgeType>
	class EdgeExpansion
    {

		
		typedef typename EdgeFinder<EdgeType>::isConstNetworkEdgeType isConstNetworkEdgeType;

		
		//typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
        typedef typename VertexType::NetworkNodeContainerType NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
		//typedef boost::ptr_map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;
        typedef typename VertexType::NetworkLinkContainerType NetworkEdgeMapType;
		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
	public:
        
        static int verboseExpand;

        
		/* Constructor **********************************************/
		EdgeExpansion(NetworkVertexMapType& networkVertexMapRef_in,
		/*         */ NetworkEdgeMapType&     networkEdgeMapRef_in) :
        /* init list */ networkVertexMapRef(networkVertexMapRef_in),
		/* init list */ networkEdgeMapRef(networkEdgeMapRef_in)
        {/*! Initializes internal references 
          */
        }
		
		/* expand ***************************************************/
		template <typename ...NodeArgTypes>
		std::pair<typename NetworkVertexMapType::iterator,bool> expand(const size_t& i, const size_t& j, const NodeArgTypes&... Args)
        {
            VerboseExpand(1,"expanding "<<i<<"->"<<j<<std::endl;);
			isConstNetworkEdgeType Lij(EdgeFinder<EdgeType>(networkEdgeMapRef).link(i,j));
			assert(Lij.first && "EXPANDING NON-EXISTING LINK.");
            const std::pair<typename NetworkVertexMapType::iterator,bool> temp(VertexInsertion<VertexType>(networkVertexMapRef).insert(EdgeRef<EdgeType>(*Lij.second),Args...));
//			const size_t newID(temp.first->first);
			VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(i,temp.first->first,EdgeRef<EdgeType>(*Lij.second));
			VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).connect(temp.first->first,j,EdgeRef<EdgeType>(*Lij.second));
			VertexConnection<VertexType,EdgeType>(networkVertexMapRef,networkEdgeMapRef).template disconnect<0>(i,j);
			return temp;
		}
		
	};
    
    template <typename VertexType, typename EdgeType>
    int EdgeExpansion<VertexType,EdgeType>::verboseExpand=0;
	
	/****************************************************************/
	/****************************************************************/
} // namespace model
#endif
