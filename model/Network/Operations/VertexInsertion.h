/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXINSERTION_H_
#define model_VERTEXINSERTION_H_

#include <memory> // for auto_ptr
#include <assert.h>
#include <boost/ptr_container/ptr_map.hpp>

namespace model {
	
	template <typename VertexType>
	class VertexInsertion{
		
		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
		//! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
	public:
//		typedef size_t ReturnType;

		/* Constructor ******************************/
		VertexInsertion(NetworkVertexMapType& networkVertexMapRef_in) : networkVertexMapRef(networkVertexMapRef_in) {}
		
		/* insert ***********************************/
		template <typename ...NodeArgTypes>
		size_t insert(const NodeArgTypes&... NodeInput){
			//! 1- Creates a new node
			std::auto_ptr<VertexType> pN (new VertexType(NodeInput...) );
			//! 2- Inserts the new node in NodeContainer
			size_t nodeID(pN->sID);
			assert(networkVertexMapRef.insert(nodeID , pN ).second && "CANNOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
			//! 3- Returns the static ID of the new node.
			return nodeID;
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif
