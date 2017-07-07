/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXINSERTION_H_
#define model_VERTEXINSERTION_H_

//#include <memory> // std::auto_ptr
#include <assert.h>
#include <utility> // std::pair
#include <model/Utilities/StaticID.h>


namespace model
{
	
    /**************************************************************************/
	template <typename VertexType>
	class VertexInsertion
    {
		
        typedef typename std::map<size_t,VertexType> NetworkVertexMapType;

        //! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
	public:

		/* Constructor ******************************/
		VertexInsertion(NetworkVertexMapType& networkVertexMapRef_in) :
        /* init */ networkVertexMapRef(networkVertexMapRef_in)
        {}

        
        /* insert ***********************************/
        template <typename ...NodeArgTypes>
        //size_t insert(const NodeArgTypes&... NodeInput)
        std::pair<typename NetworkVertexMapType::iterator,bool> insert(const NodeArgTypes&... NodeInput)
        {/*!@tparam NodeArgTypes a variable number of input parameter types.
          * @param[in] NodeInput a variable number of input parameters.
          *
          * Constructs a VertexType vertex with constructor arguments NodeInput
          * and inserts it into the network.
          */
            // 1- Creates a new node
            //std::auto_ptr<VertexType> pN (new VertexType(NodeInput...) );
            // 2- Inserts the new node in NodeContainer
            const size_t nodeID(StaticID<VertexType>::nextID());
            
            const std::pair<typename NetworkVertexMapType::iterator,bool> inserted=networkVertexMapRef.emplace(std::piecewise_construct,
                                                                                                               std::make_tuple(nodeID),
                                                                                                               std::make_tuple(NodeInput...) );            
            assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
            assert(inserted.first->first == nodeID && "KEY != nodeID");
            assert(inserted.first->second.sID == nodeID && "sID != nodeID");
            
            // 3- Returns the static ID of the new node.
//            return nodeID;
            return inserted;
        }
        
//		/* insert ***********************************/
//		template <typename ...NodeArgTypes>
//		size_t insert(const NodeArgTypes&... NodeInput)
//        {/*!@tparam NodeArgTypes a variable number of input parameter types.
//          * @param[in] NodeInput a variable number of input parameters.
//          *
//          * Constructs a VertexType vertex with constructor arguments NodeInput
//          * and inserts it into the network.
//          */
//			// 1- Creates a new node
//			std::auto_ptr<VertexType> pN (new VertexType(NodeInput...) );
//			// 2- Inserts the new node in NodeContainer
//			const size_t nodeID(pN->sID);
//            const bool success(networkVertexMapRef.insert(nodeID , pN ).second);
//            assert(success && "CANNOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");			
//            // 3- Returns the static ID of the new node.
//			return nodeID;
//		}
		
	};
	
    /**************************************************************************/
} // namespace model
#endif
