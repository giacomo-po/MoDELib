/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXFINDER_H_
#define model_VERTEXFINDER_H_

#include <utility> // for std::pair
//#include <boost/ptr_container/ptr_map.hpp> // for boost::ptr_map
#include <map> // for boost::ptr_map

namespace model
{
    
    template <typename VertexType>
    class VertexFinder
    {
        
        //		typedef boost::ptr_map<size_t,VertexType> NetworkVertexMapType;
        typedef std::map<size_t,VertexType> NetworkVertexMapType;
        //! A reference to the network vertex map
        NetworkVertexMapType& networkVertexMapRef;
        
    public:
        //typedef size_t ReturnType;
        
        typedef std::pair<bool,VertexType* const>		 isNetworkVertexType;
        typedef std::pair<bool,const VertexType* const>	 isConstNetworkVertexType;
        
        
        /* Constructor ******************************/
        VertexFinder(NetworkVertexMapType& networkVertexMapRef_in) : networkVertexMapRef(networkVertexMapRef_in) {}
        
        /************************************************************/
        // node
        isNetworkVertexType node(const size_t & k)
        {/*!If a node with static ID = k exists in the network the function returns the pair (true, node pointer).
          * Otherwise returns the pair (false,NULL).
          */
            typename NetworkVertexMapType::iterator nodeIter(networkVertexMapRef.find(k));
            //			return (nodeIter!=networkVertexMapRef.end())? std::make_pair(true,nodeIter->second) : std::make_pair(false,(VertexType* const) NULL);
            return (nodeIter!=networkVertexMapRef.end())? std::make_pair(true,&nodeIter->second) : std::make_pair(false,(VertexType*) NULL);
        }
        
        isConstNetworkVertexType node(const size_t & k) const
        {
            typename NetworkVertexMapType::const_iterator nodeIter(networkVertexMapRef.find(k));
            return (nodeIter!=networkVertexMapRef.end())? std::make_pair(true,&nodeIter->second) : std::make_pair(false,(const VertexType* const) NULL);
        }
        
    };
    
    //////////////////////////////////////////////////////////////
} // namespace model
#endif
