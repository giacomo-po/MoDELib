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
#include <map> // for boost::ptr_map

namespace model
{
    
    template <typename VertexType>
    class VertexFinder
    {
        
        typedef std::map<size_t,VertexType> NetworkVertexMapType;
        
        //! A reference to the network vertex map
        NetworkVertexMapType& networkVertexMapRef;
        
    public:
        
        typedef std::pair<bool,VertexType* const>		 isNetworkVertexType;
        typedef std::pair<bool,const VertexType* const>	 isConstNetworkVertexType;
        
        
        /**********************************************************************/
        VertexFinder(NetworkVertexMapType& networkVertexMapRef_in) :
        /* init */ networkVertexMapRef(networkVertexMapRef_in)
        {
        
        }
        
        /**********************************************************************/
        isNetworkVertexType node(const size_t & k)
        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
            typename NetworkVertexMapType::iterator nodeIter(networkVertexMapRef.find(k));
            return (nodeIter!=networkVertexMapRef.end())? std::make_pair(true,&nodeIter->second) : std::make_pair(false,(VertexType*) NULL);
        }
        
        isConstNetworkVertexType node(const size_t & k) const
        {/*!\returns A <bool,const NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
            typename NetworkVertexMapType::const_iterator nodeIter(networkVertexMapRef.find(k));
            return (nodeIter!=networkVertexMapRef.end())? std::make_pair(true,&nodeIter->second) : std::make_pair(false,(const VertexType* const) NULL);
        }
        
    };
    
} // namespace model
#endif
