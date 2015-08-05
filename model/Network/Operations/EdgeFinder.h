/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EDGEFINDER_H_
#define model_EDGEFINDER_H_

#include <utility> // for std::pair
//#include <boost/ptr_container/ptr_map.hpp> // for boost::ptr_map
#include <map>

namespace model {
	
	template <typename EdgeType, bool forceConst=false>
	class EdgeFinder
    {

        typedef std::map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;

		//! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
	public:

		typedef std::pair<bool,EdgeType* const>			isNetworkEdgeType;
		typedef std::pair<bool,const EdgeType* const>	isConstNetworkEdgeType;
		
		
		/* Constructor ******************************/
		EdgeFinder(NetworkEdgeMapType& networkEdgeMapRef_in) : networkEdgeMapRef(networkEdgeMapRef_in) {}
		
		/************************************************************/
		// Edge
		isNetworkEdgeType link(const size_t & i, const size_t & j)
        {
			/*!If the Edge (i->j) exists in the network the function returns the pair (true, Edge pointer).
			 * Otherwise returns the pair (false,NULL). 
			 */
			typename NetworkEdgeMapType::iterator edgeIter(networkEdgeMapRef.find(std::make_pair(i,j))); 
//			return (edgeIter!=networkEdgeMapRef.end())? std::make_pair(true,edgeIter->second) : std::make_pair(false,(EdgeType* const) NULL);
			return (edgeIter!=networkEdgeMapRef.end())? std::make_pair(true,&edgeIter->second) : std::make_pair(false,(EdgeType*) NULL);
		}
		
		isConstNetworkEdgeType link(const size_t & i, const size_t & j) const
        {
			typename NetworkEdgeMapType::const_iterator edgeIter(networkEdgeMapRef.find(std::make_pair(i,j))); 
			return (edgeIter!=networkEdgeMapRef.end())? std::make_pair(true,&edgeIter->second) : std::make_pair(false,(const EdgeType* const) NULL);
		}
		
	};
	
	template <typename EdgeType>
	class EdgeFinder<EdgeType,true>
    {
	
        typedef std::map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;

        //! A reference to the network Edge map
		const NetworkEdgeMapType& networkEdgeMapRef;
		
	public:
		
		typedef std::pair<bool,EdgeType* const>			isNetworkEdgeType;
		typedef std::pair<bool,const EdgeType* const>	isConstNetworkEdgeType;
		
		
		/* Constructor ******************************/
		EdgeFinder(const NetworkEdgeMapType& networkEdgeMapRef_in) : networkEdgeMapRef(networkEdgeMapRef_in) {}
		
		/************************************************************/
		// Edge
//		isNetworkEdgeType link(const size_t & i, const size_t & j){
//			/*!If the Edge (i->j) exists in the network the function returns the pair (true, Edge pointer).
//			 * Otherwise returns the pair (false,NULL). 
//			 */
//			typename NetworkEdgeMapType::iterator edgeIter(networkEdgeMapRef.find(std::make_pair(i,j))); 
//			return (edgeIter!=networkEdgeMapRef.end())? std::make_pair(true,edgeIter->second) : std::make_pair(false,(EdgeType* const) NULL);
//		}
		
		isConstNetworkEdgeType link(const size_t & i, const size_t & j) const
        {
			typename NetworkEdgeMapType::const_iterator edgeIter(networkEdgeMapRef.find(std::make_pair(i,j))); 
			return (edgeIter!=networkEdgeMapRef.end())? std::make_pair(true,&edgeIter->second) : std::make_pair(false,(const EdgeType* const) NULL);
		}
		
	};
	
	//////////////////////////////////////////////////////////////
} // namespace model
#endif
