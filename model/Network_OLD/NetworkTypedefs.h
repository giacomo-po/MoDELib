/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

public:

	typedef typename TypeTraits<Derived>::NetworkType				NetworkType;			
    typedef typename TypeTraits<Derived>::NodeType					NodeType;
	typedef typename TypeTraits<Derived>::LinkType					LinkType;
	typedef typename TypeTraits<Derived>::FlowType					FlowType;
    typedef NetworkComponent<NodeType,LinkType>                     NetworkComponentType;


    typedef typename std::map<size_t,NodeType>                          NetworkNodeContainerType;
    typedef typename VertexFinder<NodeType>::isNetworkVertexType		isNetworkNodeType;
	typedef typename VertexFinder<NodeType>::isConstNetworkVertexType	isConstNetworkNodeType;
	typedef typename EdgeFinder<LinkType>::isNetworkEdgeType            isNetworkLinkType;
	typedef typename EdgeFinder<LinkType>::isConstNetworkEdgeType		isConstNetworkLinkType;



	typedef std::pair<NodeType* const,NodeType* const>			    NodePairType;
	typedef std::pair<size_t,size_t>								LinkIDType;
    typedef std::map<LinkIDType,LinkType>                     NetworkLinkContainerType;


	typedef std::map<size_t,NetworkComponentType* const>			NetworkComponentContainerType;
    typedef typename NetworkComponent<NodeType,LinkType>::NetworkComponentNodeContainerType SubNetworkNodeContainerType;
    typedef typename NetworkComponent<NodeType,LinkType>::NetworkComponentLinkContainerType SubNetworkLinkContainerType;


	typedef std::tuple<NodeType* const ,LinkType* const,short int>				NeighborType;
    typedef std::map<size_t,NeighborType>						    	NeighborContainerType;
