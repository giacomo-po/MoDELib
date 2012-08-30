/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

public:

	/* Basic types defined by TypeTraits<Derived> ***************************************************/
	typedef typename TypeTraits<Derived>::NetworkType				NetworkType;			
	typedef typename TypeTraits<Derived>::SubNetworkType			SubNetworkType;
	typedef typename TypeTraits<Derived>::NodeType					NodeType;
	typedef typename TypeTraits<Derived>::LinkType					LinkType;
	typedef typename TypeTraits<Derived>::FlowType					FlowType;

	/* Containers and iterators used in Network *****************************************************/
	typedef typename boost::ptr_map<size_t,NodeType>				NetworkNodeContainerType;

	typedef typename VertexFinder<NodeType>::isNetworkVertexType		isNetworkNodeType;
	typedef typename VertexFinder<NodeType>::isConstNetworkVertexType	isConstNetworkNodeType;
	typedef typename EdgeFinder<LinkType>::isNetworkEdgeType			   isNetworkLinkType;
	typedef typename EdgeFinder<LinkType>::isConstNetworkEdgeType		isConstNetworkLinkType;



	typedef std::pair<NodeType* const,NodeType* const>			    NodePairType;
	typedef std::pair<size_t,size_t>								       LinkIDType;
	typedef boost::ptr_map<LinkIDType,LinkType>						 NetworkLinkContainerType;

	/* Containers and iterators used in SubNetwork **************************************************/
	typedef std::map<size_t,SubNetworkType* const>					SubNetworkContainerType;
	typedef std::map<size_t,NodeType* const>						SubNetworkNodeContainerType;
	typedef std::map<LinkIDType,LinkType* const>					SubNetworkLinkContainerType;

	/* Containers and iterators used in SubNetwork **************************************************/
	typedef boost::tuple<NodeType*,LinkType*,short int>				NeighborType;
	typedef std::map<size_t,NeighborType>							NeighborContainerType;