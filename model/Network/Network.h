/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NETWORK_H_
#define model_NETWORK_H_


#include <assert.h>

//#include <map>
#include <boost/ptr_container/ptr_map.hpp>
//#include <boost/smart_ptr/shared_ptr.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/utility.hpp>
//#include <memory> // std::shared_ptr

#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/AddressBook.h>

#include <model/Network/Operations/includeNetworkOperations.h>
#include <model/Network/Algorithms/ParallelExecute.h>

#include <model/Utilities/SequentialBinFile.h>


namespace model {
	
	
	template <typename Derived>
	class Network : boost::noncopyable,
	/*          */  protected boost::ptr_map<size_t,typename TypeTraits<Derived>::NodeType>,
	/*          */  protected boost::ptr_map<std::pair<size_t,size_t>,typename TypeTraits<Derived>::LinkType>,
//	/*          */  protected std::map<size_t,std::auto_ptr<typename TypeTraits<Derived>::NodeType> >,
//	/*          */  protected std::map<std::pair<size_t,size_t>,std::auto_ptr<typename TypeTraits<Derived>::LinkType> >,
	/*          */  public  CRTP<Derived>,
	/*          */  public  AddressBook<typename TypeTraits<Derived>::SubNetworkType,0>{
		
		/*!
		 *		\code
		 *		namespace model{
		 *			template<>
		 *			struct TypeTraits<Derived>{
		 *				typedef SomeSubNetworkType	NodeType;
		 *				typedef SomeNodeType		NodeType;
		 *				typedef SomeLinkType		LinkType;
		 *				typedef SomeFlowType		FlowType;
		 *			};	
		 *		}
		 *		\endcode
		 */
		
#include <model/Network/NetworkTypedefs.h>
		
		
		
		
		typedef std::map<size_t,SubNetworkType* const> AddressMapType;
		typedef typename AddressMapType::iterator AddressMapIteratorType;

	private:
		

		
		
	public:
		
		/************************************************************/
		// nodeOrder
		size_t nodeOrder() const {
			//! The number of nodes in the Network
			return NetworkNodeContainerType::size();
		}
		
		/************************************************************/
		// linkOrder
		size_t linkOrder() const {
			//! The number of edges in the Network
			return NetworkLinkContainerType::size();
		}

		/************************************************************/
		// node
		isNetworkNodeType node(const size_t & k)  {
			return VertexFinder<NodeType>(*this).node(k);
		}
		
		isConstNetworkNodeType node(const size_t & k) const {
			return VertexFinder<NodeType>(*this).node(k);
		}
		
		/************************************************************/
		// link
		isNetworkLinkType link(const size_t& i, const size_t& j){
			return EdgeFinder<LinkType>(*this).link(i,j);
		}
		
//		isConstNetworkLinkType link(const size_t& i, const size_t& j) const {
//			return EdgeFinder<LinkType>(*this).link(i,j);
//		}

		isConstNetworkLinkType link(const size_t& i, const size_t& j) const {
			return EdgeFinder<LinkType,true>(*this).link(i,j);
		}
		
		
		/************************************************************/
		// nodeBegin
		typename NetworkNodeContainerType::iterator nodeBegin() {
			//! An iterator to the first node in the network. 
			return NetworkNodeContainerType::begin();
		}
		
		typename NetworkNodeContainerType::const_iterator nodeBegin() const {
			//! A const iterator to the node in the network. 
			return NetworkNodeContainerType::begin();
		}
		
		/************************************************************/
		// nodeEnd
		typename NetworkNodeContainerType::iterator nodeEnd() {
			//! An iterator to the first link in the network. 
			return NetworkNodeContainerType::end();
		}
		
		typename NetworkNodeContainerType::const_iterator nodeEnd() const {
			//! An const iterator to the first link in the network. 
			return NetworkNodeContainerType::end();
		}
		
		/************************************************************/
		// linkBegin
		typename NetworkLinkContainerType::iterator linkBegin()  {
			return NetworkLinkContainerType::begin();
		}
		
		typename NetworkLinkContainerType::const_iterator linkBegin() const {
			return NetworkLinkContainerType::begin();
		}
		
		/************************************************************/
		// linkEnd
		typename NetworkLinkContainerType::iterator linkEnd() {
			return NetworkLinkContainerType::end();
		}
		
		typename NetworkLinkContainerType::const_iterator linkEnd() const {
			return NetworkLinkContainerType::end();
		}

		/* insert (a new vertex) **************************************/
		template <typename ...NodeArgTypes>
		size_t insert(const NodeArgTypes&... NodeInput){
			return VertexInsertion<NodeType>(*this).insert(NodeInput...);
		}
		
		/* connect ***************************************************/
		bool connect(const size_t& i, const size_t& j, const FlowType& f){
			return VertexConnection<NodeType,LinkType>(*this,*this).connect(i,j,f);
		}
		
		/* disconnect ************************************************/ 
		template<bool removeIsolatedNodes>
		bool disconnect(const size_t& i, const size_t& j){
			return VertexConnection<NodeType,LinkType>(*this,*this).disconnect<removeIsolatedNodes>(i,j);
		}
		
		/************************************************************/
		// remove (a node)
		template<bool removeIsolatedNodes>
		bool remove(const size_t& k){		
			return VertexConnection<NodeType,LinkType>(*this,*this).remove<removeIsolatedNodes>(k);
		}
		
		/************************************************************/
		// disconnect_if
		template<bool removeIsolatedNodes>
		size_t disconnect_if(bool (LinkType::*Lfptr)(void) const){
			return VertexConnection<NodeType,LinkType>(*this,*this).disconnect_if<removeIsolatedNodes>(Lfptr);
		}
		
		/************************************************************/
		// disconnect_if
		template<bool removeIsolatedNodes, typename T>
		size_t disconnect_if(bool (LinkType::*Lfptr)(const T &) const, const T & input){
			return VertexConnection<NodeType,LinkType>(*this,*this).disconnect_if<removeIsolatedNodes,T>(Lfptr,input);
		}
		
		/* expand ***************************************************/
		template <typename ...NodeArgTypes>
		std::pair<bool,size_t> expand(const size_t & i, const size_t & j, const NodeArgTypes&... Args){
			return EdgeExpansion<NodeType,LinkType>(*this,*this).expand(i,j,Args...);
		}
		
		
		/* multiExpand **************************************************/ // CLEAN THIS
		template <typename T>
		std::map<T,size_t>  multiExpand(const size_t& i, const size_t& j, const std::map<double,T>& expandMap){
			

			
			isNetworkLinkType Lij(this->link(i,j));
			assert(Lij.first);
			
			// REMOVE THIS SECTION FROM NETWORK LAYER
			enum {dim=3};
			typedef Eigen::Matrix<double,3,1> VectorDim;
			std::map<double,VectorDim> pointMap;
			for (typename std::map<double,T>::const_iterator iter=expandMap.begin(); iter!=expandMap.end();++iter){				
				pointMap.insert(std::make_pair(iter->first,Lij.second->get_r(iter->first)));
			}
			
			std::map<T,size_t> temp;
			typename std::map<double,T>::const_iterator iterEx=expandMap.begin();
			std::pair<bool,size_t> currNode = std::make_pair(true,i); // initialize currNode with i
			for (typename std::map<double,VectorDim>::const_iterator iter=pointMap.begin(); iter!=pointMap.end();++iter){				
				currNode = expand(currNode.second,j,iter->second);
				assert(currNode.first);
				temp.insert(std::make_pair(iterEx->second,currNode.second));
				iterEx++;
			}
			
			assert(temp.size()==expandMap.size());
			
			
			return temp;	
		}

		
		/* contract **************************************************/
		template <typename T>
		void contract(const size_t& i, const size_t& j, const T& NodeInput){
			VertexContraction<NodeType,LinkType>(*this,*this).contract(i,j,NodeInput);
		}
		

		/* contractSecond ********************************************/
		void contractSecond(const size_t& i, const size_t& j){
			VertexContraction<NodeType,LinkType>(*this,*this).contractSecond(i,j);
		} 


		/*************************************************************/
		void parallelExecute(void (LinkType::*Lfptr)(void)){
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Lfptr);
		}		

		/*************************************************************/
		void parallelExecute(void (NodeType::*Vfptr)(void)){
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Vfptr);
		}

		/*************************************************************/
		template <typename T>
		void parallelExecute(void (NodeType::*Vfptr)(const T&), const T& input)
        {
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Vfptr,input);
		}		
        
		/*************************************************************/
		// friend T& operator <<
		template <class T>
		friend T& operator << (T& os, const NetworkNodeContainerType& nnC)
        {
			for (typename NetworkNodeContainerType::const_iterator nodeIter=nnC.begin();nodeIter!=nnC.end();++nodeIter){				
				os << (*nodeIter->second) << "\n";
			}
			return os;
		}
		
		/*************************************************************/
		// friend T& operator <<
		template <class T>
		friend T& operator << (T& os, const NetworkLinkContainerType& nlC)
        {
			for (typename NetworkLinkContainerType::const_iterator linkIter=nlC.begin();linkIter!=nlC.end();++linkIter){				
				os << (*linkIter->second) << "\n";
			}
			return os;
		}
        

	};	// end Network
	/************************************************************/
	/************************************************************/
} // namespace model
#endif
