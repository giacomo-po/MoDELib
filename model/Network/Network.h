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

//#include <boost/ptr_container/ptr_map.hpp>


#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/CRTP.h>
#include <model/Utilities/AddressBook.h>
#include <model/Utilities/NonCopyable.h>

#include <model/Network/NetworkComponent.h>
#include <model/Network/Operations/includeNetworkOperations.h>
#include <model/Network/Algorithms/ParallelExecute.h>

#include <model/Utilities/SequentialBinFile.h>


namespace model
{
		
	template <typename Derived>
    class Network : public NonCopyable,
    /*          */  protected std::map<size_t,typename TypeTraits<Derived>::NodeType>,
    /*          */  protected std::map<std::pair<size_t,size_t>,typename TypeTraits<Derived>::LinkType>,
	/*          */  public  CRTP<Derived>,
	/*          */  public  AddressBook<NetworkComponent<typename TypeTraits<Derived>::NodeType,typename TypeTraits<Derived>::LinkType>,0>
    {
		
		/*!
		 *		\code
		 *		namespace model
         *      {
		 *			template<>
		 *			struct TypeTraits<Derived>
         *          {
		 *				typedef SomeNodeType		NodeType;
		 *				typedef SomeLinkType		LinkType;
		 *				typedef SomeFlowType		FlowType;
		 *			};
		 *		}
		 *		\endcode
		 */
		
#include <model/Network/NetworkTypedefs.h>
		
        
		
	public:
		
        /**********************************************************************/
		size_t nodeOrder() const
        {/*!\returns The number of vertices in the Network
          */
			return NetworkNodeContainerType::size();
		}
		
        /**********************************************************************/
		size_t linkOrder() const
        {/*!\returns The number of edges in the Network
          */
			return NetworkLinkContainerType::size();
		}
        
        /**********************************************************************/
        const NetworkNodeContainerType& nodes() const
        {/*!\returns A const reference to the node container
          */
            return *this;
        }
        
        /**********************************************************************/
        NetworkNodeContainerType& nodes()
        {/*!\returns A reference to the node container
          */
            return *this;
        }
        
        /**********************************************************************/
        const NetworkLinkContainerType& links() const
        {/*!\returns A const reference to the link container
          */
            return *this;
        }
        
        /**********************************************************************/
        NetworkLinkContainerType& links()
        {/*!\returns A reference to the link container
          */
            return *this;
        }
        
        /**********************************************************************/
		isNetworkNodeType node(const size_t & k)
        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer 
          * the the node
          */
			return VertexFinder<NodeType>(*this).node(k);
		}
		
		isConstNetworkNodeType node(const size_t & k) const
        {/*!\returns A <bool,const NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
			return VertexFinder<NodeType>(*this).node(k);
		}
		
        /**********************************************************************/
		isNetworkLinkType link(const size_t& i, const size_t& j)
        {
			return EdgeFinder<LinkType>(*this).link(i,j);
		}
        
		isConstNetworkLinkType link(const size_t& i, const size_t& j) const
        {
			return EdgeFinder<LinkType,true>(*this).link(i,j);
		}
		
        /**********************************************************************/
		typename NetworkNodeContainerType::iterator nodeBegin()
        {//!\returns An iterator to the first node in the network.
			return NetworkNodeContainerType::begin();
		}
		
		typename NetworkNodeContainerType::const_iterator nodeBegin() const
        {//!\returns A const iterator to the first node in the network.
			return NetworkNodeContainerType::begin();
		}
		
        /**********************************************************************/
		typename NetworkNodeContainerType::iterator nodeEnd()
        {/*! @param[out] An iterator referring to the past-the-end vertex in
          *  the network.
          */
			return NetworkNodeContainerType::end();
		}
		
		typename NetworkNodeContainerType::const_iterator nodeEnd() const
        {/*! @param[out] A const iterator referring to the past-the-end vertex
          *  in the network.
          */
			return NetworkNodeContainerType::end();
		}
		
        /**********************************************************************/
		typename NetworkLinkContainerType::iterator linkBegin()
        {
			return NetworkLinkContainerType::begin();
		}
		
		typename NetworkLinkContainerType::const_iterator linkBegin() const
        {
			return NetworkLinkContainerType::begin();
		}
		
        /**********************************************************************/
		typename NetworkLinkContainerType::iterator linkEnd()
        {
			return NetworkLinkContainerType::end();
		}
		
		typename NetworkLinkContainerType::const_iterator linkEnd() const
        {
			return NetworkLinkContainerType::end();
		}
        
        /**********************************************************************/
        typename NetworkComponentContainerType::iterator componentBegin()
        {
            return this->ABbegin();
        }
        
        /**********************************************************************/
        typename NetworkComponentContainerType::iterator componentEnd()
        {
            return this->ABend();
        }
        
        /**********************************************************************/
		template <typename ...NodeArgTypes>
        std::pair<typename NetworkNodeContainerType::iterator,bool> insertVertex(const NodeArgTypes&... nodeInput)
        {/*! @param[in] nodeInput
          *  Inserts a new vertex in the Network using nodeInput as variable
          *  constructor arguments
          */
			return VertexInsertion<NodeType>(*this).insert(nodeInput...);
		}
		
        /**********************************************************************/
        template <typename ...EdgeArgTypes>
        bool connect(const size_t& i, const size_t& j, const EdgeArgTypes&... edgeArgs)
        {/*! @param[in] i the StaticID of the source vertex
          *  @param[in] j the StaticID of the sink vertex
          *  @param[in] f the flow
          *  Connects source vertex i to sink vertex j with flow f
          */
            return VertexConnection<NodeType,LinkType>(*this,*this).connect(i,j,edgeArgs...);
		}
		
        /**********************************************************************/
		template<bool removeIsolatedNodes>
		bool disconnect(const size_t& i, const size_t& j)
        {/*! @param[in] i the StaticID of the source vertex
          *  @param[in] j the StaticID of the sink vertex
          *  Disconnects the edge i->j
          */
			return VertexConnection<NodeType,LinkType>(*this,*this).template disconnect<removeIsolatedNodes>(i,j);
		}
		
        /**********************************************************************/
		template<bool removeIsolatedNodes>
		bool removeVertex(const size_t& k)
        {/*! @param[in] k the StaticID of the vertex
          *  Removes the k-th Vertex from the Network
          */
			return VertexConnection<NodeType,LinkType>(*this,*this).template remove<removeIsolatedNodes>(k);
		}
		
        /**********************************************************************/
		template<bool removeIsolatedNodes>
		size_t disconnect_if(bool (LinkType::*Lfptr)(void) const)
        {
			return VertexConnection<NodeType,LinkType>(*this,*this).template disconnect_if<removeIsolatedNodes>(Lfptr);
		}
		
        /**********************************************************************/
		template<bool removeIsolatedNodes, typename T>
		size_t disconnect_if(bool (LinkType::*Lfptr)(const T &) const, const T & input)
        {
			return VertexConnection<NodeType,LinkType>(*this,*this).template disconnect_if<removeIsolatedNodes,T>(Lfptr,input);
		}
		
        /**********************************************************************/
		template <typename ...NodeArgTypes>
		std::pair<typename NetworkNodeContainerType::iterator,bool> expand(const size_t & i, const size_t & j, const NodeArgTypes&... Args)
        {
			return EdgeExpansion<NodeType,LinkType>(*this,*this).expand(i,j,Args...);
		}
		
		/* multiExpand **************************************************/ // CLEAN THIS
		template <typename T>
		std::map<T,size_t>  multiExpand(const size_t& i, const size_t& j, const std::map<double,T>& expandMap)
        {
			
            
			
			isNetworkLinkType Lij(this->link(i,j));
			assert(Lij.first);
			
			// REMOVE THIS SECTION FROM NETWORK LAYER
			enum {dim=3};
			typedef Eigen::Matrix<double,3,1> VectorDim;
			std::map<double,VectorDim> pointMap;
			for (typename std::map<double,T>::const_iterator iter=expandMap.begin(); iter!=expandMap.end();++iter)
            {
				pointMap.insert(std::make_pair(iter->first,Lij.second->get_r(iter->first)));
			}
			
			std::map<T,size_t> temp;
			typename std::map<double,T>::const_iterator iterEx=expandMap.begin();
			std::pair<bool,size_t> currNode = std::make_pair(true,i); // initialize currNode with i
			for (typename std::map<double,VectorDim>::const_iterator iter=pointMap.begin(); iter!=pointMap.end();++iter)
            {
				currNode = expand(currNode.second,j,iter->second);
				assert(currNode.first);
				temp.insert(std::make_pair(iterEx->second,currNode.second));
				iterEx++;
			}
			
			assert(temp.size()==expandMap.size());
			
			
			return temp;
		}
		
        /**********************************************************************/
		template <typename ...NodeArgTypes>
		void contract(const size_t& i, const size_t& j, const NodeArgTypes&... nodeInput)
        {
			VertexContraction<NodeType,LinkType>(*this,*this).contract(i,j,nodeInput...);
		}
		
        /**********************************************************************/
		void contractSecond(const size_t& i, const size_t& j)
        {
			VertexContraction<NodeType,LinkType>(*this,*this).contractSecond(i,j);
		}
        
        /**********************************************************************/
		void parallelExecute(void (LinkType::*Lfptr)(void))
        {
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Lfptr);
		}
        
        /**********************************************************************/
		void parallelExecute(void (NodeType::*Vfptr)(void))
        {
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Vfptr);
		}
        
        /**********************************************************************/
		template <typename T>
		void parallelExecute(void (NodeType::*Vfptr)(const T&), const T& input)
        {
			ParallelExecute<NodeType,LinkType>(*this,*this).execute(Vfptr,input);
		}
        
        /**********************************************************************/
		template <class T>
		friend T& operator << (T& os, const NetworkNodeContainerType& nnC)
        {
			for (typename NetworkNodeContainerType::const_iterator nodeIter=nnC.begin();nodeIter!=nnC.end();++nodeIter)
            {
				os << (*nodeIter->second) << "\n";
            }
            return os;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const NetworkLinkContainerType& nlC)
        {
            for (typename NetworkLinkContainerType::const_iterator linkIter=nlC.begin();linkIter!=nlC.end();++linkIter)
            {
                os << (*linkIter->second) << "\n";
            }
            return os;
        }
        
        
    };	// end Network
    
} // namespace model
#endif
