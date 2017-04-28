/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNetwork_H_
#define model_LoopNetwork_H_

#include <iostream>
#include <assert.h>
#include <memory>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <utility>
#include <algorithm>

#include <model/LoopNetwork/LoopNode.h>
//#include <model/LoopNetwork/NetworkNode.h>
#include <model/LoopNetwork/LoopLink.h>
#include <model/LoopNetwork/NetworkLink.h>
#include <model/LoopNetwork/NetworkLinkObserver.h>
#include <model/LoopNetwork/Loop.h>
#include <model/LoopNetwork/LoopObserver.h>

namespace model
{
    template<typename Derived>
    class LoopNetwork :
//    /*               */ private std::map<size_t,typename TypeTraits<Derived>::NodeType>,
    /*               */ private std::map<size_t,std::shared_ptr<typename TypeTraits<Derived>::NodeType>>,
    /*               */ private std::multimap<std::pair<size_t,size_t>,LoopLink<typename TypeTraits<Derived>::LinkType>>,
    /*               */ public NetworkLinkObserver<typename TypeTraits<Derived>::LinkType>,
    /*               */ public NodeObserver<typename TypeTraits<Derived>::NodeType>,
    /*               */ private std::map<size_t,const typename TypeTraits<Derived>::LoopType* const>
    {
        typedef typename TypeTraits<Derived>::NodeType NodeType;
//        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<typename TypeTraits<Derived>::LinkType> LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        
        //        typedef std::map<size_t,NodeType>     LoopNodeContainerType;

        typedef std::multimap<std::pair<size_t,size_t>,LoopLinkType> LoopLinkContainerType;
        typedef std::map<size_t,const LoopType* const> LoopContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename NetworkLinkObserverType::IsConstNetworkLinkType IsConstNetworkLinkType;
        typedef LoopObserver<LoopType> LoopObserverType;
        
        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef std::map<size_t,SharedNodePtrType>     LoopNodeContainerType;

        typedef std::pair<bool,SharedNodePtrType>          IsNodeType;
        typedef std::pair<bool,const std::shared_ptr<const NodeType>>	 IsConstNodeType;
        
        typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>> ExpandTupleType;

        
        /**********************************************************************/
        LoopNodeContainerType& danglingNodes()
        {
            return *this;
        }
        
        //        /**********************************************************************/
        //        const LoopNodeContainerType& danglingNodes() const
        //        {
        //            return *this;
        //        }
        
        /**********************************************************************/
        IsNodeType danglingNode(const size_t & k)
        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
            //            std::cout<<"finding "<<k<<std::endl;
            typename LoopNodeContainerType::iterator nodeIter(danglingNodes().find(k));
            //            return (nodeIter!=danglingNodes().end())? std::make_pair(true,&nodeIter->second) : std::make_pair(false,(NodeType*) NULL);
            return (nodeIter!=danglingNodes().end())? std::make_pair(true,nodeIter->second) : std::make_pair(false,SharedNodePtrType(nullptr));
        }
        
        /**********************************************************************/
        Derived&   derived()
        {//!\returns A reference to the derived object (CRTP)
            return *static_cast<Derived*>(this);
        }
        
        const Derived&   derived() const
        {//!\returns A  reference to the const derived object (CRTP)
            return *static_cast<const Derived*>(this);
        }
        
        /**********************************************************************/
        void connect(const SharedNodePtrType& n0,
                     const SharedNodePtrType& n1,
                     const std::shared_ptr<LoopType>& tempLoop)
        {
            
            assert(n0->sID!=n1->sID && "Cannot connect a node to itself");
            
            loopLinks().emplace(std::piecewise_construct,
                                std::make_tuple(n0->sID,n1->sID),
                                std::make_tuple(n0,n1, tempLoop) );
        }
        

        
    public:



        
        
        /**********************************************************************/
        void clearDanglingNodes()
        {
            danglingNodes().clear();
        }

        /**********************************************************************/
        LoopLinkContainerType& loopLinks()
        {
            return *this;
        }

        
        /**********************************************************************/
        const LoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::shared_ptr<LinkType> pLink(const SharedNodePtrType& nI, const SharedNodePtrType& nJ) const
        {
            typename LoopLinkContainerType::const_iterator iterIJ(loopLinks().find(std::make_pair(nI->sID,nJ->sID)));
            if(iterIJ!=loopLinks().end())
            {
                return iterIJ->second.pLink;
            }
            else
            {
                typename LoopLinkContainerType::const_iterator iterJI(loopLinks().find(std::make_pair(nJ->sID,nI->sID)));
                if(iterJI!=loopLinks().end())
                {
                    return iterJI->second.pLink;
                }
                else
                {
                    return nI->sID<nJ->sID? std::make_shared<LinkType>(nI,nJ) : std::make_shared<LinkType>(nJ,nI);
                }
            }
        }
        
        /**********************************************************************/
        const LoopContainerType& loops() const
        {//!\returns the loop container
            return *this;
        }
        

        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        std::pair<typename LoopNodeContainerType::iterator,bool> insertDanglingNode(const NodeArgTypes&... nodeInput)
        {/*! @param[in] nodeInput
          *  Inserts a new vertex in the Network using nodeInput as variable
          *  constructor arguments
          */
            const size_t nodeID(StaticID<NodeType>::nextID());

            const std::pair<typename LoopNodeContainerType::iterator,bool> inserted=danglingNodes().emplace(std::piecewise_construct,
                                                                                                    std::make_tuple(nodeID),
                                                                                                    std::make_tuple(new NodeType(nodeInput...)) );

            assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
            assert(inserted.first->first == nodeID && "KEY != nodeID");
            assert(inserted.first->second->sID == nodeID && "sID != nodeID");
            return inserted;
            
        }
                
        /**********************************************************************/
        template <typename ...LoopArgTypes>
        void insertLoop(const std::vector<size_t> nodeIDs,const LoopArgTypes&... loopInput)
        {
            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(derived(),loopInput...);
            
            for(size_t k=0;k<nodeIDs.size();++k)
            {
                const size_t next=k+1<nodeIDs.size()? k+1 : 0;
                
                IsNodeType n0=danglingNode(nodeIDs[k]);
                IsNodeType n1=danglingNode(nodeIDs[next]);
                assert(n0.first && "Node not found");
                assert(n1.first && "Node not found");
                
                connect(n0.second,n1.second,tempLoop);
            }
            
            assert(tempLoop->isLoop() && "Not a loop.");
        }
        
        

        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        SharedNodePtrType expand(const size_t& a, const size_t& b, const NodeArgTypes&... Args)
        {
            
            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
            
            const size_t i=std::min(a,b);
            const size_t j=std::max(a,b);
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Expanding non-existing link");
            
            std::cout<<"Expanding "<<i<<","<<j<<std::endl;
            
            // Create new node
            SharedNodePtrType newNode=SharedNodePtrType(new NodeType(Args...));
//            const size_t newID=newNode->sID;
            
            // Store what needs to be connected (i,new,j,Loop), or (j,new,i,Loop). This also holds temporarily disconnected nodes
            std::deque<ExpandTupleType> expandDeq;
            for(const auto& llink : Lij.second->loopLinks())
            {
                expandDeq.emplace_back(llink->source,llink->sink,llink->pLoop);
            }
            
            // Delete all LoopLinks of type i->j or j->i
            loopLinks().erase(std::make_pair(i,j));
            loopLinks().erase(std::make_pair(j,i));
            
            for (const auto& tup : expandDeq)
            {
//                std::cout<<std::get<0>(tup)->sID<<"->"<<std::get<1>(tup)->sID<<", loop "<<std::get<2>(tup)->sID<<std::endl;
                
                connect(std::get<0>(tup),newNode, std::get<2>(tup));
//                loopLinks().emplace(std::piecewise_construct,
//                                    std::make_tuple(std::get<0>(tup)->sID,newID),
//                                    std::make_tuple(std::get<0>(tup),newNode, std::get<2>(tup)) );

                                connect(newNode,std::get<1>(tup), std::get<2>(tup));
//                loopLinks().emplace(std::piecewise_construct,
//                                    std::make_tuple(newID,std::get<1>(tup)->sID),
//                                    std::make_tuple(newNode,std::get<1>(tup), std::get<2>(tup)) );
            }
            
            return newNode;
        }
        
        /**********************************************************************/
        void contract(const size_t& a, const size_t& b)
        {/*! Performs a contract operation in which the nodes of a NetworkLink
          * (a,b) are merged. Node a survives the merge, while node b is 
          * destroyed
          */
            
            const size_t i=std::min(a,b);
            const size_t j=std::max(a,b);
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Contracting non-existing link");
            
            
            if(i==a && j==b) // j is destroyed
            {
                for(const auto& loopLink:Lij.second->loopLinks())
                {
                
                }
            }
            else // i is destroyed
            {
            
            
            }
            
//            FINISH HERE
//            typename LinkType::LoopLinkContainerType loopLinksIJ=Lij.second->loopLinks();
//            for(const auto& llink : loopLinksIJ)
//            {
//                if(llink->source->sID==i && llink->sink->sID==j) // same direction
//                {
//                    loopLinks().emplace(std::piecewise_construct,
//                                        std::make_tuple(i,newID),
//                                        std::make_tuple(llink->source,&newNode.first->second, llink->pLoop) );
//                    
//                    loopLinks().emplace(std::piecewise_construct,
//                                        std::make_tuple(newID,j),
//                                        std::make_tuple(&newNode.first->second,llink->sink,llink->pLoop) );
//                }
//                else if(llink->source->sID==j && llink->sink->sID==i) // opposite direction
//                {
//                    loopLinks().emplace(std::piecewise_construct,
//                                        std::make_tuple(j,newID),
//                                        std::make_tuple(llink->source,&newNode.first->second, llink->pLoop) );
//                    
//                    loopLinks().emplace(std::piecewise_construct,
//                                        std::make_tuple(newID,i),
//                                        std::make_tuple(&newNode.first->second,llink->sink,llink->pLoop) );
//                }
//                else
//                {
//                    assert(0 && "LoopLink must connect i->j or viceversa");
//                }
//            }
//            
//            loopLinks().erase(std::make_pair(i,j));
//            loopLinks().erase(std::make_pair(j,i));
            
        }
        
        /**********************************************************************/
        void checkLoops() const
        {
            std::cout<<"Checking "<<LoopObserverType::loops().size()<<" loops"<<std::endl;
            for(const auto& loop : LoopObserverType::loops())
            {
                assert(loop.second->isLoop() && "note a closed loop");
            }
        
        }
        
        void printLoopLinks() const
        {
            for(const auto& link : loopLinks())
            {
                std::cout<<"link "<<link.second.source->sID<<"->"<<link.second.sink->sID
                <<" (prev "<<link.second.prev->source->sID<<"->"<<link.second.prev->sink->sID<<")"
                <<" (next "<<link.second.next->source->sID<<"->"<<link.second.next->sink->sID<<")"
                <<std::endl;
            }
        }
        
        void printNodes() const
        {
            for(const auto& node : this->nodes())
            {
                std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
            }
        }
    };
    
    
}
#endif
