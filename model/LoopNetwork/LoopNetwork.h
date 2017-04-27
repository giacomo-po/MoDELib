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
        
        typedef std::map<size_t,std::shared_ptr<NodeType>>     LoopNodeContainerType;
        //        typedef std::map<size_t,NodeType>     LoopNodeContainerType;

        typedef std::multimap<std::pair<size_t,size_t>,LoopLinkType> LoopLinkContainerType;
        typedef std::map<size_t,const LoopType* const> LoopContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename NetworkLinkObserverType::IsConstNetworkLinkType IsConstNetworkLinkType;
        typedef LoopObserver<LoopType> LoopObserverType;
        
        typedef std::pair<bool,std::shared_ptr<NodeType> const>          IsNodeType;
        typedef std::pair<bool,const std::shared_ptr<const NodeType>>	 IsConstNodeType;
        
        
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
            return (nodeIter!=danglingNodes().end())? std::make_pair(true,nodeIter->second) : std::make_pair(false,std::shared_ptr<NodeType>(nullptr));
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
        std::shared_ptr<LinkType> pLink(const NodeType* const nI, const NodeType* const nJ) const
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
        LoopContainerType& loops()
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
//            const std::pair<typename LoopNodeContainerType::iterator,bool> inserted=danglingNodes().emplace(std::piecewise_construct,
//                                                                                                    std::make_tuple(nodeID),
//                                                                                                    std::make_tuple(nodeInput...) );

            const std::pair<typename LoopNodeContainerType::iterator,bool> inserted=danglingNodes().emplace(std::piecewise_construct,
                                                                                                    std::make_tuple(nodeID),
                                                                                                    std::make_tuple(new NodeType(nodeInput...)) );

            assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
            assert(inserted.first->first == nodeID && "KEY != nodeID");
            assert(inserted.first->second->sID == nodeID && "sID != nodeID");
            return inserted;
            
        }
        
//        'model::LoopLink<model::Dlink>::LoopLink(std::shared_ptr<const model::Dnode>, std::shared_ptr<const model::Dnode>, std::shared_ptr<model::Dloop>)'
        
        /**********************************************************************/
        template <typename ...LoopArgTypes>
        void insertLoop(const std::vector<size_t> nodeIDs,const LoopArgTypes&... loopInput)
        {
            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(derived(),loopInput...);
            
            for(size_t k=0;k<nodeIDs.size();++k)
            {
                const size_t next=k+1<nodeIDs.size()? k+1 : 0;
                
//                std::cout<<nodeIDs[k]<<" "<<nodeIDs[next]<<std::endl;

                IsNodeType n0=danglingNode(nodeIDs[k]);
                IsNodeType n1=danglingNode(nodeIDs[next]);
                
                
                connect(n0,n1,tempLoop);
//                std::cout<<"done connecting"<<std::endl;
//                loopLinks().emplace(std::piecewise_construct,
//                                    std::make_tuple(n0.second->sID,n1.second->sID),
//                                    std::make_tuple(n0.second,n1.second, tempLoop) );
//            
//                danglingNodes().erase(n0.second->sID); // remove temporary shared_ptr
//                danglingNodes().erase(n1.second->sID); // remove temporary shared_ptr
            }
            
            assert(tempLoop->isLoop() && "Not a loop.");
        }
        
        
        /**********************************************************************/
        void connect(const IsNodeType& n0,
                     const IsNodeType& n1,
                     const std::shared_ptr<LoopType>& tempLoop)
        {
            assert(n0.first && "Node not found");
            assert(n1.first && "Node not found");

            loopLinks().emplace(std::piecewise_construct,
                                std::make_tuple(n0.second->sID,n1.second->sID),
                                std::make_tuple(n0.second,n1.second, tempLoop) );
            
//            danglingNodes().erase(n0.second->sID); // remove temporary shared_ptr
//            danglingNodes().erase(n1.second->sID); // remove temporary shared_ptr

        }
        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        std::pair<typename LoopNodeContainerType::iterator,bool> expand(const size_t& a, const size_t& b, const NodeArgTypes&... Args)
        {
            
            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
            
            const size_t i=std::min(a,b);
            const size_t j=std::max(a,b);
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Expanding non-existing link");
            
            // Insert new node
            std::pair<typename LoopNodeContainerType::iterator,bool> newNode=insertDanglingNode(Args...);
            const size_t newID=newNode.first->first;
            
            typename LinkType::LoopLinkContainerType loopLinksIJ=Lij.second->loopLinks();
            for(const auto& llink : loopLinksIJ)
            {
                if(llink->source->sID==i && llink->sink->sID==j) // same direction
                {
                    loopLinks().emplace(std::piecewise_construct,
                                        std::make_tuple(i,newID),
                                        std::make_tuple(llink->source,newNode.first->second, llink->pLoop) );

                    loopLinks().emplace(std::piecewise_construct,
                                        std::make_tuple(newID,j),
                                        std::make_tuple(newNode.first->second,llink->sink,llink->pLoop) );
                }
                else if(llink->source->sID==j && llink->sink->sID==i) // opposite direction
                {
                    loopLinks().emplace(std::piecewise_construct,
                                        std::make_tuple(j,newID),
                                        std::make_tuple(llink->source,newNode.first->second, llink->pLoop) );

                    loopLinks().emplace(std::piecewise_construct,
                                        std::make_tuple(newID,i),
                                        std::make_tuple(newNode.first->second,llink->sink,llink->pLoop) );
                }
                else
                {
                    assert(0 && "LoopLink must connect i->j or viceversa");
                }
            }
            
            danglingNodes().erase(newID);
            loopLinks().erase(std::make_pair(i,j));
            loopLinks().erase(std::make_pair(j,i));
            
            return newNode;
        }
        
        /**********************************************************************/
        void contractSecond(const size_t& a, const size_t& b)
        {
            
            const size_t i=std::min(a,b);
            const size_t j=std::max(a,b);
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Contracting non-existing link");
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
        void checkLoops()
        {
            std::cout<<"Checking "<<LoopObserverType::loops().size()<<" loops"<<std::endl;
            for(const auto& loop : LoopObserverType::loops())
            {
                assert(loop.second->isLoop() && "note a closed loop");
            }
        
        }
        
    };
    
    
}
#endif
