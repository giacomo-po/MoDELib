/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNetwork_H_
#define model_LoopNetwork_H_


template <typename T>
class LoopNode;

template <typename T>
class LoopLink;

template <typename T>
class NetworkLink;


#include <iostream>
#include <assert.h>
#include <memory>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <utility>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <CRTP.h>
#include <MPIcout.h>
//#include <NetworkNode.h>
//#include <LoopLink.h>
//#include <NetworkLink.h>
////#include <NetworkLinkObserver.h>
#include <Loop.h>
#include <LoopNode.h>
#include <LoopLink.h>

#include <NetworkNode.h>
#include <WeakPtrFactories.h>

////#include <LoopObserver.h>
//#include <NetworkComponent.h>
//#include <NetworkComponentObserver.h>
#include <EqualIteratorRange.h>
//#include <KeyConstructableSharedPtrFactory.h>

#define VerboseLoopNetwork(N,x) if(verboseLevel>=N){model::cout<<x;}

namespace model
{
    
    template<typename FlowType>
    struct NullFlow
    {
        static bool isZero(const FlowType&)
        {
            return false;
        }
    };
    
    
    

    
    
    template<typename Derived>
    class LoopNetwork : public CRTP<Derived>
    //    /*               */ public NetworkComponentObserver<NetworkComponent<typename TypeTraits<Derived>::NodeType,typename TypeTraits<Derived>::LinkType>>,
    /*               */,private std::map<typename TypeTraits<Derived>::LoopLinkType::KeyType,typename TypeTraits<Derived>::LoopLinkType>
    /*               */,public WeakPtrFactory<Derived,typename TypeTraits<Derived>::LoopType>
    /*               */,public WeakPtrFactory<Derived,typename TypeTraits<Derived>::LoopNodeType>
    /*               */,public WeakPtrFactory<Derived,typename TypeTraits<Derived>::NetworkNodeType>
    /*               */,public KeyConstructableWeakPtrFactory<Derived,typename TypeTraits<Derived>::NetworkLinkType>
    {
        
    public:
        
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<Derived>::NetworkNodeType NetworkNodeType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<Derived>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::FlowType FlowType;
        typedef WeakPtrFactory<Derived,LoopType> LoopContainerType;
        typedef std::map<typename LoopLinkType::KeyType,LoopLinkType> LoopLinkContainerType;
        typedef WeakPtrFactory<Derived,LoopNodeType> LoopNodeContainerType;
        typedef WeakPtrFactory<Derived,NetworkNodeType> NetworkNodeContainerType;
        typedef KeyConstructableWeakPtrFactory<Derived,NetworkLinkType> NetworkLinkContainerType;
        typedef typename LoopNodeContainerType::SharedPtrType SharedLoopNodePtrType;
        typedef typename NetworkNodeContainerType::SharedPtrType SharedNetworkNodePtrType;

        const LoopContainerType& loops() const
        {
            return *this;
        }
        
        LoopContainerType& loops()
        {
            return *this;
        }
        
        const LoopNodeContainerType& loopNodes() const
        {
            return *this;
        }
        
        LoopNodeContainerType& loopNodes()
        {
            return *this;
        }
        
        const NetworkNodeContainerType& networkNodes() const
        {
            return *this;
        }
        
        NetworkNodeContainerType& networkNodes()
        {
            return *this;
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
        NetworkLinkContainerType& networkLinks()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        const NetworkLinkContainerType& networkLinks() const
        {
            return *this;
        }
        
    private:
        
        std::pair<SharedLoopNodePtrType,SharedLoopNodePtrType> cutLoop(const SharedLoopNodePtrType& Ni,const SharedLoopNodePtrType& Nj)
        {
            const auto loop(Ni->loop());
            VerboseLoopNetwork(1,"cutting loop "<<loop->tag()<<" at "<<Ni->tag()<<","<<Nj->tag()<<std::endl);
            if(Ni!=Nj)
            {
                const auto linkAtI=loop->linkStartingAt(Ni);
                const auto linkAtJ=loop->linkStartingAt(Nj);
                
                if(linkAtI.first && linkAtJ.first)
                {
                    if(linkAtI.second->sink!=Nj && linkAtJ.second->sink!=Ni)
                    {
                        auto currLink(linkAtI.second);
                        if(currLink->sink!=Nj)
                        {
                            std::shared_ptr<LoopType> newLoop(loops().clone(loop->key));
                            
                            const auto cloneNi(loopNodes().clone(Ni->key,newLoop,networkNodes().clone(Ni->networkNode->key)));
                            const auto cloneNj(loopNodes().clone(Nj->key,newLoop,networkNodes().clone(Nj->networkNode->key)));
                            
                            std::vector<std::tuple<std::shared_ptr<LoopNodeType>,std::shared_ptr<LoopNodeType>>> linksToDisconect;
                            while(currLink->source!=Nj)
                            {
                                linksToDisconect.emplace_back(currLink->source,currLink->sink);
                                currLink=currLink->next;
                            }
                            
                            
                            for(const auto& tup : linksToDisconect)
                            {
                                disconnect(std::get<0>(tup),std::get<1>(tup),loop);
                                //                                        connect(std::get<0>(tup),std::get<1>(tup),newLoop);
                            }
                            
                            for(const auto& tup : linksToDisconect)
                            {
                                const std::shared_ptr<LoopNodeType> source(std::get<0>(tup)==Ni? cloneNi : std::get<0>(tup));
                                const std::shared_ptr<LoopNodeType> sink(std::get<1>(tup)==Nj? cloneNj : std::get<1>(tup));
                                source->resetLoop(newLoop);
                                sink->resetLoop(newLoop);
                                connect(source,sink,newLoop);
                            }
                            
                            //                                    connect(std::get<1>(linksToDisconect.back()),std::get<0>(linksToDisconect.front()),newLoop); // cloneNj->cloneNi
                            connect(cloneNj,cloneNi,newLoop); // cloneNj->cloneNi
                            connect(std::get<0>(linksToDisconect.front()),std::get<1>(linksToDisconect.back()),loop); // Ni->Nj
                            return std::make_pair(cloneNi,cloneNj);
                        }
                    }
                }
            }
            return std::make_pair(SharedLoopNodePtrType(nullptr),SharedLoopNodePtrType(nullptr));
        }
        
        /**********************************************************************/
        void connect(const SharedLoopNodePtrType& n0,
                     const SharedLoopNodePtrType& n1,
                     const std::shared_ptr<LoopType>& tempLoop)
        {
            VerboseLoopNetwork(1,"connecting "<<n0->sID<<"->"<<n1->sID<<" ("<<tempLoop->sID<<")"<<std::endl);
            
            if(n0->sID!=n1->sID)
            {
                const auto key=LoopLinkType::loopLinkKey(tempLoop,n0,n1); // even for LoopLinks we store with key (minID,maxID)
                const auto linkIter(loopLinks().find(key));
                if(linkIter==loopLinks().end())
                {// link not found
                    loopLinks().try_emplace(key,n0,n1, tempLoop);
                }
                else
                {// link found
                    if(linkIter->second.source->sID==n1->sID && linkIter->second.sink->sID==n0->sID)
                    {// If opposite direction disconnect
                        loopLinks().erase(linkIter);
                    }
                }
            }
        }

        /**********************************************************************/
        size_t disconnect(const SharedLoopNodePtrType& n0, const SharedLoopNodePtrType& n1, const std::shared_ptr<LoopType>& loop)
        {
            VerboseLoopNetwork(1,"disconnecting "<<n0->sID<<"->"<<n1->sID<<" ("<<loop->sID<<")"<<std::endl);
            return loopLinks().erase(LoopLinkType::loopLinkKey(loop,n0,n1));
        }

    public:
        //
        static int verboseLevel;

        /**********************************************************************/
        void insertLoop(const std::shared_ptr<LoopType>& loop,
                        const std::vector<std::shared_ptr<LoopNodeType>> nodes)
        {/*!@param[in] nodes in the loop
          * @param[loopInput] Loop constructor arguments
          *
          * Inserts a Loop connecting the sequence of nodes,
          * The loop constructor arguments loopInput
          * are forwarded to the loop constructor.
          */
            
            for(size_t k=0;k<nodes.size();++k)
            {
                const size_t next=k+1<nodes.size()? k+1 : 0;
                connect(nodes[k],nodes[next],loop);
            }
            
            assert(loop->isLoop() && "Not a loop.");
            
            loop->update();
        }
        
        /**********************************************************************/
        void deleteLoop(const size_t& loopID)
        {
            for(typename LoopLinkContainerType::const_iterator loopIter=loopLinks().begin();
                /*                                          */ loopIter!=loopLinks().end();)
            {
                if(loopIter->second.loop()->sID==loopID)
                {
                    loopIter=loopLinks().erase(loopIter);
                }
                else
                {
                    loopIter++;
                }
            }
        }
        
        
        SharedNetworkNodePtrType expandNetworkLink(const std::shared_ptr<NetworkLinkType>& Lij, const SharedNetworkNodePtrType& networkNode)
        {
            VerboseLoopNetwork(1,"expanding NetworkLink "<<Lij->tag()<<" at node "<<networkNode->tag()<<std::endl);
            const auto loopLinks(Lij->loopLinks()); // copy loop links
            for(const auto& loopLink : loopLinks)
            {
                expandLoopLink(*loopLink,loopNodes().create(loopLink->loop,networkNode));
            }
            return networkNode;
        }
        
        /**********************************************************************/
        SharedLoopNodePtrType expandLoopLink(const LoopLinkType& Lij, const SharedLoopNodePtrType& loopNode)
        {
            VerboseLoopNetwork(1,"expanding LoopLink "<<Lij.tag()<<" at node "<<loopNode->tag()<<std::endl);
            const SharedLoopNodePtrType source(Lij.source);
            const SharedLoopNodePtrType sink(Lij.sink);
            const std::shared_ptr<LoopType> loop(Lij.loop);
            disconnect(source,sink,loop);
            connect(source,loopNode,loop);
            connect(loopNode,sink,loop);
            return loopNode;
        }

        /**********************************************************************/
        void removeLoopNode(const SharedLoopNodePtrType& node)
        {
            if(node)
            {
                const auto loop(node->loop());
                if(loop)
                {
                    const SharedLoopNodePtrType prev(loopNodes().get(node->prev.first->key));
                    const SharedLoopNodePtrType next(loopNodes().get(node->next.first->key));
                    if(prev)
                    {
                        disconnect(prev,node,loop);
                    }
                    if(next)
                    {
                        disconnect(node,next,loop);
                    }
                    if(prev && next)
                    {
                        connect(prev,next,loop);
                    }
                }
            }
        }

        /**********************************************************************/
        bool contractLoopNodes(const LoopNodeType* const Ni,const LoopNodeType* const Nj)
        {
            if(Ni!=nullptr && Nj!=nullptr)
            {
                const auto sharedNi(loopNodes().get(Ni->key));
                const auto sharedNj(loopNodes().get(Nj->key));
                assert(sharedNi && sharedNj);
                return contractLoopNodes(sharedNi,sharedNj);
            }
            return false;
        }
        
        /**********************************************************************/
        bool contractLoopNodes(const SharedLoopNodePtrType& Ni,const SharedLoopNodePtrType& Nj)
        {
            if(Nj->isContractableTo(Ni.get()))
            {
                const auto clones(cutLoop(Ni,Nj));
                removeLoopNode(Nj);
                if(clones.first && clones.second)
                {
                    removeLoopNode(clones.second);
                }
                return true;
            }
            return false;
        }
        
        /**********************************************************************/
        bool contractNetworkNodes(const SharedNetworkNodePtrType& Ni,const SharedNetworkNodePtrType& Nj)
        {
            if(Nj->isContractableTo(Ni))
            {
                std::deque<std::pair<const LoopNodeType* const,const LoopNodeType* const>> loopNodesToContract;
                std::deque<std::pair<const LoopNodeType* const,      LoopNodeType* const>> loopNodesToReassign;

                for(const auto& loopNodeI : Ni->loopNodes())
                {
                        for(const auto& loopNodeJ : Nj->loopNodes())
                        {
                            if(loopNodeJ->isContractableTo(loopNodeI))
                            {// this implies loopNodeI and loopNodeJ in same loop
                                loopNodesToContract.emplace_back(loopNodeI,loopNodeJ);
                            }
                            else
                            {// two cases: 1) loopNodeI and loopNodeJ in different loops 2) loopNodeI and loopNodeJ in same loop but notContractable
                                loopNodesToReassign.emplace_back(loopNodeI,loopNodeJ);
                            }
                        }
                }
                
                for(const auto pair : loopNodesToContract)
                {
                    contractLoopNodes(pair.first,pair.second);
                }
                
                for(auto& pair : loopNodesToReassign)
                {
                    pair.second->networkNode->removeLoopNode(pair.second);
                    pair.second->networkNode=pair.first->networkNode;
                    pair.second->networkNode->addLoopNode(pair.second);
                    if(pair.second->prev.second)
                    {
                        pair.second->prev.second->resetNetworkLink();
                    }
                    if(pair.second->next.second)
                    {
                        pair.second->next.second->resetNetworkLink();
                    }
                }
            }
            
            return false;
        }
        
        /**********************************************************************/
        void printLoops() const
        {
            for(const auto& loop : loops())
            {
                if(!loop.second.expired())
                {
                    loop.second.lock()->printLoop();
                }
            }
        }
        
        /**********************************************************************/
        void printLoopNodes() const
        {
            std::cout<<"Printing LoopNodes "<<loopNodes().size()<<std::endl;
            for(const auto& node : loopNodes())
            {
                if(!node.second.expired())
                {
                    std::cout<<node.second.lock()->tag()<<std::endl;
                }
            }
        }
        
        /**********************************************************************/
        void printNetworkNodes() const
        {
            std::cout<<"Printing NetworkNodes "<<networkNodes().size()<<std::endl;
            for(const auto& node : networkNodes())
            {
                if(!node.second.expired())
                {
                    std::cout<<node.second.lock()->tag()<<std::endl;
                }
            }
        }

    };
    
    template<typename Derived>
    int LoopNetwork<Derived>::verboseLevel=0;
    
}
#endif
