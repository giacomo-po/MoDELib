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

#include <model/Utilities/CRTP.h>
#include <model/MPI/MPIcout.h>
#include <model/LoopNetwork/LoopNode.h>
//#include <model/LoopNetwork/NetworkNode.h>
#include <model/LoopNetwork/LoopLink.h>
#include <model/LoopNetwork/NetworkLink.h>
#include <model/LoopNetwork/NetworkLinkObserver.h>
#include <model/LoopNetwork/Loop.h>
#include <model/LoopNetwork/LoopObserver.h>

#define VerboseLoopNetwork(N,x) if(verboseLevel>=N){model::cout<<x;}

namespace model
{
    template<typename Derived>
    class LoopNetwork : public CRTP<Derived>,
    //    /*               */ private std::map<size_t,typename TypeTraits<Derived>::NodeType>,
    /*               */ private std::map<size_t,std::shared_ptr<typename TypeTraits<Derived>::NodeType>>,
    /*               */ private std::multimap<std::pair<size_t,size_t>,LoopLink<typename TypeTraits<Derived>::LinkType>>,
    /*               */ public NetworkLinkObserver<typename TypeTraits<Derived>::LinkType>,
    /*               */ public LoopObserver<typename TypeTraits<Derived>::LoopType>,
    /*               */ public NodeObserver<typename TypeTraits<Derived>::NodeType>,
    /*               */ private std::map<size_t,const typename TypeTraits<Derived>::LoopType* const>
    {
        
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        //        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<typename TypeTraits<Derived>::LinkType> LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::FlowType FlowType;
        
        //        typedef std::map<size_t,NodeType>     LoopNodeContainerType;
        
        typedef std::multimap<std::pair<size_t,size_t>,LoopLinkType> LoopLinkContainerType;
        //        typedef std::map<size_t,const LoopType* const> LoopContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        typedef typename NetworkLinkObserverType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename NetworkLinkObserverType::IsConstNetworkLinkType IsConstNetworkLinkType;
        typedef LoopObserver<LoopType> LoopObserverType;
        
        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef std::map<size_t,SharedNodePtrType>     LoopNodeContainerType;
        
        typedef std::pair<bool,SharedNodePtrType>          IsNodeType;
        typedef std::pair<bool,const std::shared_ptr<const NodeType>>	 IsConstNodeType;
        
        typedef std::pair<bool,LoopLinkType*> IsLoopLinkType;
        
        typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>> ExpandTupleType;
        
        typedef std::tuple<SharedNodePtrType, // A node, preceding B in loop being merged
        /*              */ SharedNodePtrType, // B node in loop being merged, which will become node E
        /*              */ SharedNodePtrType, // C node in loop being merged, which will become node F
        /*              */ SharedNodePtrType, // D node, after C in loop being merged
        /*              */ std::shared_ptr<LoopType>, // loop being merged
        /*              */ SharedNodePtrType, // E node, will replace B
        /*              */ SharedNodePtrType, // F node, will replace C
        /*              */ std::shared_ptr<LoopType> // loop being merged
        > AMTtupleType;
        
        
        //        typedef std::map<int,AMTtupleType> MergeMapType; // map is sorted by operation priority
        typedef std::deque<AMTtupleType> AMTdequeType;
        
        /**********************************************************************/
        LoopNodeContainerType& danglingNodes()
        {
            return *this;
        }
        
        const LoopNodeContainerType& danglingNodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        IsNodeType danglingNode(const size_t & k)
        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
            typename LoopNodeContainerType::iterator nodeIter(danglingNodes().find(k));
            return (nodeIter!=danglingNodes().end())? std::make_pair(true,nodeIter->second) : std::make_pair(false,SharedNodePtrType(nullptr));
        }
        
        
        
        /**********************************************************************/
        void connect(const SharedNodePtrType& n0,
                     const SharedNodePtrType& n1,
                     const std::shared_ptr<LoopType>& tempLoop)
        {
            VerboseLoopNetwork(1,"connecting "<<n0->sID<<"->"<<n1->sID<<std::endl);
//            assert(n0->sID!=n1->sID && "Cannot connect a node to itself");
            
            if(n0->sID!=n1->sID)
            {
                //check if n1->n0 exist for the same loop
                //            const auto iterPair = loopLinks().equal_range(std::pair<size_t,size_t>(n1->sID,n0->sID));
                const auto key=LoopLinkType::getKey(n0,n1);
                const auto iterPair = loopLinks().equal_range(key);
                
                typename LoopLinkContainerType::const_iterator loopIter=iterPair.second;
                for (loopIter=iterPair.first;loopIter!=iterPair.second;++loopIter)
                {
                    if(loopIter->second.loop().get()==tempLoop.get() && // link in same loop
                       loopIter->second.source()->sID==n1->sID && loopIter->second.sink()->sID==n0->sID) // opposite link exists
                    {
                        break;
                    }
                    
                }
                
                if(loopIter!=iterPair.second)
                {// opposite link exists, disconnect it
                    loopLinks().erase(loopIter);
                }
                else
                {// opposite link does not exists, connect n0->n1
                    loopLinks().emplace(std::piecewise_construct,
                                        //                                    std::make_tuple(n0->sID,n1->sID),
                                        std::make_tuple(key.first,key.second),
                                        std::make_tuple(n0,n1, tempLoop) );
                    
                }
            }

        }
        
        /**********************************************************************/
        IsLoopLinkType loopLink(const SharedNodePtrType& n0, const SharedNodePtrType& n1, const std::shared_ptr<LoopType>& loop)
        {
            
            IsLoopLinkType temp=IsLoopLinkType(false,nullptr);
            
            const auto key=LoopLinkType::getKey(n0,n1);
            const auto iterPair = loopLinks().equal_range(key);
            for(typename LoopLinkContainerType::iterator loopIter=iterPair.first;
                /*                                          */ loopIter!=iterPair.second;
                /*                                          */ loopIter++)
            {
                if(loopIter->second.loop().get()==loop.get() &&
                   loopIter->second.source()->sID==n0->sID && loopIter->second.sink()->sID==n1->sID)
                {
                    temp=IsLoopLinkType(true,&loopIter->second);
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        size_t disconnect(const SharedNodePtrType& n0, const SharedNodePtrType& n1, const std::shared_ptr<LoopType>& loop)
        {
            VerboseLoopNetwork(1,"disconnecting "<<n0->sID<<"->"<<n1->sID<<std::endl);
            size_t nDisconnected=0;
            const auto key=LoopLinkType::getKey(n0,n1);
            const auto iterPair = loopLinks().equal_range(key);
            
            for(typename LoopLinkContainerType::const_iterator loopIter=iterPair.first;
                /*                                          */ loopIter!=iterPair.second;
                /*                                          */ loopIter++)
            {
                if(loopIter->second.loop().get()==loop.get() &&
                   loopIter->second.source()->sID==n0->sID && loopIter->second.sink()->sID==n1->sID)
                {
                    loopLinks().erase(loopIter);
                    nDisconnected=1;
                    nDisconnected+=disconnect(n0,n1,loop); // make sure that recursive calls don't find another link
                    assert(nDisconnected==1 && "More than one LoopLink with same key and same loop exist.");
                    break;
                }
            }
            return nDisconnected;
        }
        
        
    public:
        
        static int verboseLevel;
        
        
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
        
        const NetworkLinkContainerType networkLinks() const
        {
            return NetworkLinkObserverType::links();
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
        template <typename ...NodeArgTypes>
        std::pair<typename LoopNodeContainerType::iterator,bool> insertDanglingNode(const NodeArgTypes&... nodeInput)
        {/*! @param[in] nodeInput
          *\returns
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
        void insertLoop(const std::vector<size_t> nodeIDs,
                        const FlowType& f,
                        const LoopArgTypes&... loopInput)
        {/*!@param[in] nodeIDs IDs of the nodes in the loop
          * @param[in] f the loop flow
          * @param[loopInput] additional Loop constructor arguments
          *
          * Inserts a Loop connecting the sequence of nodes with IDs nodeIDs,
          * flow f. The additional loop constructor arguments loopInput
          * are forwarded to the loop constructor.
          */
            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->derived(),f,loopInput...);
            
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
        void removeLoop(const std::shared_ptr<LoopType>& pL)
        {/*!\param[in] pL a shared_ptr to the loop to be removed
          *
          * Disconnects all segments in loop pL, therefore removing the loop itself.
          */
            for(const auto& lLink : pL->linkSequence())
            {
                disconnect(pL->source,pL->sink,pL);
            }
        }
        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        SharedNodePtrType expand(const size_t& a, const size_t& b, const NodeArgTypes&... Args)
        {
            VerboseLoopNetwork(1,"expanding "<<a<<","<<b<<std::endl);
            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
            
            const auto key=LoopLinkType::getKey(a,b);
            const size_t& i=key.first;
            const size_t& j=key.second;
            //            const size_t i=std::min(a,b);
            //            const size_t j=std::max(a,b);
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Expanding non-existing link");
            
            
            //            std::cout<<;
            
            // Create new node
            SharedNodePtrType newNode=SharedNodePtrType(new NodeType(Args...));
            //            const size_t newID=newNode->sID;
            
            // Store what needs to be connected (i,new,j,Loop), or (j,new,i,Loop). This also holds temporarily disconnected nodes
            std::deque<ExpandTupleType> expandDeq;
            for(const auto& llink : Lij.second->loopLinks())
            {
                expandDeq.emplace_back(llink->source,llink->sink,llink->loop());
            }
            
            // Delete all LoopLinks of type i->j or j->i
            //            const auto key=LoopLinkType::getKey(a,b);
            loopLinks().erase(key);
            
            //            loopLinks().erase(std::make_pair(i,j));
            //            loopLinks().erase(std::make_pair(j,i));
            
            for (const auto& tup : expandDeq)
            {
                //                std::cout<<std::get<0>(tup)->sID<<"->"<<std::get<1>(tup)->sID<<", loop "<<std::get<2>(tup)->sID<<std::endl;
                
                connect(std::get<0>(tup),newNode, std::get<2>(tup));
                connect(newNode,std::get<1>(tup), std::get<2>(tup));
            }
            
            return newNode;
        }
        
//        /**********************************************************************/
//        void contract(const size_t& a, const size_t& b)
//        {/*! Performs a contract operation in which the nodes of a NetworkLink
//          * (a,b) are merged. Node a survives the merge, while node b is
//          * destroyed
//          */
//            
//            THIS FUNCTION IS WRONG. IF THERE IS A LOOP ATTACHED TO EITHER I OR J (NOT BOTH), THEN THAT NODE SURVIVES!
//            
//            VerboseLoopNetwork(1,"contracting "<<a<<","<<b<<std::endl);
//            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
//            
//            const auto key=LoopLinkType::getKey(a,b);
//            const size_t& i=key.first;
//            const size_t& j=key.second;
//            //            const size_t i=std::min(a,b);
//            //            const size_t j=std::max(a,b);
//            
//            
//            // Find NetworkLink i->j
//            const IsConstNetworkLinkType Lij(this->link(i,j));
//            assert(Lij.first && "Contracting non-existing link");
//            
//            typedef std::pair<size_t,size_t> DisconnectPairType;
//            std::deque<DisconnectPairType> disconnectDeq;
//            
//            
//            typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>> ReconnectTupleType;
//            std::deque<ReconnectTupleType> reconnectDeq;
//            
//            for(const auto& loopLink : Lij.second->loopLinks())
//            {
//                if(loopLink->source()->sID==b)
//                {
//                    assert(loopLink->source()->sID==b && loopLink->sink()->sID==a);
//                    disconnectDeq.emplace_back(loopLink->prev->source()->sID,b);
//                    disconnectDeq.emplace_back(b,a);
//                    reconnectDeq.emplace_back(loopLink->prev->source(),loopLink->sink(),loopLink->loop());
//                }
//                else
//                {
//                    assert(loopLink->source()->sID==a && loopLink->sink()->sID==b);
//                    disconnectDeq.emplace_back(a,b);
//                    disconnectDeq.emplace_back(b,loopLink->next->sink()->sID);
//                    reconnectDeq.emplace_back(loopLink->source(),loopLink->next->sink(),loopLink->loop());
//                }
//                
//            }
//            
//            for (const auto& pair : disconnectDeq)
//            {
//                loopLinks().erase(pair);
//            }
//            
//            for (const auto& tup : reconnectDeq)
//            {
//                connect(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup));
//            }
//            
//        }
        
        
        /**********************************************************************/
        void cutLoop(const size_t& L,
                     const size_t& i,const size_t &j)
        {
            const auto loop=this->loop(L);
            if(loop.first && i!=j)
            {
//                const auto nodeSeq=loop.seocond->nodeIDSequence();
                const auto linkAtI=loop.second->linkStartingAt(i);
                const auto linkAtJ=loop.second->linkStartingAt(j);
                
                if(linkAtI.first && linkAtJ.first)
                {
                    if(linkAtI.second->sink()->sID!=j && linkAtJ.second->sink()->sID!=i)
                    {
                        VerboseLoopNetwork(1,"cutting loop "<<L<<" at "<<i<<" and "<<j<<std::endl);

                        std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->derived(),linkAtI.second->flow());
                        linkAtI.second->resetLoop(tempLoop,i,j);
                        
                        linkAtI.second->prev->next=nullptr;
                        linkAtI.second->prev=nullptr;

                        //
                        linkAtJ.second->prev->next=nullptr;
                        linkAtJ.second->prev=nullptr;

                        
                        connect(linkAtJ.second->source(),linkAtI.second->source(),tempLoop);
                        connect(linkAtI.second->source(),linkAtJ.second->source(),linkAtJ.second->loop());

                    }
                    
                }
//                const auto linkSeq=loop.second->linkSequence();
            }
        }
        
        /**********************************************************************/
        bool contractSecond(const size_t& a,const size_t &b)
        {
            VerboseLoopNetwork(1,"contracting "<<a<<","<<b<<std::endl);
            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");

            
            bool success=false;
            const auto nA=this->sharedNode(a);
            const auto nB=this->sharedNode(b);
            
            if(nA.first && nB.first)
            {
                

                
                // Collect IDs of all loops passing through both a and b
                std::set<size_t> loopIDs;
                for(const auto& loopLinkA : nA.second->loopLinks())
                {
                    for(const auto& loopLinkB : nB.second->loopLinks())
                    {
                        
                        if(loopLinkA->loop()->sID==loopLinkB->loop()->sID)
                        {
                            loopIDs.insert(loopLinkA->loop()->sID);
                        }
                    }
                }
                
                // Cut those loops
                for(const auto loopID : loopIDs)
                {
                    cutLoop(loopID,a,b);
                }
                
                
                // Store links connected to b
                typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>,bool> ContractTupleType;
                typedef std::deque<ContractTupleType> ContractDequeType;
                ContractDequeType contractDeq;
                
                for(const auto& loopLink : nB.second->loopLinks())
                {
                
                    if(loopLink->source()->sID==b)
                    {
                        contractDeq.emplace_back(loopLink->source(),loopLink->sink(),loopLink->loop(),0);
                    }
                    else if(loopLink->sink()->sID==b)
                    {
                        contractDeq.emplace_back(loopLink->source(),loopLink->sink(),loopLink->loop(),1);
                    }
                    else
                    {
                        assert(0 && "source or sink must be b.");
                    }

                    
                }
                
                // Disconnect links from b
                for(const auto& tup : contractDeq)
                {
                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup));
                }

                // Re-connect replacing b with a
                for(const auto& tup : contractDeq)
                {
                    if(std::get<3>(tup))
                    {
                        connect(std::get<0>(tup),nA.second,std::get<2>(tup));
                    }
                    else
                    {
                        connect(nA.second,std::get<1>(tup),std::get<2>(tup));
                    }
                }

                
                
            }
            
            return success;
        }
        
        
        
        
//        /**********************************************************************/
//        void annihilateTransfer(const size_t& a, const size_t& b,const size_t& a1, const size_t& b1)
//        {
//            annihilate(a,b,a1,b1);
//            transfer(a,b,a1,b1);
//        }
//        
//        /**********************************************************************/
//        void annihilate(const size_t& a, const size_t& b,const size_t& a1, const size_t& b1)
//        {
//            VerboseLoopNetwork(1,"annihilating ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
//            
//            if(a!=a1 || b!=b1) // avoid trivia merge
//            {
//                // Find NetworkLink i->j
//                const auto key=LoopLinkType::getKey(a,b);
//                const size_t& i=key.first;
//                const size_t& j=key.second;
//                const IsConstNetworkLinkType Lij(this->link(i,j));
//                
//                // Find NetworkLink i1->j1
//                const auto key1=LoopLinkType::getKey(a1,b1);
//                const size_t& i1=key1.first;
//                const size_t& j1=key1.second;
//                const IsConstNetworkLinkType Li1j1(this->link(i1,j1));
//                
//                if(Lij.first && Li1j1.first)
//                {
//                    const SharedNodePtrType nA((Lij.second->source->sID==a)? Lij.second->source : Lij.second->sink);
//                    const SharedNodePtrType nB((Lij.second->source->sID==b)? Lij.second->source : Lij.second->sink);
//                    
//                    AMTdequeType deq;
//                    
//                    // For each loop link in Li1j1, store merge operations in a MergeMapType
//                    for (const auto& loopLink1 : Li1j1.second->loopLinks())
//                    {
//                        
//                        size_t nAnnihilations=0;
//                        
//                        for (const auto& loopLink : Lij.second->loopLinks())
//                        {
//                            if(loopLink1->loop().get()==loopLink->loop().get()) // annihilation
//                            {
//                                
//                                if(loopLink->source()->sID==b && loopLink->sink()->sID==a &&
//                                   loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1)
//                                {
//                                    
//                                    deq.emplace_back(loopLink1->prev->source(),
//                                                     loopLink1->source(),
//                                                     loopLink1->sink(),
//                                                     loopLink1->next->sink(),
//                                                     loopLink1->loop(),
//                                                     nA,
//                                                     nB,
//                                                     loopLink->loop() // unused
//                                                     );
//                                    
//                                }
//                                else if(loopLink->source()->sID==a && loopLink->sink()->sID==b &&
//                                        loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1)
//                                {
//                                    
//                                    deq.emplace_back(loopLink1->prev->source(),
//                                                     loopLink1->source(),
//                                                     loopLink1->sink(),
//                                                     loopLink1->next->sink(),
//                                                     loopLink1->loop(),
//                                                     nB,
//                                                     nA,
//                                                     loopLink->loop() // unused
//                                                     );
//                                }
//                                else
//                                {
//                                    assert(0 && "Invalid annihilation, same link direction.");
//                                }
//                                nAnnihilations++;
//                            }
//                            
//                        }
//                        
//                        assert(nAnnihilations<=1);
//                    }
//                    
//                    // Perform annihilate
//                    for(const auto tup : deq)
//                    {
//                        if(std::get<0>(tup)->sID==std::get<6>(tup)->sID || std::get<0>(tup)->sID==std::get<5>(tup)->sID)
//                        {// part of the loop is pinched off, no need to create another loop
//                            disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                            disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                            disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                            connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                            disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                        }
//                        else if(std::get<3>(tup)->sID==std::get<6>(tup)->sID || std::get<3>(tup)->sID==std::get<5>(tup)->sID)
//                        {// part of the loop is pinched off, no need to create another loop
//                            disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                            disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                            disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                            connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                            disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                        }
//                        else
//                        {
//                            // Break the loop in two parts
//                            disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                            disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                            
//                            // links 0->1 and 2->3 should still exist, but now they are disconnected, so reset the loop ptr on on of them
//                            IsLoopLinkType link01=loopLink(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                            assert(link01.first && "LoopLink must exist");
//                            
//                            IsLoopLinkType link23=loopLink(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                            assert(link23.first && "LoopLink must exist");
//                            
//                            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->derived(),std::get<4>(tup)->flow());
//                            link01.second->resetLoop(tempLoop);
//                            
//                            // Finish disconnecting
//                            disconnect(std::get<0>(tup),std::get<1>(tup),tempLoop); // after resetLoop link01 is in loop tempLoop
//                            disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                            
//                            connect(std::get<0>(tup),std::get<5>(tup),tempLoop);
//                            connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                            
//                        }
//                    }
//                    
//                }
//                
//            }
//            
//            //                assert(Li1j1.first && "Merging non-existing link");
//            
//            
//            VerboseLoopNetwork(1,"done annhilating ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            
//        }
        
//        /**********************************************************************/
//        void transfer(const size_t& a, const size_t& b,const size_t& a1, const size_t& b1)
//        {
//            VerboseLoopNetwork(1,"transfering ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
//            
//            if(a!=a1 || b!=b1) // avoid trivia merge
//            {
//                // Find NetworkLink i->j
//                const auto key=LoopLinkType::getKey(a,b);
//                const size_t& i=key.first;
//                const size_t& j=key.second;
//                const IsConstNetworkLinkType Lij(this->link(i,j));
//                
//                // Find NetworkLink i1->j1
//                const auto key1=LoopLinkType::getKey(a1,b1);
//                const size_t& i1=key1.first;
//                const size_t& j1=key1.second;
//                const IsConstNetworkLinkType Li1j1(this->link(i1,j1));
//                
//                
//                if(Lij.first && Li1j1.first)
//                {
//                    std::cout<<"I'm here!"<<std::endl;
//                    const SharedNodePtrType nA((Lij.second->source->sID==a)? Lij.second->source : Lij.second->sink);
//                    const SharedNodePtrType nB((Lij.second->source->sID==b)? Lij.second->source : Lij.second->sink);
//                    
//                    AMTdequeType deq;
//                    
//                    // For each loop link in Li1j1, store merge operations in a MergeMapType
//                    for (const auto& loopLink1 : Li1j1.second->loopLinks())
//                    {
//                        
//                        if(loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1)
//                        {
//                            deq.emplace_back(loopLink1->prev->source(),
//                                             loopLink1->source(),
//                                             loopLink1->sink(),
//                                             loopLink1->next->sink(),
//                                             loopLink1->loop(),
//                                             nA,
//                                             nB,
//                                             loopLink1->loop() // unused
//                                             );
//                            
//                        }
//                        else if(loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1)
//                        {
//                            deq.emplace_back(loopLink1->prev->source(),
//                                             loopLink1->source(),
//                                             loopLink1->sink(),
//                                             loopLink1->next->sink(),
//                                             loopLink1->loop(),
//                                             nB,
//                                             nA,
//                                             loopLink1->loop() // unused
//                                             );
//                        }
//                        else
//                        {
//                            assert(0);
//                        }
//                        
//                    }
//                    
//                    // Perform merging
//                    for(const auto tup : deq)
//                    {
//                        disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                        connect(std::get<5>(tup),std::get<6>(tup),std::get<4>(tup));
//                        
//                        if(std::get<1>(tup)->sID!=std::get<5>(tup)->sID)
//                        {
//                            disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                            connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                        }
//                        
//                        if(std::get<2>(tup)->sID!=std::get<6>(tup)->sID)
//                        {
//                            disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                            connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                        }
//                        
//                    }
//                    
//                }
//                
//            }
//            
//            VerboseLoopNetwork(1,"done transfer ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            
//        }
        
//        /**********************************************************************/
//        void mergeOLD(const size_t& a, const size_t& b,const size_t& a1, const size_t& b1)
//        {
//            VerboseLoopNetwork(1,"merging ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            assert(danglingNodes().empty() && "You must call clearDanglingNodes() after inserting all loops.");
//            
//            if(a!=a1 || b!=b1) // avoid trivia merge
//            {
//                // Find NetworkLink i->j
//                //                const size_t i=std::min(a,b);
//                //                const size_t j=std::max(a,b);
//                const auto key=LoopLinkType::getKey(a,b);
//                const size_t& i=key.first;
//                const size_t& j=key.second;
//                const IsConstNetworkLinkType Lij(this->link(i,j));
//                assert(Lij.first && "Merging non-existing link");
//                const SharedNodePtrType nA((Lij.second->source->sID==a)? Lij.second->source : Lij.second->sink);
//                const SharedNodePtrType nB((Lij.second->source->sID==b)? Lij.second->source : Lij.second->sink);
//                
//                // Find NetworkLink i1->j1
//                const auto key1=LoopLinkType::getKey(a1,b1);
//                const size_t& i1=key1.first;
//                const size_t& j1=key1.second;
//                //                const size_t i1=std::min(a1,b1);
//                //                const size_t j1=std::max(a1,b1);
//                const IsConstNetworkLinkType Li1j1(this->link(i1,j1));
//                assert(Li1j1.first && "Merging non-existing link");
//                //          const  SharedNodePtrType nA1((Li1j1.second->source()->sID==a1)? Li1j1.second->source : Li1j1.second->sink);
//                //          const  SharedNodePtrType nB1((Li1j1.second->source()->sID==b1)? Li1j1.second->source : Li1j1.second->sink);
//                
//                
//                //                MergeDequeType mergeDeq;
//                AMTdequeType deq;
//                
//                // For each loop link in Li1j1, store merge operations in a MergeMapType
//                for (const auto& loopLink1 : Li1j1.second->loopLinks())
//                {
//                    //                    mergeDeq.push_back(MergeMapType());
//                    
//                    size_t nAnnihilations=0;
//                    
//                    for (const auto& loopLink : Lij.second->loopLinks())
//                    {
//                        if(loopLink1->loop().get()==loopLink->loop().get()) // annihilation
//                        {
//                            
//                            if(loopLink->source()->sID==b && loopLink->sink()->sID==a &&
//                               loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1)
//                            {
//                                //                                mergeDeq.back().emplace(std::piecewise_construct,
//                                //                                                        std::make_tuple(0),
//                                //                                                        std::make_tuple(loopLink1->prev->source(),
//                                //                                                                        loopLink1->source(),
//                                //                                                                        loopLink1->sink(),
//                                //                                                                        loopLink1->next->sink(),
//                                //                                                                        loopLink1->loop(),
//                                //                                                                        nA,
//                                //                                                                        nB,
//                                //                                                                        loopLink->loop())
//                                //                                                        );
//                                
//                                deq.emplace_back(loopLink1->prev->source(),
//                                                 loopLink1->source(),
//                                                 loopLink1->sink(),
//                                                 loopLink1->next->sink(),
//                                                 loopLink1->loop(),
//                                                 nA,
//                                                 nB,
//                                                 loopLink->loop()
//                                                 );
//                                
//                            }
//                            else if(loopLink->source()->sID==a && loopLink->sink()->sID==b &&
//                                    loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1)
//                            {
//                                //                                mergeDeq.back().emplace(std::piecewise_construct,
//                                //                                                        std::make_tuple(0),
//                                //                                                        std::make_tuple(loopLink1->prev->source(),
//                                //                                                                        loopLink1->source(),
//                                //                                                                        loopLink1->sink(),
//                                //                                                                        loopLink1->next->sink(),
//                                //                                                                        loopLink1->loop(),
//                                //                                                                        nB,
//                                //                                                                        nA,
//                                //                                                                        loopLink->loop())
//                                //                                                        );
//                                
//                                deq.emplace_back(loopLink1->prev->source(),
//                                                 loopLink1->source(),
//                                                 loopLink1->sink(),
//                                                 loopLink1->next->sink(),
//                                                 loopLink1->loop(),
//                                                 nB,
//                                                 nA,
//                                                 loopLink->loop()
//                                                 );
//                            }
//                            else
//                            {
//                                assert(0 && "Invalid merge, same link direction.");
//                            }
//                            nAnnihilations++;
//                        }
//                        //                        else // merge or transfer
//                        //                        {
//                        //                            if(loopLink->source()->sID==a && loopLink->sink()->sID==b &&
//                        //                               loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1 &&
//                        //                               loopLink->flow() == loopLink1->flow()) // opposite direction, same flow
//                        //                            {
//                        //                                mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                        std::make_tuple(1),
//                        //                                                        std::make_tuple(loopLink1->prev->source(),
//                        //                                                                        loopLink1->source(),
//                        //                                                                        loopLink1->sink(),
//                        //                                                                        loopLink1->next->sink(),
//                        //                                                                        loopLink1->loop(),
//                        //                                                                        nB,
//                        //                                                                        nA,
//                        //                                                                        loopLink->loop())
//                        //                                                        );
//                        //                            }
//                        //                            else if(loopLink->source()->sID==b && loopLink->sink()->sID==a &&
//                        //                                    loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1 &&
//                        //                                    loopLink->flow() == loopLink1->flow()) // opposite direction, same flow
//                        //                            {
//                        //                                mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                        std::make_tuple(1),
//                        //                                                        std::make_tuple(loopLink1->prev->source(),
//                        //                                                                        loopLink1->source(),
//                        //                                                                        loopLink1->sink(),
//                        //                                                                        loopLink1->next->sink(),
//                        //                                                                        loopLink1->loop(),
//                        //                                                                        nA,
//                        //                                                                        nB,
//                        //                                                                        loopLink->loop())
//                        //                                                        );
//                        //                            }
//                        //                            else if(loopLink->source()->sID==a && loopLink->sink()->sID==b &&
//                        //                                    loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1 &&
//                        //                                    loopLink->flow() == -loopLink1->flow()) // same direction, opposite flow
//                        //                            {
//                        //                                loopLink1->loop()->flip(); // flip second loop
//                        //                                mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                        std::make_tuple(1),
//                        //                                                        std::make_tuple(loopLink1->prev->source(),
//                        //                                                                        loopLink1->source(),
//                        //                                                                        loopLink1->sink(),
//                        //                                                                        loopLink1->next->sink(),
//                        //                                                                        loopLink1->loop(),
//                        //                                                                        nB,
//                        //                                                                        nA,
//                        //                                                                        loopLink->loop())
//                        //                                                        );
//                        //                            }
//                        //                            else if(loopLink->source()->sID==b && loopLink->sink()->sID==a &&
//                        //                                    loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1 &&
//                        //                                    loopLink->flow() == -loopLink1->flow()) // same direction, opposite flow
//                        //                            {
//                        //                                loopLink1->loop()->flip(); // flip second loop
//                        //                                mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                        std::make_tuple(1),
//                        //                                                        std::make_tuple(loopLink1->prev->source(),
//                        //                                                                        loopLink1->source(),
//                        //                                                                        loopLink1->sink(),
//                        //                                                                        loopLink1->next->sink(),
//                        //                                                                        loopLink1->loop(),
//                        //                                                                        nA,
//                        //                                                                        nB,
//                        //                                                                        loopLink->loop())
//                        //                                                        );
//                        //                            }
//                        //                            else // different flow, perform transfer
//                        //                            {
//                        //                                if(loopLink1->source()->sID==a1 && loopLink1->sink()->sID==b1)
//                        //                                {
//                        //                                    mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                            std::make_tuple(2),
//                        //                                                            std::make_tuple(loopLink1->prev->source(),
//                        //                                                                            loopLink1->source(),
//                        //                                                                            loopLink1->sink(),
//                        //                                                                            loopLink1->next->sink(),
//                        //                                                                            loopLink1->loop(),
//                        //                                                                            nA,
//                        //                                                                            nB)
//                        //                                                            );
//                        //
//                        //                                }
//                        //                                else if(loopLink1->source()->sID==b1 && loopLink1->sink()->sID==a1)
//                        //                                {
//                        //                                    mergeDeq.back().emplace(std::piecewise_construct,
//                        //                                                            std::make_tuple(2),
//                        //                                                            std::make_tuple(loopLink1->prev->source(),
//                        //                                                                            loopLink1->source(),
//                        //                                                                            loopLink1->sink(),
//                        //                                                                            loopLink1->next->sink(),
//                        //                                                                            loopLink1->loop(),
//                        //                                                                            nB,
//                        //                                                                            nA)
//                        //                                                            );
//                        //                                }
//                        //                                else
//                        //                                {
//                        //                                    assert(0);
//                        //                                }
//                        //                            }
//                        //                        }
//                        
//                    }
//                    
//                    assert(nAnnihilations<=1);
//                }
//                
//                // Perform merging
//                for(const auto tup : deq)
//                {
//                    //                    for(const auto& pair : mergeMap)
//                    //                    {
//                    //                        const auto& pair(mergeMap)
//                    //                        const int&   op(mergeMap.begin()->first);   // only perform first operation for each element of mergeDeq
//                    //                        const auto& tup(mergeMap.begin()->second);  // only perform first operation for each element of mergeDeq
//                    
//                    //                    assert(op==0);
//                    loopLinkAnnihilate(tup);
//                    
//                    //                        switch (op)
//                    //                        {
//                    //                            case 0: // annhilation
//                    //                            {
//                    //
//                    ////                                if(std::get<0>(tup)->sID==std::get<6>(tup)->sID || std::get<0>(tup)->sID==std::get<5>(tup)->sID)
//                    ////                                {// part of the loop is pinched off, no need to create another loop
//                    ////                                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////                                    connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                    ////                                }
//                    ////                                else if(std::get<3>(tup)->sID==std::get<6>(tup)->sID || std::get<3>(tup)->sID==std::get<5>(tup)->sID)
//                    ////                                {// part of the loop is pinched off, no need to create another loop
//                    ////                                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////                                    connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                    ////                                }
//                    ////                                else
//                    ////                                {
//                    ////                                    // Break the loop in two parts
//                    ////                                    disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    ////                                    disconnect(std::get<6>(tup),std::get<5>(tup),std::get<4>(tup)); // note 6->5, not 5->6 as in merge
//                    ////
//                    ////                                    // links 0->1 and 2->3 should still exist, but now they are disconnected, so reset the loop ptr on on of them
//                    ////                                    IsLoopLinkType link01=loopLink(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    ////                                    assert(link01.first && "LoopLink must exist");
//                    ////
//                    ////                                    IsLoopLinkType link23=loopLink(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////                                    assert(link23.first && "LoopLink must exist");
//                    ////
//                    ////                                    std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->derived(),std::get<4>(tup)->flow());
//                    ////                                    link01.second->resetLoop(tempLoop);
//                    ////
//                    ////                                    // Finish disconnecting
//                    ////                                    disconnect(std::get<0>(tup),std::get<1>(tup),tempLoop); // after resetLoop link01 is in loop tempLoop
//                    ////                                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////
//                    ////                                    connect(std::get<0>(tup),std::get<5>(tup),tempLoop);
//                    ////                                    connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                    ////
//                    ////                                }
//                    //                                break;
//                    //                            }
//                    //
//                    //                            case 1: // loop merge
//                    //                            {
//                    //                                disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    //                                connect(std::get<5>(tup),std::get<6>(tup),std::get<4>(tup));
//                    //
//                    //                                if(std::get<1>(tup)->sID!=std::get<5>(tup)->sID)
//                    //                                {
//                    //                                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    //                                    connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                    //                                }
//                    //
//                    //                                if(std::get<2>(tup)->sID!=std::get<6>(tup)->sID)
//                    //                                {
//                    //                                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //                                    connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //                                }
//                    //
//                    //                                break;
//                    //                            }
//                    //
//                    //                            case 2: // simple merge
//                    //                            {
//                    //                                disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    //                                connect(std::get<5>(tup),std::get<6>(tup),std::get<4>(tup));
//                    //
//                    //                                if(std::get<1>(tup)->sID!=std::get<5>(tup)->sID)
//                    //                                {
//                    //                                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    //                                    connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                    //                                }
//                    //
//                    //                                if(std::get<2>(tup)->sID!=std::get<6>(tup)->sID)
//                    //                                {
//                    //                                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //                                    connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //                                }
//                    //
//                    //                                break;
//                    //                            }
//                    //
//                    //                            default:
//                    //                                assert(0 && "More than one loop links found in other network link.");
//                    //                                break;
//                    //                        }
//                    
//                    
//                    //                if(std::get<7>(tup)) // annihilation
//                    //                {
//                    //                    disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    //                    disconnect(std::get<5>(tup),std::get<6>(tup),std::get<4>(tup));
//                    //
//                    //                    std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->derived(),std::get<4>(tup)->flow());
//                    //
//                    //                    std::get<8>(tup)->resetLoop(tempLoop);
//                    //
//                    //                    connect(std::get<5>(tup),std::get<2>(tup),tempLoop);
//                    //                    connect(std::get<1>(tup),std::get<6>(tup),std::get<4>(tup));
//                    //                }
//                    //                else
//                    //                {
//                    //                    disconnect(std::get<0>(tup),std::get<1>(tup),std::get<4>(tup));
//                    //                    disconnect(std::get<1>(tup),std::get<2>(tup),std::get<4>(tup));
//                    //                    disconnect(std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //
//                    //                    connect(std::get<0>(tup),std::get<5>(tup),std::get<4>(tup));
//                    //                    connect(std::get<5>(tup),std::get<6>(tup),std::get<4>(tup));
//                    //                    connect(std::get<6>(tup),std::get<3>(tup),std::get<4>(tup));
//                    //                }
//                    //                    }
//                }
//                
//            }
//            
//            VerboseLoopNetwork(1,"done merging ("<<a<<","<<b<<") and ("<<a1<<","<<b1<<")"<<std::endl);
//            
//        }
        
        
        /**********************************************************************/
        void checkLoops() const
        {
            std::cout<<"Checking "<<LoopObserverType::loops().size()<<" loops"<<std::endl;
            for(const auto& loop : LoopObserverType::loops())
            {
                assert(loop.second->isLoop() && "not a closed loop");
            }
            
        }
        
        /**********************************************************************/
        void printLoops() const
        {
            for(const auto& loop : LoopObserverType::loops())
            {
                loop.second->printLoop();
            }
            
        }
        
        void printLoopLinks() const
        {
            for(const auto& link : loopLinks())
            {
                std::cout<<"link "<<link.second.source()->sID<<"->"<<link.second.sink()->sID
                <<" (prev "<<link.second.prev->source()->sID<<"->"<<link.second.prev->sink()->sID<<")"
                <<" (next "<<link.second.next->source()->sID<<"->"<<link.second.next->sink()->sID<<")"
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
    
    template<typename Derived>
    int LoopNetwork<Derived>::verboseLevel=1;
    
}
#endif
