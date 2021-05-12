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
#include <LoopNode.h>
//#include <NetworkNode.h>
#include <LoopLink.h>
#include <NetworkLink.h>
#include <NetworkLinkObserver.h>
#include <Loop.h>
#include <LoopObserver.h>
#include <NetworkComponent.h>
#include <NetworkComponentObserver.h>
#include <EqualIteratorRange.h>


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
    class LoopNetwork : public CRTP<Derived>,
    //    /*               */ private std::map<size_t,typename TypeTraits<Derived>::NodeType>,
    /*               */ public LoopObserver<typename TypeTraits<Derived>::LoopType>,
    /*               */ public NetworkComponentObserver<NetworkComponent<typename TypeTraits<Derived>::NodeType,typename TypeTraits<Derived>::LinkType>>,
    /*               */ public NodeObserver<typename TypeTraits<Derived>::NodeType>,
    /*               */ public NetworkLinkObserver<typename TypeTraits<Derived>::LinkType>,
//    /*               */ private std::map<size_t,std::shared_ptr<typename TypeTraits<Derived>::NodeType>>,
//    /*               */ private std::multimap<std::pair<size_t,size_t>,
//    /*                                     */ LoopLink<typename TypeTraits<Derived>::LinkType>,
//    /*                                     */ std::less<std::pair<size_t,size_t>>
//    /*                                     */ >
    /*               */ private std::map<typename LoopLink<typename TypeTraits<Derived>::LinkType>::KeyType,
    /*                                     */ LoopLink<typename TypeTraits<Derived>::LinkType>>
    {
        
    public:
        
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<typename TypeTraits<Derived>::LinkType> LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::FlowType FlowType;
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
        
        //        typedef std::multimap<std::pair<size_t,size_t>,LoopLinkType> LoopLinkContainerType;
//        typedef std::multimap<std::pair<size_t,size_t>,
//        /*                 */ LoopLink<typename TypeTraits<Derived>::LinkType>,
//        /*                 */ std::less<std::pair<size_t,size_t>>
//        /*                 */ > LoopLinkContainerType;
        typedef std::map<typename LoopLink<typename TypeTraits<Derived>::LinkType>::KeyType,
        /*                                     */ LoopLink<typename TypeTraits<Derived>::LinkType>> LoopLinkContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        typedef typename NetworkLinkObserverType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename NetworkLinkObserverType::IsConstNetworkLinkType IsConstNetworkLinkType;
        typedef LoopObserver<LoopType> LoopObserverType;
        
        typedef NodeObserver<NodeType> NodeObserverType;
        //        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef typename NodeObserverType::SharedNodePtrType SharedNodePtrType;
//        typedef std::map<size_t,SharedNodePtrType>     DanglingNodeContainerType;
        
        //        typedef std::pair<bool,SharedNodePtrType>          IsNodeType;
        typedef typename NodeObserverType::IsSharedNodeType          IsSharedNodeType;
        //        typedef std::pair<bool,const std::shared_ptr<const NodeType>>	 IsConstNodeType;
        
        typedef std::pair<bool,LoopLinkType*> IsLoopLinkType;
        
        typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>> ExpandTupleType;
        
        //        typedef std::tuple<SharedNodePtrType, // A node, preceding B in loop being merged
        //        /*              */ SharedNodePtrType, // B node in loop being merged, which will become node E
        //        /*              */ SharedNodePtrType, // C node in loop being merged, which will become node F
        //        /*              */ SharedNodePtrType, // D node, after C in loop being merged
        //        /*              */ std::shared_ptr<LoopType>, // loop being merged
        //        /*              */ SharedNodePtrType, // E node, will replace B
        //        /*              */ SharedNodePtrType, // F node, will replace C
        //        /*              */ std::shared_ptr<LoopType> // loop being merged
        //        > AMTtupleType;
        //
        //
        //        //        typedef std::map<int,AMTtupleType> MergeMapType; // map is sorted by operation priority
        //        typedef std::deque<AMTtupleType> AMTdequeType;
        
    private:
        
        
        
        

        
        
        
        /**********************************************************************/
        void connect(const SharedNodePtrType& n0,
                     const SharedNodePtrType& n1,
                     const std::shared_ptr<LoopType>& tempLoop)
        {
            VerboseLoopNetwork(1,"connecting "<<n0->sID<<"->"<<n1->sID<<" ("<<tempLoop->sID<<")"<<std::endl);
            //            assert(n0->sID!=n1->sID && "Cannot connect a node to itself");
            
            if(n0->sID!=n1->sID)
            {
                //check if n1->n0 exist for the same loop
                //            const auto iterPair = loopLinks().equal_range(std::pair<size_t,size_t>(n1->sID,n0->sID));

                const auto key=LoopLinkType::loopLinkKey(tempLoop,n0,n1); // even for LoopLinks we store with key (minID,maxID)
//                const auto iterPair = loopLinks().equal_range(key);
                
                const auto linkIter(loopLinks().find(key));
                if(linkIter==loopLinks().end())
                {// link not found
                    loopLinks().try_emplace(key,n0,n1, tempLoop);
                }
                else
                {// link found
                    if(linkIter->second.source()->sID==n1->sID && linkIter->second.sink()->sID==n0->sID)
                    {// If opposite direction disconnect
                        loopLinks().erase(linkIter);
                    }
                }
                
//                const auto key=LoopLinkType::loopLinkKey(n0,n1); // even for LoopLinks we store with key (minID,maxID)
//                const auto iterPair = loopLinks().equal_range(key);
//
//                typename LoopLinkContainerType::const_iterator loopIter=iterPair.second;
//                for (loopIter=iterPair.first;loopIter!=iterPair.second;++loopIter)
//                {
//                    if(loopIter->second.loop().get()==tempLoop.get() && // link in same loop
//                       loopIter->second.source()->sID==n1->sID && loopIter->second.sink()->sID==n0->sID) // opposite link exists
//                    {
//                        break;
//                    }
//
//                }
//
//                if(loopIter!=iterPair.second)
//                {// opposite link exists, disconnect it
//                    loopLinks().erase(loopIter);
//                }
//                else
//                {// opposite link does not exists, connect n0->n1
//                    loopLinks().emplace(std::piecewise_construct,
//                                        std::make_tuple(key.first,key.second),
//                                        std::make_tuple(n0,n1, tempLoop)
//                                        );
//                }
            }
            
        }
        
        /**********************************************************************/
        IsLoopLinkType loopLink(const SharedNodePtrType& n0, const SharedNodePtrType& n1, const std::shared_ptr<LoopType>& loop)
        {
            
            IsLoopLinkType temp=IsLoopLinkType(false,nullptr);
            
            const auto key=LoopLinkType::loopLinkKey(loop,n0,n1);
            const auto linkIter(loopLinks().find(key));
            if(linkIter!=loopLinks().end())
            {// link found
                temp=IsLoopLinkType(true,&linkIter->second);
            }
//            const auto iterPair = loopLinks().equal_range(key);
//            for(typename LoopLinkContainerType::iterator loopIter=iterPair.first;
//                /*                                          */ loopIter!=iterPair.second;
//                /*                                          */ loopIter++)
//            {
//                if(loopIter->second.loop().get()==loop.get() &&
//                   loopIter->second.source()->sID==n0->sID && loopIter->second.sink()->sID==n1->sID)
//                {
//                    temp=IsLoopLinkType(true,&loopIter->second);
//                    break;
//                }
//            }
            return temp;
        }
        
        /**********************************************************************/
        size_t disconnect(const SharedNodePtrType& n0, const SharedNodePtrType& n1, const std::shared_ptr<LoopType>& loop)
        {
            VerboseLoopNetwork(1,"disconnecting "<<n0->sID<<"->"<<n1->sID<<" ("<<loop->sID<<")"<<std::endl);
//            size_t nDisconnected=0;
//            const auto key=LoopLinkType::loopLinkKey(n0,n1);
//            const auto iterPair = loopLinks().equal_range(key);
//
//            for(typename LoopLinkContainerType::const_iterator loopIter=iterPair.first;
//                /*                                          */ loopIter!=iterPair.second;
//                /*                                          */ loopIter++)
//            {
//                if(loopIter->second.loop().get()==loop.get() &&
//                   loopIter->second.source()->sID==n0->sID && loopIter->second.sink()->sID==n1->sID)
//                {
//                    loopLinks().erase(loopIter);
//                    nDisconnected=1;
//                    nDisconnected+=disconnect(n0,n1,loop); // make sure that recursive calls don't find another link
//                    assert(nDisconnected==1 && "More than one LoopLink with same key and same loop exist.");
//                    break;
//                }
//            }
//            return nDisconnected;
            return loopLinks().erase(LoopLinkType::loopLinkKey(loop,n0,n1));
        }
        
//        /**********************************************************************/
//        size_t disconnect(const size_t& sourceID, const size_t& sinkID, const size_t& loopID)
//        {
//            VerboseLoopNetwork(1,"disconnecting "<<sourceID<<"->"<<sinkID<<" ("<<loopID<<")"<<std::endl);
////            size_t nDisconnected=0;
////            const auto key=LoopLinkType::loopLinkKey(sourceID,sinkID);
////            const auto iterPair = loopLinks().equal_range(key);
////
////            for(typename LoopLinkContainerType::const_iterator loopIter=iterPair.first;
////                /*                                          */ loopIter!=iterPair.second;
////                /*                                          */ loopIter++)
////            {
////                if(loopIter->second.loop()->sID==loopID &&
////                   loopIter->second.source()->sID==sourceID && loopIter->second.sink()->sID==sinkID)
////                {
////                    loopLinks().erase(loopIter);
////                    nDisconnected=1;
////                    nDisconnected+=disconnect(sourceID,sinkID,loopID); // make sure that recursive calls don't find another link
////                    assert(nDisconnected==1 && "More than one LoopLink with same key and same loop exist.");
////                    break;
////                }
////            }
////            return nDisconnected;
//            return loopLinks().erase(LoopLinkType::loopLinkKey(loopID,sourceID,sinkID));
//        }
        
        
        /**********************************************************************/
        bool replaceNodeInLoop(const std::shared_ptr<NodeType>& nodeToBeReplaced,
                               const std::shared_ptr<LoopType>& loopToModify,
                               const std::deque<std::shared_ptr<NodeType>>& nodesToBeInserted)
        {
            
            const auto linkByLoopID(nodeToBeReplaced->linksByLoopID());
            const auto loopIter(linkByLoopID.find(loopToModify->sID));
            bool temp(false);
            if(loopIter!=linkByLoopID.end())
            {
                if(loopIter->second.size()!=2)
                {
                    std::cout<<"FATAL ERROR. LoopNetwork::replaceNodeInLoop "<<nodeToBeReplaced->sID<<std::endl;
                    std::cout<<"loop "<<loopIter->first<<" has "<<loopIter->second.size()<<" links"<<std::endl;
                    for(const auto& temp : loopIter->second)
                    {
                        std::cout<<temp->tag()<<std::endl;
                    }
                    std::cout<<"EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                LoopLinkType* const link0(*loopIter->second.begin());
                LoopLinkType* const link1(*loopIter->second.rbegin());
                
                const std::shared_ptr<NodeType> startNode(link0->sink()->sID==nodeToBeReplaced->sID? link0->source() : link1->source());
                const std::shared_ptr<NodeType>   endNode(link0->sink()->sID==nodeToBeReplaced->sID? link1   ->sink(): link0  ->sink());

                disconnect(startNode,nodeToBeReplaced,loopToModify);
                disconnect(nodeToBeReplaced,endNode,loopToModify);
                
                
                for(size_t k=0;k<nodesToBeInserted.size();++k)
                {
                    std::shared_ptr<NodeType> currentSource(k==0? startNode : nodesToBeInserted[k-1]);
                    std::shared_ptr<NodeType> currentSink(k==nodesToBeInserted.size()-1? endNode : nodesToBeInserted[k]);
                    connect(currentSource,currentSink,loopToModify);
                }

                assert(loopToModify->isLoop());
                temp=true;
            }
            return temp;
        }
        
        /**********************************************************************/
        bool remove(const std::shared_ptr<NodeType>& node)
        {/*!\param[in] nodeID the StaticID of the node to be removed
          * \returns true if the node is succesfully removed.
          */
            
            
            
            //            const auto isNode=this->node(nodeID);
            //            if(isNode.first)
            //            {
            const auto linkByLoopID=node->linksByLoopID();
            std::vector<std::tuple<SharedNodePtrType,SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>>> disconnectVector;
            //                std::vector<std::tuple<const SharedNodePtrType&,const SharedNodePtrType&,const SharedNodePtrType&,const std::shared_ptr<LoopType>&>> disconnectVector;
            //                std::vector<std::tuple<size_t,size_t,size_t,size_t>> disconnectVector;
            
            
            for(auto& pair : linkByLoopID)
            {// store what needs to be disconnected and re-connected
                const auto& set(pair.second);
                if(set.size()!=2)
                {
                    std::cout<<"FATAL ERROR. LoopNetwork::remove "<<node->sID<<std::endl;
                    std::cout<<"loop "<<pair.first<<" has "<<set.size()<<" links"<<std::endl;
                    for(const auto& temp : set)
                    {
                        std::cout<<temp->tag()<<std::endl;
                    }
                    std::cout<<"EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                LoopLinkType* const link0(*set.begin());
                LoopLinkType* const link1(*set.rbegin());
                
                if(link0->source().get()==node.get() && link1->sink().get()==node.get())
                {
                    SharedNodePtrType n0(link1->source());
                    SharedNodePtrType n1(link1->sink());
                    SharedNodePtrType n2(link0->sink());
                    std::shared_ptr<LoopType> loop(link0->loop());
                    disconnectVector.emplace_back(n0,n1,n2,loop);
                    //                        disconnectVector.emplace_back(n0->sID,n1->sID,n2->sID,loop->sID);
                    
                    //                        disconnect(n0,n1,loop);
                    //                        disconnect(n1,n2,loop);
                    //                        connect(n0,n2,loop);
                }
                else if(link0->sink().get()==node.get() && link1->source().get()==node.get())
                {
                    SharedNodePtrType n0(link0->source());
                    SharedNodePtrType n1(link0->sink());
                    SharedNodePtrType n2(link1->sink());
                    std::shared_ptr<LoopType> loop(link1->loop());
                    disconnectVector.emplace_back(n0,n1,n2,loop);
                    //                        disconnectVector.emplace_back(n0->sID,n1->sID,n2->sID,loop->sID);
                    
                    //                        disconnect(n0,n1,loop);
                    //                        disconnect(n1,n2,loop);
                    //                        connect(n0,n2,loop);
                }
                else
                {
                    std::cout<<"link0 "<<link0->tag()<<", nodeID="<<node->sID<<std::endl;
                    std::cout<<"link1 "<<link1->tag()<<", nodeID="<<node->sID<<std::endl;
                    assert(false && "Links must connect to node being removed");
                }
            }
            
            for(const auto& tup : disconnectVector)
            {// perform disconnection and re-connection
                disconnect(std::get<0>(tup),std::get<1>(tup),std::get<3>(tup));
                disconnect(std::get<1>(tup),std::get<2>(tup),std::get<3>(tup));
                if(std::get<1>(tup) && std::get<2>(tup) && std::get<3>(tup))
                {
                    connect(std::get<0>(tup),std::get<2>(tup),std::get<3>(tup));
                    if(std::get<3>(tup)->links().size())
                    {
//                        connect(std::get<0>(tup),std::get<2>(tup),std::get<3>(tup));
                        assert(std::get<3>(tup)->isLoop());
                    }
                }
            }
            
            //                for(auto iter=disconnectVector.begin();iter!=disconnectVector.end();)
            //                {// perform disconnection and re-connection
            //                    disconnect(std::get<0>(*iter),std::get<1>(*iter),std::get<3>(*iter));
            //                    disconnect(std::get<1>(*iter),std::get<2>(*iter),std::get<3>(*iter));
            //                    connect(std::get<0>(*iter),std::get<2>(*iter),std::get<3>(*iter));
            //                    iter=disconnectVector.erase(iter); // erase element in disconnectVector to destroy shared_ptr(s)
            //                }
            
            //            }
            
            
//            return !this->node(node->sID).first;
            return true;
        }
        
    public:
        
        static int verboseLevel;
        
        const std::shared_ptr<NetworkComponentType> commonNetworkComponent;
        
//        commonNetworkComponent
        LoopNetwork(const std::shared_ptr<NetworkComponentType>& nc) :
        /* init */ commonNetworkComponent(nc)
        {

        }
        
//        /**********************************************************************/
//        DanglingNodeContainerType& danglingNodes() __attribute__ ((deprecated)) // INSERT USING VECTOR OF SHARED PTR INSTEAD
//        {
//            return *this;
//        }
//
//        const DanglingNodeContainerType& danglingNodes() const __attribute__ ((deprecated)) // INSERT USING VECTOR OF SHARED PTR INSTEAD
//        {
//            return *this;
//        }
        
//        /**********************************************************************/
//        IsSharedNodeType danglingNode(const size_t & k) __attribute__ ((deprecated)) // INSERT USING VECTOR OF SHARED PTR INSTEAD
//        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
//          * if node k is in the network, in which case pair.second is a pointer
//          * the the node
//          */
//            typename DanglingNodeContainerType::iterator nodeIter(danglingNodes().find(k));
//            return (nodeIter!=danglingNodes().end())? std::make_pair(true,nodeIter->second) : std::make_pair(false,SharedNodePtrType(nullptr));
//        }
//
//        /**********************************************************************/
//        void clearDanglingNodes() __attribute__ ((deprecated)) // INSERT USING VECTOR OF SHARED PTR INSTEAD
//        {
//            danglingNodes().clear();
//        }
        
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
        template <typename ...LoopArgTypes>
        std::shared_ptr<LoopType> insertLoop(const std::vector<std::shared_ptr<NodeType>> nodes,
                                             const LoopArgTypes&... loopInput)
        {/*!@param[in] nodes in the loop
          * @param[loopInput] Loop constructor arguments
          *
          * Inserts a Loop connecting the sequence of nodes,
          * The loop constructor arguments loopInput
          * are forwarded to the loop constructor.
          */
            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->p_derived(),loopInput...);
            
            for(size_t k=0;k<nodes.size();++k)
            {
                const size_t next=k+1<nodes.size()? k+1 : 0;
                connect(nodes[k],nodes[next],tempLoop);
            }
            
            assert(tempLoop->isLoop() && "Not a loop.");
            
            tempLoop->update();
            return tempLoop;
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
            
//            assert(this->loops().find(loopID)==this->loops().end() && "LOOP WAS NOT REMOVED");
        }
        
//        /**********************************************************************/
//        template <typename ...LoopArgTypes>
//        std::shared_ptr<LoopType>& reconnectLoop(std::vector<LoopType>& lp,
//                                             const std::vector<std::shared_ptr<NodeType>> nodes)
//        {/*!@param[in] nodes in the loop
//          *
//          * Inserts a Loop connecting the sequence of nodes,
//          */
//            deleteLoop(lp->sID);
//            
//            for(size_t k=0;k<nodes.size();++k)
//            {
//                const size_t next=k+1<nodes.size()? k+1 : 0;
//                connect(nodes[k],nodes[next],lp);
//            }
//            
//            assert(lp->isLoop() && "Not a loop.");
//            
//            return lp;
//        }


        /**********************************************************************/
        SharedNodePtrType expand(const LinkType* const Lij, const SharedNodePtrType& newNode)
        {
            VerboseLoopNetwork(1,"expanding "<<Lij->tag()<<std::endl);
            
            // Store what needs to be connected (i,new,j,Loop), or (j,new,i,Loop). This also holds temporarily disconnected nodes
            std::deque<ExpandTupleType> expandDeq;
            for(const auto& llink : Lij->loopLinks())
            {
                expandDeq.emplace_back(llink->source(),llink->sink(),llink->loop());
            }
            
            for(const auto& tup : expandDeq)
            {
                disconnect(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup));
                connect(std::get<0>(tup),newNode, std::get<2>(tup));
                connect(newNode,std::get<1>(tup), std::get<2>(tup));
            }
            return newNode;
        }
        
        /**********************************************************************/
        SharedNodePtrType expand(const IsConstNetworkLinkType Lij, const SharedNodePtrType& newNode)
        {
            assert(Lij.first && "Expanding non-existing link");
            return expand(Lij.second,newNode);
        }
        
        /**********************************************************************/
        SharedNodePtrType expand(const size_t& a, const size_t& b, const SharedNodePtrType& newNode)
        {
            const IsConstNetworkLinkType Lij(this->link(std::min(a,b),std::max(a,b)));
            return expand(Lij,newNode);
        }
        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        SharedNodePtrType expand(const size_t& a, const size_t& b, const NodeArgTypes&... Args)
        {
            const IsConstNetworkLinkType Lij(this->link(std::min(a,b),std::max(a,b)));
            return expand(Lij,SharedNodePtrType(new NodeType(*Lij.second,Args...)));
        }
        
        /**********************************************************************/
        SharedNodePtrType expandLoopLink(const LoopLinkType* const Lij, const SharedNodePtrType& newNode)
        {
            const SharedNodePtrType source(Lij->source());
            const SharedNodePtrType sink(Lij->sink());
            const std::shared_ptr<LoopType> loop(Lij->loop());
            disconnect(source,sink,loop);
            connect(source,newNode,loop);
            connect(newNode,sink,loop);
            return newNode;
        }
        
        
        
//        /**********************************************************************/
//        void cutLoop(const size_t& L,
//                     const size_t& i,const size_t &j)
//        {
//
//            VerboseLoopNetwork(1,"cutting loop "<<L<<" at "<<i<<","<<j<<std::endl);
//
//
//            const auto loop=this->loop(L);
//
//            if(loop.first && i!=j)
//            {
//                const auto linkAtI=loop.second->linkStartingAt(i);
//                const auto linkAtJ=loop.second->linkStartingAt(j);
//
//                if(linkAtI.first && linkAtJ.first)
//                {
//                    if(linkAtI.second->sink()->sID!=j && linkAtJ.second->sink()->sID!=i)
//                    {
//                        VerboseLoopNetwork(1,"cutting loop "<<L<<" at "<<i<<" and "<<j<<std::endl);
//
//
//                        std::vector<std::tuple<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>,std::shared_ptr<LoopType>>> linksToDisconect;
//
//                        auto currLink(linkAtI.second);
//                        while(currLink->source()->sID!=j)
//                        {
//                            linksToDisconect.emplace_back(currLink->source(),currLink->sink(),currLink->loop());
//                            currLink=currLink->next;
//                        }
//
//                        //std::shared_ptr<LoopType> tempLoop=loop.second->clone();
//                        std::shared_ptr<LoopType> tempLoop(new LoopType(*loop.second));
//
//                        for(const auto& tup : linksToDisconect)
//                        {
//                            disconnect(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup));
//                            connect(std::get<0>(tup),std::get<1>(tup),tempLoop);
//                        }
//
//
//
////                        linkAtI.second->resetLoop(tempLoop,i,j);
////                        linkAtI.second->prev->next=nullptr;
////                        linkAtI.second->prev=nullptr;
////                        linkAtJ.second->prev->next=nullptr;
////                        linkAtJ.second->prev=nullptr;
//
////                        connect(linkAtJ.second->source(),linkAtI.second->source(),tempLoop);
////                        connect(linkAtI.second->source(),linkAtJ.second->source(),linkAtJ.second->loop());
//
//                        connect(std::get<1>(linksToDisconect.back()),std::get<0>(linksToDisconect.front()),tempLoop);
//                        connect(std::get<0>(linksToDisconect.front()),std::get<1>(linksToDisconect.back()),linkAtJ.second->loop());
//
//
//                    }
//                }
//            }
//
//
//        }
        
        /**********************************************************************/
        void cutLoop(const std::shared_ptr<LoopType>& loop,
                     const size_t& i,const size_t &j)
        {
            VerboseLoopNetwork(1,"cutting loop "<<loop->sID<<" at "<<i<<","<<j<<std::endl);
            if(i!=j)
            {
                const auto linkAtI=loop->linkStartingAt(i);
                const auto linkAtJ=loop->linkStartingAt(j);
                
                if(linkAtI.first && linkAtJ.first)
                {
                    if(linkAtI.second->sink()->sID!=j && linkAtJ.second->sink()->sID!=i)
                    {
                        auto currLink(linkAtI.second);
                        if(currLink->sink()->sID!=j)
                        {
                            std::vector<std::tuple<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>>> linksToDisconect;
                            while(currLink->source()->sID!=j)
                            {
                                linksToDisconect.emplace_back(currLink->source(),currLink->sink());
                                currLink=currLink->next;
                            }
                            
                            std::shared_ptr<LoopType> newLoop(new LoopType(*loop));
                            
                            for(const auto& tup : linksToDisconect)
                            {
                                disconnect(std::get<0>(tup),std::get<1>(tup),loop);
                                connect(std::get<0>(tup),std::get<1>(tup),newLoop);
                            }
                            
                            connect(std::get<1>(linksToDisconect.back()),std::get<0>(linksToDisconect.front()),newLoop);
                            connect(std::get<0>(linksToDisconect.front()),std::get<1>(linksToDisconect.back()),loop);
                        }
                    }
                }
            }
        }
        
        /**********************************************************************/
        bool contractSecond(const SharedNodePtrType& nA,const SharedNodePtrType& nB)
        {            
            const size_t a(nA->sID);
            const size_t b(nB->sID);
            VerboseLoopNetwork(1,"contracting "<<a<<","<<b<<std::endl);
            
            // Collect IDs of all loops passing through both a and b
//            std::set<size_t> loopIDs;
            std::set<std::shared_ptr<LoopType>> cutLoops;
            for(const auto& loopLinkA : nA->loopLinks())
            {
                for(const auto& loopLinkB : nB->loopLinks())
                {
                    if(loopLinkA->loop()==loopLinkB->loop())
                    {
                        cutLoops.insert(loopLinkA->loop());
                    }
                }
            }
            
            // Cut those loops
            for(const auto& loop : cutLoops)
            {
                cutLoop(loop,a,b);
            }
            
            
//            std::deque<std::tuple<const LoopLinkType* const,const LoopLinkType* const,int>> mergeDeque;
//            for(const auto& neighborA : nA->neighbors())
//            {
//                for(const auto& neighborB : nB->neighbors())
//                {
//                    if(std::get<0>(neighborA.second)==std::get<0>(neighborB.second))
//                    {
//                        for(const auto& linkA : std::get<1>(neighborA.second)->loopLinks())
//                        {
//                            for(const auto& linkB : std::get<1>(neighborB.second)->loopLinks())
//                            {
//                                if(linkA->loop()!=linkB->loop())
//                                {
//                                    if(linkA->source()==linkB->source() && NullFlow<FlowType>::isZero(linkA->loop()->flow()+linkB->loop()->flow()))
//                                    {// links in same sense
//                                        if(linkA->loop()->isMergeable(linkB->loop()))
//                                        {
//                                            mergeDeque.emplace_back(linkA,linkB,0);
//                                            break;
//                                        }
//                                    }
//                                    if(linkA->sink()==linkB->sink() && NullFlow<FlowType>::isZero(linkA->loop()->flow()+linkB->loop()->flow()))
//                                    {// links in same sense
//                                        if(linkA->loop()->isMergeable(linkB->loop()))
//                                        {
//                                            mergeDeque.emplace_back(linkA,linkB,1);
//                                            break;
//                                        }
//                                    }
//                                    
//                                    if(linkA->source()==linkB->sink() && NullFlow<FlowType>::isZero(linkA->loop()->flow()-linkB->loop()->flow()))
//                                    {// links in opposite sense
//                                        if(linkA->loop()->isMergeable(linkB->loop()))
//                                        {
//                                            mergeDeque.emplace_back(linkA,linkB,2);
//                                            break;
//                                        }
//                                    }
//                                    if(linkA->sink()==linkB->source() && NullFlow<FlowType>::isZero(linkA->loop()->flow()-linkB->loop()->flow()))
//                                    {// links in opposite sense
//                                        if(linkA->loop()->isMergeable(linkB->loop()))
//                                        {
//                                            mergeDeque.emplace_back(linkA,linkB,3);
//                                            break;
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//            
//            VerboseLoopNetwork(2,"contracting, mergeDeque.size()="<<mergeDeque.size()<<std::endl);
//
//            for(const auto& tup : mergeDeque)
//            {
//                auto sourceA(std::get<0>(tup)->source());
//                auto sinkA(std::get<0>(tup)->sink());
//                auto loopA(std::get<0>(tup)->loop());
//                disconnect(sourceA,sinkA,loopA);
//                
//                switch (std::get<2>(tup))
//                {
//                    case 0:
//                    {
//                        VerboseLoopNetwork(2,"contracting, merging case 0"<<std::endl);
//                        connect(std::get<1>(tup)->sink(),sinkA,loopA);
//                        auto currentLinkB(std::get<1>(tup));
//                        while(currentLinkB->source()!=std::get<1>(tup)->sink())
//                        {
//                            currentLinkB=currentLinkB->prev;
//                            connect(currentLinkB->sink(),currentLinkB->source(),loopA);
//                        }
//                        break;
//                    }
//                    case 1:
//                    {
//                        VerboseLoopNetwork(2,"contracting, merging case 1"<<std::endl);
//                        connect(sourceA,std::get<1>(tup)->source(),loopA);
//                        auto currentLinkB(std::get<1>(tup));
//                        while(currentLinkB->source()!=std::get<1>(tup)->sink())
//                        {
//                            currentLinkB=currentLinkB->prev;
//                            connect(currentLinkB->sink(),currentLinkB->source(),loopA);
//                        }
//                        break;
//                    }
//                    case 2:
//                    {
//                        VerboseLoopNetwork(2,"contracting, merging case 2"<<std::endl);
//                        connect(std::get<1>(tup)->source(),sinkA,loopA);
//                        auto currentLinkB(std::get<1>(tup));
//                        while(currentLinkB->sink()!=std::get<1>(tup)->source())
//                        {
//                            currentLinkB=currentLinkB->next;
//                            connect(currentLinkB->source(),currentLinkB->sink(),loopA);
//                        }
//                        break;
//                    }
//                    case 3:
//                    {
//                        VerboseLoopNetwork(2,"contracting, merging case 3"<<std::endl);
//                        connect(sourceA,std::get<1>(tup)->sink(),loopA);
//                        auto currentLinkB(std::get<1>(tup));
//                        while(currentLinkB->sink()!=std::get<1>(tup)->source())
//                        {
//                            currentLinkB=currentLinkB->next;
//                            connect(currentLinkB->source(),currentLinkB->sink(),loopA);
//                        }
//                        break;
//                    }
//                    default:
//                    {
//                        break;
//                    }
//                }
//                
//                VerboseLoopNetwork(2,"contracting, deleting merged loop"<<std::endl);
//                deleteLoop(std::get<1>(tup)->loop()->sID);
//            }
            
            
            // Store links connected to b
            typedef std::tuple<SharedNodePtrType,SharedNodePtrType,std::shared_ptr<LoopType>,bool> ContractTupleType;
            typedef std::deque<ContractTupleType> ContractDequeType;
            ContractDequeType contractDeq;
            
            for(const auto& loopLink : nB->loopLinks())
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
                    connect(std::get<0>(tup),nA,std::get<2>(tup));
                }
                else
                {
                    connect(nA,std::get<1>(tup),std::get<2>(tup));
                }
            }
            
            return true;
        }
        
        
        /**********************************************************************/
        bool contractSecond(const size_t& a,const size_t &b)
        {
            bool success=false;
            const auto nA=this->sharedNode(a);
            const auto nB=this->sharedNode(b);
            
            if(nA.first && nB.first)
            {
                success=contractSecond(nA.second,nB.second);
            }
            
            return success;
        }

        /**********************************************************************/
        bool remove(const size_t& nodeID)
        {/*!\param[in] nodeID the StaticID of the node to be removed
          * \returns true if the node is succesfully removed.
          */
            bool success=false;
            const auto isNode=this->sharedNode(nodeID);
            if(isNode.first)
            {
                success=remove(isNode.second);
            }
            return success;
        }
        
        /**********************************************************************/
        void parallelExecute(void (LinkType::*Lfptr)(void))
        {
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
            EqualIteratorRange<typename NetworkLinkObserverType::LinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
            
#pragma omp parallel for
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                {
                    (linkIter->second->*Lfptr)();
                }
            }
#else
            for (auto& linkIter : this->links())
            {
                (linkIter.second->*Lfptr)();
            }
#endif
        }
        
        /* execute ************************************************************/
        template <typename T>
        void parallelExecute(void (LinkType::*Lfptr)(const T&), const T & input)
        {
#ifdef _OPENMP
            const size_t nThreads = omp_get_max_threads();
            EqualIteratorRange<typename NetworkLinkObserverType::LinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
            
#pragma omp parallel for
            for (size_t thread=0;thread<eir.size();thread++)
            {
                for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
                {
                    (linkIter->second->*Lfptr)(input);
                }
            }
#else
            for (auto& linkIter : this->links())
            {
                (linkIter.second->*Lfptr)(input);
            }
#endif
        }
        
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
        
        /**********************************************************************/
        void printLoopLinks() const
        {
            for(const auto& link : loopLinks())
            {
                std::cout<<"link "<<link.second.tag()
                <<" (prev "<<link.second.prev->tag()<<")"
                <<" (next "<<link.second.next->tag()<<")"
                <<std::endl;
            }
        }
        
        /**********************************************************************/
        void printNodes() const
        {
            for(const auto& node : this->nodes())
            {
                std::cout<<"node "<<node.second->sID<<" neighbor-size="<<node.second->loopLinks().size()<<std::endl;
            }
        }
    };
    
    template<typename Derived>
    int LoopNetwork<Derived>::verboseLevel=0;
    
}
#endif
