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

#include <model/Utilities/CRTP.h>
#include <model/MPI/MPIcout.h>
#include <model/LoopNetwork/LoopNode.h>
//#include <model/LoopNetwork/NetworkNode.h>
#include <model/LoopNetwork/LoopLink.h>
#include <model/LoopNetwork/NetworkLink.h>
#include <model/LoopNetwork/NetworkLinkObserver.h>
#include <model/LoopNetwork/Loop.h>
#include <model/LoopNetwork/LoopObserver.h>
#include <model/LoopNetwork/NetworkComponent.h>
#include <model/LoopNetwork/NetworkComponentObserver.h>
#include <model/Threads/EqualIteratorRange.h>


#define VerboseLoopNetwork(N,x) if(verboseLevel>=N){model::cout<<x;}

namespace model
{
    template<typename Derived>
    class LoopNetwork : public CRTP<Derived>,
    //    /*               */ private std::map<size_t,typename TypeTraits<Derived>::NodeType>,
    /*               */ public LoopObserver<typename TypeTraits<Derived>::LoopType>,
    /*               */ public NetworkComponentObserver<NetworkComponent<typename TypeTraits<Derived>::NodeType,typename TypeTraits<Derived>::LinkType>>,
    /*               */ public NodeObserver<typename TypeTraits<Derived>::NodeType>,
    /*               */ public NetworkLinkObserver<typename TypeTraits<Derived>::LinkType>,
    /*               */ private std::map<size_t,std::shared_ptr<typename TypeTraits<Derived>::NodeType>>,
    /*               */ private std::multimap<std::pair<size_t,size_t>,
    /*                                     */ LoopLink<typename TypeTraits<Derived>::LinkType>,
    /*                                     */ std::less<std::pair<size_t,size_t>>,
    /*                                     */ Eigen::aligned_allocator<std::pair<std::pair<size_t,size_t>, LoopLink<typename TypeTraits<Derived>::LinkType>>>
    /*                                     */ >
    {
        
    public:
        
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<typename TypeTraits<Derived>::LinkType> LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::FlowType FlowType;
        
//        typedef std::multimap<std::pair<size_t,size_t>,LoopLinkType> LoopLinkContainerType;
        typedef std::multimap<std::pair<size_t,size_t>,
        /*                 */ LoopLink<typename TypeTraits<Derived>::LinkType>,
        /*                 */ std::less<std::pair<size_t,size_t>>,
        /*                 */ Eigen::aligned_allocator<std::pair<std::pair<size_t,size_t>, LoopLink<typename TypeTraits<Derived>::LinkType>>>
        /*                 */ > LoopLinkContainerType;
        typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
        typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;
        typedef typename NetworkLinkObserverType::IsNetworkLinkType IsNetworkLinkType;
        typedef typename NetworkLinkObserverType::IsConstNetworkLinkType IsConstNetworkLinkType;
        typedef LoopObserver<LoopType> LoopObserverType;
        
        typedef NodeObserver<NodeType> NodeObserverType;
//        typedef std::shared_ptr<NodeType> SharedNodePtrType;
        typedef typename NodeObserverType::SharedNodePtrType SharedNodePtrType;
        typedef std::map<size_t,SharedNodePtrType>     DanglingNodeContainerType;
        
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
        DanglingNodeContainerType& danglingNodes()
        {
            return *this;
        }
        
        const DanglingNodeContainerType& danglingNodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        IsSharedNodeType danglingNode(const size_t & k)
        {/*!\returns A <bool,NodeType* const> pair, where pair.first is true
          * if node k is in the network, in which case pair.second is a pointer
          * the the node
          */
            typename DanglingNodeContainerType::iterator nodeIter(danglingNodes().find(k));
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
        std::shared_ptr<LinkType> sharedLink(const SharedNodePtrType& nI, const SharedNodePtrType& nJ) const
        {/*!\param[in] nI node I
          * \param[in] nJ node J
          * \returns a shared_ptr to a NetworkLink nI->nJ (or nJ->nI). If no link
          * exists, a new link is created and returned via shared_ptr.
          */
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
                    return nI->sID<nJ->sID? std::shared_ptr<LinkType>(new LinkType(nI,nJ)) : std::shared_ptr<LinkType>(new LinkType(nJ,nI));
//                    return nI->sID<nJ->sID? std::make_shared<LinkType>(nI,nJ) : std::make_shared<LinkType>(nJ,nI);
                }
            }
        }
        
        /**********************************************************************/
        template <typename ...NodeArgTypes>
        std::pair<typename DanglingNodeContainerType::iterator,bool> insertDanglingNode(const NodeArgTypes&... nodeInput)
        {/*! @param[in] nodeInput
          *\returns
          *  Inserts a new vertex in the Network using nodeInput as variable
          *  constructor arguments
          */
            const size_t nodeID(StaticID<NodeType>::nextID());
            
            const std::pair<typename DanglingNodeContainerType::iterator,bool> inserted=danglingNodes().emplace(std::piecewise_construct,
                                                                                                            std::make_tuple(nodeID),
                                                                                                            std::make_tuple(new NodeType(this->p_derived(),nodeInput...)) );
            
            assert(inserted.second && "COULD NOT INSERT NETWORK VERTEX IN VERTEX CONTAINER.");
            assert(inserted.first->first == nodeID && "KEY != nodeID");
            assert(inserted.first->second->sID == nodeID && "sID != nodeID");
            return inserted;
            
        }
        
        /**********************************************************************/
        template <typename ...LoopArgTypes>
        std::shared_ptr<LoopType> insertLoop(const std::vector<size_t> nodeIDs,
//                        const FlowType& f,
                        const LoopArgTypes&... loopInput)
        {/*!@param[in] nodeIDs IDs of the nodes in the loop
          * @param[in] f the loop flow
          * @param[loopInput] additional Loop constructor arguments
          *
          * Inserts a Loop connecting the sequence of nodes with IDs nodeIDs,
          * flow f. The additional loop constructor arguments loopInput
          * are forwarded to the loop constructor.
          */
            std::shared_ptr<LoopType> tempLoop=std::make_shared<LoopType>(this->p_derived(),loopInput...);
            
            for(size_t k=0;k<nodeIDs.size();++k)
            {
                const size_t next=k+1<nodeIDs.size()? k+1 : 0;
                
                IsSharedNodeType n0=danglingNode(nodeIDs[k]);
                IsSharedNodeType n1=danglingNode(nodeIDs[next]);
                assert(n0.first && "Node not found");
                assert(n1.first && "Node not found");
                
                connect(n0.second,n1.second,tempLoop);
            }
            
            assert(tempLoop->isLoop() && "Not a loop.");
            
            return tempLoop;
        }
        
        /**********************************************************************/
        void deleteLoop(const std::shared_ptr<LoopType>& pL)
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
            
            // Find NetworkLink i->j
            const IsConstNetworkLinkType Lij(this->link(i,j));
            assert(Lij.first && "Expanding non-existing link");
            
            // Create new node
            SharedNodePtrType newNode=SharedNodePtrType(new NodeType(*Lij.second,Args...));
            
            // Store what needs to be connected (i,new,j,Loop), or (j,new,i,Loop). This also holds temporarily disconnected nodes
            std::deque<ExpandTupleType> expandDeq;
            for(const auto& llink : Lij.second->loopLinks())
            {
                expandDeq.emplace_back(llink->source(),llink->sink(),llink->loop());
            }
            
            // Delete all LoopLinks of type i->j or j->i
            loopLinks().erase(key);
            
            for (const auto& tup : expandDeq)
            {
                connect(std::get<0>(tup),newNode, std::get<2>(tup));
                connect(newNode,std::get<1>(tup), std::get<2>(tup));
            }
            
            return newNode;
        }
        
        /**********************************************************************/
        void cutLoop(const size_t& L,
                     const size_t& i,const size_t &j)
        {
            const auto loop=this->loop(L);
            
            if(loop.first && i!=j)
            {
                const auto linkAtI=loop.second->linkStartingAt(i);
                const auto linkAtJ=loop.second->linkStartingAt(j);
                
                if(linkAtI.first && linkAtJ.first)
                {
                    if(linkAtI.second->sink()->sID!=j && linkAtJ.second->sink()->sID!=i)
                    {
                        VerboseLoopNetwork(1,"cutting loop "<<L<<" at "<<i<<" and "<<j<<std::endl);

                        //std::shared_ptr<LoopType> tempLoop=loop.second->clone();
                        std::shared_ptr<LoopType> tempLoop(new LoopType(*loop.second));


                        linkAtI.second->resetLoop(tempLoop,i,j);
                        
                        linkAtI.second->prev->next=nullptr;
                        linkAtI.second->prev=nullptr;

                        linkAtJ.second->prev->next=nullptr;
                        linkAtJ.second->prev=nullptr;

                        connect(linkAtJ.second->source(),linkAtI.second->source(),tempLoop);
                        connect(linkAtI.second->source(),linkAtJ.second->source(),linkAtJ.second->loop());
                    }
                }
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
                
                success=true;
            }
            
            return success;
        }
        
        /* execute ************************************************************/
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
                std::cout<<"link "<<link.second.source()->sID<<"->"<<link.second.sink()->sID
                <<" (prev "<<link.second.prev->source()->sID<<"->"<<link.second.prev->sink()->sID<<")"
                <<" (next "<<link.second.next->source()->sID<<"->"<<link.second.next->sink()->sID<<")"
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
