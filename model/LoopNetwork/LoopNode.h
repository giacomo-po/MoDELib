/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopNode_H_
#define model_LoopNode_H_

#include <iostream>
//#include <list>
#include <set>
#include <memory>
#include <assert.h>
#include <tuple>
#include <limits.h>
//#include <iterator>

#include <StaticID.h>
#include <CRTP.h>
//#include <iomanip>
//#include <vector>
//#include <Eigen/Dense>
#include <NetworkComponent.h>
#include <MPIcout.h>
#include <NodeObserver.h>
#include <LoopLink.h>

#define VerboseLoopNode(N,x) if(verboseLevel>=N){model::cout<<x;}


namespace model
{
    
//    template<typename LoopLinkType>
//    struct LoopLinkNodeConnection
//    {// USE THIS CLASS TO STORE LOOPLINKS IN A NODE. MAKE MAP WITH KEY LOOPID, AND THIS AS VALUE
//      
//        LoopType* const loop;
//        LoopLinkType*  inLink;
//        LoopLinkType* outLink;
//        
////        LoopLinkNodeConnection
//        
//    };
    
    
    template<typename Derived>
    class LoopNode : public  StaticID<Derived>,
    /*            */ public  CRTP<Derived>,
    /*            */ private std::set<LoopLink<typename TypeTraits<Derived>::LinkType>*>,
    /*            */ private std::map<size_t,std::tuple<Derived* const ,typename TypeTraits<Derived>::LinkType* const>>
    {
        
    public:
        
        typedef Derived NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef NodeObserver<Derived> NodeObserverType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        typedef std::map<size_t,LoopLinkContainerType> LinkByLoopContainerType;
        typedef std::tuple<Derived* const ,LinkType* const>				NeighborType;
        typedef std::map<size_t,NeighborType>						    NeighborContainerType;

        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;

        friend class NetworkLink<LinkType>; // allow NetworkLink to call private NetworkNode::formNetworkComponent

        LoopNetworkType* const loopNetwork;

        
    private:
        
//        NodeObserverType* const loopNetwork;
        std::shared_ptr<NetworkComponentType> psn;

        /**********************************************************************/
        void resetPSN()
        {
            //! 1- Removes this from the current NetworkComponent
            this->psn->remove(this->p_derived());
            //! 2- Creates a new NetworkComponent containing this
            this->psn.reset(new NetworkComponentType(this->p_derived()));
            //! 3- Transmits 'formNetworkComponent' to the neighbors
            typedef void (Derived::*node_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
            node_member_function_pointer_type Nmfp(&Derived::formNetworkComponent);
            //			Nmfp=&Derived::formNetworkComponent;
            typedef void (LinkType::*link_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
            link_member_function_pointer_type Lmfp(&LinkType::formNetworkComponent);
            //			Lmfp=&LinkType::formNetworkComponent;
            depthFirstExecute(Nmfp,Lmfp,this->psn);
        }
        
        /**********************************************************************/
        void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
        {/*!@param[in] psnOther a shared_ptr to another NetworkComponent
          *
          * If the shared_ptr to the NetworkComponent of *this is different from
          * psnOther, the former is reset (this may destroy the NetworkComponent)
          * and reassigned to psnOther.
          */
            if (psn!=psnOther)
            {
                psn->remove(this->p_derived());
                psn=psnOther;		// redirect psn to the new NetworkComponent
                psn->add(this->p_derived());	// add this in the new NetworkComponent
            }
        }
        
    public:


        
        static int verboseLevel;
        
        
        
    
        /**********************************************************************/
        LoopNode(LoopNetworkType* const ln) :
        /* init list */ loopNetwork(ln),
        /* init list */ psn(new NetworkComponentType(this->p_derived()))
        {
            VerboseLoopNode(1,"Constructing LoopNode "<<tag()<<std::endl);

//            std::cout<<"Constructing LoopNode "<<this->sID<<std::endl;
            loopNetwork->addNode(this->p_derived());
            
//            const bool success=neighbors().emplace(std::make_pair(this->sID, NeighborType(this->p_derived(),(LinkType*) NULL,0) )).second;
//            assert(success && "CANNOT INSERT SELF IN NEIGHBORHOOD.");

        }
        
        /**********************************************************************/
        LoopNode(const LoopNode&) =delete;

        
        /**********************************************************************/
        ~LoopNode()
        {
            VerboseLoopNode(1,"Destroying LoopNode "<<tag()<<std::endl);

            
            loopNetwork->removeNode(this->p_derived());
            
            
            assert(loopLinks().empty());
            
            this->psn->remove(this->p_derived());

            
//            const int success=neighbors().erase(this->sID);
//            assert(success==1 && "CANNOT ERESE SELF FROM NEIGHBORHOOD.");

        }
        
        /**********************************************************************/
        const LoopNetworkType& network() const
        {
            return *loopNetwork;
        }
        
        /**********************************************************************/
        LoopNetworkType& network()
        {
            return *loopNetwork;
        }
        
        /**********************************************************************/
        size_t snID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            return psn->snID(this->p_derived());
        }
        
        /**********************************************************************/
        size_t gID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            return network().globalNodeID(this->sID);
        }
    
        /**********************************************************************/
        const LoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        
        /**********************************************************************/
        LoopLinkContainerType& loopLinks()
        {
            return *this;
        }
        
        /**********************************************************************/
        LinkByLoopContainerType linksByLoopID() const
        {
            LinkByLoopContainerType temp;
            for(const auto& link : loopLinks())
            {
                temp[link->loop()->sID].insert(link);
            }
            return temp;
        }
        
        /**********************************************************************/
        std::set<const LoopType*> loops() const
        {
            std::set<const LoopType*> temp;
            for(const auto& link : loopLinks())
            {
                temp.insert(link->loop().get());
            }
            return temp;
        }
        
        /**********************************************************************/
        const NeighborContainerType& neighbors() const
        {
            return *this;
        }
        
        /**********************************************************************/
        NeighborContainerType& neighbors()
        {
            return *this;
        }
        
        /**********************************************************************/
        bool isSimple() const
        {/*!\returns true if neighbors().size()==2
          */
            return neighbors().size()==2;
        }
        
        /**********************************************************************/
        void addToNeighborhood(LinkType* const pL)
        {/*!@param[in] pL a pointer to a LinkType edge
          */
            if (pL->source->sID==this->sID)
            {// this vertex is the source of edge *pL, so sink of *pL is the neighbor
                const NeighborType temp(pL->sink.get(),pL);
                const bool success=neighbors().emplace(pL->sink->sID,temp).second;
                assert(success && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else if (pL->sink->sID==this->sID)
            {// this vertex is the sink of edge *pL, so source of *pL is the neighbor
                const NeighborType temp(pL->source.get(),pL);
                const bool success=neighbors().emplace(pL->source->sID,temp).second;
                assert(success  && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else
            {
                assert(0 && "CANNOT INSERT NON-INCIDENT EDGE");
            }
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(LinkType* const pL)
        {
            if (pL->source->sID==this->sID)
            {
                const size_t key=pL->sink->sID;
                const int success=neighbors().erase(key);
                assert(success==1);
            }
            else if (pL->sink->sID==this->sID)
            {
                const size_t key=pL->source->sID;
                const int success=neighbors().erase(key);
                assert(success==1);
            }
            else
            {
                assert(0 && "CANNOT REMOVE FROM NEIGHBORS");
            }
        }
        
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
//            std::cout<<"LoopNode "<<this->sID<<" adding Link "<<pL->source->sID<<"->"<<pL->sink->sID<<std::endl;
            
            assert((pL->source().get()==this || pL->sink().get()==this) && "LoopLink connected to wrong node.");

            const bool inserted=loopLinks().insert(pL).second; // Store the connecting link

            assert(inserted);
            
            
            for(const auto& other : loopLinks())
            {// Check loop consistency based on prev/next of each loopLink

                if(pL->loop().get()==other->loop().get() && pL!=other)
                {//two distinct loopLinks of the same loop
                    if(pL->source().get()==other->sink().get())
                    {
                        assert(pL->prev==nullptr || pL->prev==other);
                        assert(other->next==nullptr || other->next==pL);
                        pL->prev=other;
                        other->next=pL;
                    }
                    if (pL->sink().get()==other->source().get())
                    {
                        assert(pL->next==nullptr || pL->next==other);
                        assert(other->prev==nullptr || other->prev==pL);
                        pL->next=other;
                        other->prev=pL;
                    }
                    if(pL->source().get()==other->source().get() || pL->sink().get()==other->sink().get())
                    {
                        std::cout<<"LoopNode "<<this->sID<<std::endl;
                        std::cout<<"adding link "<<pL->source()->sID<<"->"<<pL->sink()->sID<<" (loop "<<pL->loop()->sID<<")"<<std::endl;
                        std::cout<<"existing link "<<other->source()->sID<<"->"<<other->sink()->sID<<" (loop "<<other->loop()->sID<<")"<<std::endl;
                        assert(0 && "Not a Loop");
                    }
                }
            }
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const int erased=loopLinks().erase(pL);
            assert(erased==1);
            
            if(pL->next!=nullptr)
            {
                pL->next->prev=nullptr;
                pL->next=nullptr;
            }

            if(pL->prev!=nullptr)
            {
                pL->prev->next=nullptr;
                pL->prev=nullptr;
            }

        }
        
        /**********************************************************************/
        const std::shared_ptr<NetworkComponentType> & pSN() const
        {/*!\returns A const reference to the shared-pointer to the
          * NetworkComponent containing this.
          */
            return psn;
        }
        
        /**********************************************************************/
        std::string tag() const
        {/*!\returns the string "i" where i is this->sID
          */
            return std::to_string(this->sID) ;
        }
        
        /**********************************************************************/
        template <typename T>
        void depthFirstExecute(void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX)
        {
            std::set<size_t> searchedNodes;
            std::set<std::pair<size_t,size_t> > searchedLinks;
            depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N);
        }
        
        /**********************************************************************/
        template <typename T>
        void depthFirstExecute(std::set<size_t>& searchedNodes,
                               std::set<std::pair<size_t,size_t> >& searchedLinks,
                               void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&),
                               const T & input,
                               const size_t& N = ULONG_MAX)
        {
            (this->p_derived()->*Nfptr)(input); // execute Nfptr on this node
            if (N!=0)
            {
                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
                for (const auto& neighborIter : neighbors())
                {
//                    if (!std::get<2>(neighborIter.second)==0)
//                    {
                        if (searchedLinks.find(std::get<1>(neighborIter.second)->nodeIDPair)==searchedLinks.end())
                        {  // neighbor not searched
                            (std::get<1>(neighborIter.second)->*Lfptr)(input); // execute Lfptr on connecting link
                            assert(searchedLinks.insert(std::get<1>(neighborIter.second)->nodeIDPair).second && "CANNOT INSERT CURRENT LINK IN SEARCHED LINKS"); // this node has been searched
                        }
//                    }
                    if (searchedNodes.find(std::get<0>(neighborIter.second)->sID)==searchedNodes.end())
                    {  // neighbor not searched
                        std::get<0>(neighborIter.second)->depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N-1); // continue executing on neighbor
                    }
                }
            }
        }
        
        
        /**********************************************************************/
        bool depthFirstSearch (const size_t& ID, const size_t& N = ULONG_MAX) const
        {
            std::set<size_t> searchedNodes;
            return depthFirstSearch(searchedNodes,ID,N);
        }
        
        /**********************************************************************/
        bool depthFirstSearch (std::set<size_t>& searchedNodes, const size_t& ID, const size_t& N = ULONG_MAX) const
        {
            bool reached(this->sID==ID);
            if (N!=0 && !reached)
            {
                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
                for (const auto& neighborIter : neighbors())
                {
                    if (searchedNodes.find(std::get<0>(neighborIter.second)->sID)==searchedNodes.end())
                    {  // neighbor not searched
                        reached=std::get<0>(neighborIter.second)->depthFirstSearch(searchedNodes,ID,N-1); // ask if neighbor can reach
                        if (reached)
                        {
                            break;
                        }
                    }
                }
            }
            return reached;
        }
        
    };
    
    template<typename Derived>
    int LoopNode<Derived>::verboseLevel=0;

    
    
}
#endif
