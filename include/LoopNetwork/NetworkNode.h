/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkNode_H_
#define model_NetworkNode_H_

#include <iostream>
#include <set>
#include <memory>
#include <assert.h>
#include <tuple>
#include <limits.h>

#include <StaticID.h>
#include <CRTP.h>
#include <NetworkBase.h>


#ifndef NDEBUG
#define VerboseNetworkNode(N,x) if(verboseLevel>=N){std::cout<<redColor<<x<<defaultColor;}
#else
#define VerboseNetworkNode(N,x)
#endif

namespace model
{
    
    
    template<typename Derived>
    class NetworkNode : public StaticID<Derived>
    /*            */,public CRTP<Derived>
    /*            */,public NetworkBase<Derived,size_t>
    /*            */,public std::set<typename TypeTraits<Derived>::LoopNodeType*>
    /*            */,private std::map<size_t,std::tuple<Derived* const ,typename TypeTraits<Derived>::NetworkLinkType* const>>
    {
        
    public:
        
        typedef Derived NetworkNodeType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef typename std::set<LoopLinkType*> LoopLinkContainerType;
                typedef typename TypeTraits<Derived>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef std::set<LoopNodeType*> LoopNodeContainerType;
        typedef NetworkBase<Derived,size_t> NetworkBaseType;
        typedef std::tuple<Derived* const ,NetworkLinkType* const>                NeighborType;
        typedef std::map<size_t,NeighborType>                            NeighborContainerType;

//        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;

        
    private:
        
//        NEED TO ADD NEXT/PREV

//        NodeObserverType* const loopNetwork;
//        std::shared_ptr<NetworkComponentType> psn;

//        /**********************************************************************/
//        void resetPSN()
//        {
//            //! 1- Removes this from the current NetworkComponent
//            this->psn->remove(this->p_derived());
//            //! 2- Creates a new NetworkComponent containing this
//            this->psn.reset(new NetworkComponentType(this->p_derived()));
//            //! 3- Transmits 'formNetworkComponent' to the neighbors
//            typedef void (Derived::*node_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
//            node_member_function_pointer_type Nmfp(&Derived::formNetworkComponent);
//            //            Nmfp=&Derived::formNetworkComponent;
//            typedef void (LinkType::*link_member_function_pointer_type)(const std::shared_ptr<NetworkComponentType>&);
//            link_member_function_pointer_type Lmfp(&LinkType::formNetworkComponent);
//            //            Lmfp=&LinkType::formNetworkComponent;
//            depthFirstExecute(Nmfp,Lmfp,this->psn);
//        }
        
//        /**********************************************************************/
//        void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
//        {/*!@param[in] psnOther a shared_ptr to another NetworkComponent
//          *
//          * If the shared_ptr to the NetworkComponent of *this is different from
//          * psnOther, the former is reset (this may destroy the NetworkComponent)
//          * and reassigned to psnOther.
//          */
//            if (psn!=psnOther)
//            {
//                psn->remove(this->p_derived());
//                psn=psnOther;        // redirect psn to the new NetworkComponent
//                psn->add(this->p_derived());    // add this in the new NetworkComponent
//            }
//        }
        
    public:

        
        static int verboseLevel;
        
        
        
    
        /**********************************************************************/
        NetworkNode(LoopNetworkType* const loopNetwork_in) :
        /* init */ NetworkBaseType(loopNetwork_in,&loopNetwork_in->networkNodes(),this->sID)
//        /* init list */,psn(network().commonNetworkComponent? network().commonNetworkComponent : std::shared_ptr<NetworkComponentType>(new NetworkComponentType(this->p_derived())))
        {
            VerboseNetworkNode(1,"Constructing NetworkNode "<<tag()<<std::endl);

//            if(network().commonNetworkComponent)
//            {
//                psn->add(this->p_derived());
//            }
        }
        
        /**********************************************************************/
        ~NetworkNode()
        {
            VerboseNetworkNode(1,"Destroying NetworkNode "<<tag()<<std::endl);
            assert(neighbors().empty());
            
//            this->psn->remove(this->p_derived());
            
        }
        
//        size_t networkID() const
//        {
//            return std::distance(this->network().networkNodes().begin(),this->network().networkNodes().find(this->key));
//        }
         LoopLinkContainerType outLoopLinks() const
        {
            LoopLinkContainerType temp;
            for (const auto &ln : loopNodes())
            {
                if (ln->next.second)
                {
                    temp.insert(ln->next.second);
                }
            }
            return temp;
        }

        LoopLinkContainerType inLoopLinks() const
        {
            LoopLinkContainerType temp;
            for (const auto &ln : loopNodes())
            {
                if (ln->prev.second)
                {
                    temp.insert(ln->prev.second);
                }
            }
            return temp;
        }

        std::set<size_t> loopIDs() const
        {
            std::set<size_t> temp;
            for (const auto &loopIter : loops())
            {
                temp.insert(loopIter->sID);
            }
            return temp;
        }

        std::set<LoopType*> loops() const
        {
            std::set<LoopType*> temp;
            for (const auto& ln : loopNodes())
            {
                temp.insert(ln->loop().get());
            }
            return temp;
        }

        const LoopNodeContainerType& loopNodes() const
        {
            return *this;
        }
        
        LoopNodeContainerType& loopNodes()
        {
            return *this;
        }
        
        void addLoopNode(LoopNodeType* const pN)
        {
            const bool success(loopNodes().insert(pN).second);
            if(!success)
            {
                throw std::runtime_error("Duplicate LoopNode could not be added");
            }
//            assert(success && "Duplicate LoopNode could not be added");
        }
        
        void removeLoopNode(LoopNodeType* const pN)
        {
            const size_t erased(loopNodes().erase(pN));
            if(erased!=1)
            {
                throw std::runtime_error("LoopNode could not be erased");
            }
//            assert(erased==1 && "LoopNode could not be erased");
        }
        
        bool isContractableTo(const std::shared_ptr<Derived>&) const
        {
            return true;
        }
        
//        /**********************************************************************/
//        size_t snID() const
//        {/*!\returns The NetworkComponent::snID() of the component
//          * containing this.
//          */
//            return psn->snID(this->p_derived());
//        }
        
        /**********************************************************************/
        size_t gID() const
        {/*!\returns The NetworkComponent::snID() of the component
          * containing this.
          */
            return this->network().globalNodeID(this->sID);
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
        
//        /**********************************************************************/
//        bool isSimple() const
//        {/*!\returns true if neighbors().size()==2
//          */
//            return neighbors().size()==2;
//        }
        
        /**********************************************************************/
        void addToNeighborhood(NetworkLinkType* const pL)
        {/*!@param[in] pL a pointer to a LinkType edge
          */
            if (pL->source->sID==this->sID)
            {// this vertex is the source of edge *pL, so sink of *pL is the neighbor
                const NeighborType temp(pL->sink.get(),pL);
                // std::cout<<" Addng to neighborhood "<<pL->tag()<<std::endl;
                // std::cout<<" Current Neighborhood "<<std::endl;
                // for (const auto& neigh : neighbors())
                // {
                //     std::cout<<std::get<0>(neigh.second)->tag()<<"=>"<<std::get<1>(neigh.second)->tag()<<std::endl;
                // }
                const bool success=neighbors().emplace(pL->sink->sID,temp).second;
                if(!success)
                {
                    throw std::runtime_error("CANNOT INSERT IN NEIGHBORHOOD.");
                }
//                assert(success && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else if (pL->sink->sID==this->sID)
            {// this vertex is the sink of edge *pL, so source of *pL is the neighbor
                const NeighborType temp(pL->source.get(),pL);
                const bool success=neighbors().emplace(pL->source->sID,temp).second;
                if(!success)
                {
                    throw std::runtime_error("CANNOT INSERT IN NEIGHBORHOOD.");
                }
//                assert(success  && "CANNOT INSERT IN NEIGHBORHOOD.");
            }
            else
            {
                throw std::runtime_error("CANNOT INSERT NON-INCIDENT EDGE");
//
//                assert(0 && "CANNOT INSERT NON-INCIDENT EDGE");
            }
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(NetworkLinkType* const pL)
        {
            if (pL->source->sID==this->sID)
            {
                const size_t key=pL->sink->sID;
                const int success=neighbors().erase(key);
                if(success!=1)
                {
                    throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
                }
//                assert(success==1);
            }
            else if (pL->sink->sID==this->sID)
            {
                const size_t key=pL->source->sID;
                const int success=neighbors().erase(key);
                if(success!=1)
                {
                    throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
                }
//                assert(success==1);
            }
            else
            {
                throw std::runtime_error("CANNOT REMOVE FROM NEIGHBORHOOD.");
//
//                assert(0 && "CANNOT REMOVE FROM NEIGHBORS");
            }
        }
        
        
        
//        /**********************************************************************/
//        const std::shared_ptr<NetworkComponentType> & pSN() const
//        {/*!\returns A const reference to the shared-pointer to the
//          * NetworkComponent containing this.
//          */
//            return psn;
//        }
        
        /**********************************************************************/
        std::string tag() const
        {/*!\returns the string "i" where i is this->sID
          */
            return std::to_string(this->sID) ;
        }
        
//        /**********************************************************************/
//        template <typename T>
//        void depthFirstExecute(void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&), const T & input, const size_t& N = ULONG_MAX)
//        {
//            std::set<size_t> searchedNodes;
//            std::set<std::pair<size_t,size_t> > searchedLinks;
//            depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N);
//        }
//
//        /**********************************************************************/
//        template <typename T>
//        void depthFirstExecute(std::set<size_t>& searchedNodes,
//                               std::set<std::pair<size_t,size_t> >& searchedLinks,
//                               void (Derived::*Nfptr)(const T&),void (LinkType::*Lfptr)(const T&),
//                               const T & input,
//                               const size_t& N = ULONG_MAX)
//        {
//            (this->p_derived()->*Nfptr)(input); // execute Nfptr on this node
//            if (N!=0)
//            {
//                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
//                for (const auto& neighborIter : neighbors())
//                {
////                    if (!std::get<2>(neighborIter.second)==0)
////                    {
//                        if (searchedLinks.find(std::get<1>(neighborIter.second)->nodeIDPair)==searchedLinks.end())
//                        {  // neighbor not searched
//                            (std::get<1>(neighborIter.second)->*Lfptr)(input); // execute Lfptr on connecting link
//                            assert(searchedLinks.insert(std::get<1>(neighborIter.second)->nodeIDPair).second && "CANNOT INSERT CURRENT LINK IN SEARCHED LINKS"); // this node has been searched
//                        }
////                    }
//                    if (searchedNodes.find(std::get<0>(neighborIter.second)->sID)==searchedNodes.end())
//                    {  // neighbor not searched
//                        std::get<0>(neighborIter.second)->depthFirstExecute(searchedNodes,searchedLinks,Nfptr,Lfptr,input, N-1); // continue executing on neighbor
//                    }
//                }
//            }
//        }
        
        
//        /**********************************************************************/
//        bool depthFirstSearch (const size_t& ID, const size_t& N = ULONG_MAX) const
//        {
//            std::set<size_t> searchedNodes;
//            return depthFirstSearch(searchedNodes,ID,N);
//        }
//
//        /**********************************************************************/
//        bool depthFirstSearch (std::set<size_t>& searchedNodes, const size_t& ID, const size_t& N = ULONG_MAX) const
//        {
//            bool reached(this->sID==ID);
//            if (N!=0 && !reached)
//            {
//                assert(searchedNodes.insert(this->sID).second && "CANNOT INSERT CURRENT NODE IN SEARCHED NODES"); // this node has been searched
//                for (const auto& neighborIter : neighbors())
//                {
//                    if (searchedNodes.find(std::get<0>(neighborIter.second)->sID)==searchedNodes.end())
//                    {  // neighbor not searched
//                        reached=std::get<0>(neighborIter.second)->depthFirstSearch(searchedNodes,ID,N-1); // ask if neighbor can reach
//                        if (reached)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//            return reached;
//        }
        
    };
    
    template<typename Derived>
    int NetworkNode<Derived>::verboseLevel=0;

    
    
}
#endif
