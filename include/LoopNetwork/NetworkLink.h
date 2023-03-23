/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLink_H_
#define model_NetworkLink_H_

#include <StaticID.h>
#include <set>
#include <CRTP.h>
#include <NetworkBase.h>
#include <LoopLink.h>

//#include <NetworkLinkObserver.h>
//#include <NetworkComponent.h>


#ifndef NDEBUG
#define VerboseNetworkLink(N,x) if(verboseLevel>=N){std::cout<<magentaColor<<x<<defaultColor;}
#else
#define VerboseNetworkLink(N,x)
#endif

namespace model
{
    template<typename Derived>
    class NetworkLink : public StaticID<Derived>
    /*               */,public CRTP<Derived>
    /*               */,public NetworkBase<Derived,std::pair<size_t,size_t>>
    /*               */,private std::set<typename TypeTraits<Derived>::LoopLinkType*>
    {
        
    public:
        
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<Derived>::NetworkNodeType NetworkNodeType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef NetworkBase<Derived,std::pair<size_t,size_t>> NetworkBaseType;

        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        
//        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;

        
//        friend class LoopNode<NodeType>; // allow NetworkNode to call private NetworkLink::formNetworkComponent

    private:
        

//        std::shared_ptr<NetworkComponentType> psn;
        
        
        /* formNetworkComponent ***********************************************/
//        void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
//        {
//            if (psn!=psnOther)
//            {
//                psn->remove(this->p_derived());
//                psn=psnOther;        // redirect psn to the new NetworkComponent
//                psn->add(this->p_derived());    // add this in the new NetworkComponent
//            }
//        }
//
//        /**********************************************************************/
//        void changeSN(const NetworkComponent<NodeType,LinkType>& SN)
//        {
//
//            const std::map<size_t,NodeType* const> tempNodeMap(SN);
//            const std::map<std::pair<size_t,size_t>,LinkType* const> tempLinkMap(SN);
//
//            for (typename std::map<size_t,NodeType* const>::const_iterator vIter=tempNodeMap.begin();vIter!=tempNodeMap.end();++vIter)
//            {
//                vIter->second->formNetworkComponent(psn);
//            }
//
//            for (typename std::map<std::pair<size_t,size_t>,LinkType* const>::const_iterator lIter=tempLinkMap.begin();lIter!=tempLinkMap.end();++lIter)
//            {
//                lIter->second->formNetworkComponent(psn);
//            }
//
//        }
//
//        /**********************************************************************/
//        void makeTopologyChange()
//        {
//            // Add this to NetworkLinkObserver
//            loopNetwork->addLink(this->p_derived());
//
//            // Add this to neighobors of source and sink
//
//            //! 2 - Joins source and sink NetworkComponents
//            if (source->pSN()==sink->pSN()) // source and sink are already in the same NetworkComponent
//            {
//                psn=source->pSN();                // redirect psn to the source psn
//                psn->add(this->p_derived());    // add this to the existing NetworkComponent
//            }
//            else // source and sink are in different NetworkComponents
//            {
//                // find the size of the source and sink
//                size_t sourceSNsize(source->pSN()->nodeOrder());
//                size_t   sinkSNsize(  sink->pSN()->nodeOrder());
//                if (sourceSNsize>=sinkSNsize)
//                {
//                    psn=source->pSN();                       // redirect psn to the source psn
//                    psn->add(this->p_derived());           // add this to the source NetworkComponent
//                    changeSN(*(sink->psn.get()));
//                }
//                else
//                {
//                    psn=sink->pSN();                       // redirect psn to the sink psn
//                    psn->add(this->p_derived());           // add this to the source NetworkComponent
//                    changeSN(*(source->psn.get()));
//                }
//            }
//        }
        
    public:
        
        const std::shared_ptr<NetworkNodeType> source;
        const std::shared_ptr<NetworkNodeType> sink;
        
        static int verboseLevel;

        
        /**********************************************************************/
        NetworkLink(LoopNetworkType* const loopNetwork,
                    const std::shared_ptr<NetworkNodeType>& nI,
                    const std::shared_ptr<NetworkNodeType>& nJ) :
        /* init */ NetworkBaseType(loopNetwork,&loopNetwork->networkLinks(),getKey(nI,nJ))
        /* init */,source(nI->sID<nJ->sID? nI : nJ)
        /* init */,sink(nI->sID<nJ->sID? nJ : nI)
        {
            VerboseNetworkLink(1,"Constructing NetworkLink "<<tag()<<std::endl);

            if (source==sink)
            {
                throw std::runtime_error(" Source and Sinks are equal for the node "+std::to_string(source->sID));
//                std::cout<<" Source and Sinks are equal for the node "<<source->sID<<std::endl;
            }
  //          assert(source!=sink && "Source and Sink cannot be same for a network link");
            source->addToNeighborhood(this->p_derived());
            sink->addToNeighborhood(this->p_derived());

            //            makeTopologyChange();
        }
        
        /**********************************************************************/
        ~NetworkLink()
        {
            VerboseNetworkLink(1,"Destroying NetworkLink "<<tag()<<std::endl);

            if (source!=sink)
            {
                source->removeFromNeighborhood(this->p_derived());
                sink->removeFromNeighborhood(this->p_derived());
            }
            
//            this->psn->remove(this->p_derived());
            
//            //! 3- If Now Source and Sink are disconnected then reset the NetworkComponent in the sink
//            if(!network().commonNetworkComponent)
//            {
//                const bool sourceCanReachSink(source->depthFirstSearch(sink->sID));
//
//                if (!sourceCanReachSink)
//                {
//                    sink -> resetPSN();
//                }
//            }

            
        }

        std::set<LoopType *> loops() const
        {
            std::set<LoopType *> temp;
            for (const auto &ll : loopLinks())
            {
                temp.insert(ll->loop.get());
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

        /**********************************************************************/
        static typename NetworkBaseType::KeyType getKey(const LoopLinkType* const pL)
        {
            return getKey(pL->source->networkNode,pL->sink->networkNode);
        }

        /**********************************************************************/
        static typename NetworkBaseType::KeyType getKey(const std::shared_ptr<NetworkNodeType>& Ni,
                                                        const std::shared_ptr<NetworkNodeType>& Nj)
        {
            return getKey(Ni->sID,Nj->sID);
        }
        
        /**********************************************************************/
        static typename NetworkBaseType::KeyType getKey(const size_t& i,
                               const size_t& j)
        {
            return typename NetworkBaseType::KeyType(std::min(i,j),std::max(i,j));
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
        void addLoopLink(LoopLinkType* const pL)
        {
            const bool success=loopLinks().insert(pL).second;
            if(!success)
            {
                throw std::runtime_error("Could not insert LoopLink in NetworkLink");
            }
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const size_t erased=loopLinks().erase(pL);
            if(erased!=1)
            {
                throw std::runtime_error("Could not erase LoopLink from NetworkLink");
            }
        }
        
        /**********************************************************************/
        std::string tag() const
        {
            return std::to_string(source->sID)+"->"+std::to_string(sink->sID);
        }
        
    };
    
    template<typename Derived>
    int NetworkLink<Derived>::verboseLevel=0;

}
#endif
