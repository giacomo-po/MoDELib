/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLink_H_
#define model_NetworkLink_H_

#include <set>
#include <CRTP.h>
#include <LoopLink.h>
#include <NetworkLinkObserver.h>
#include <NetworkComponent.h>


namespace model
{
    template<typename Derived>
    class NetworkLink : public StaticID<Derived>,
    /*               */ public CRTP<Derived>,
    /*               */ private std::set<LoopLink<Derived>*>
    {
        
    public:
        
//        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
//        typedef typename TypeTraits<Derived>::FlowType FlowType;
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        
        typedef LoopLink<Derived> LoopLinkType;
        //        typedef std::set<const LoopLinkType*> LoopLinkContainerType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        
        typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;

        
        friend class LoopNode<NodeType>; // allow NetworkNode to call private NetworkLink::formNetworkComponent

        LoopNetworkType* const loopNetwork;
        
    private:
        

        std::shared_ptr<NetworkComponentType> psn;
//        FlowType _flow;
        
        
        /* formNetworkComponent ***********************************************/
        void formNetworkComponent(const std::shared_ptr<NetworkComponentType> & psnOther)
        {
            if (psn!=psnOther)
            {
                psn->remove(this->p_derived());
                psn=psnOther;		// redirect psn to the new NetworkComponent
                psn->add(this->p_derived());    // add this in the new NetworkComponent
            }
        }
        
        /**********************************************************************/
        void changeSN(const NetworkComponent<NodeType,LinkType>& SN)
        {
            
            const std::map<size_t,NodeType* const> tempNodeMap(SN);
            const std::map<std::pair<size_t,size_t>,LinkType* const> tempLinkMap(SN);
            
            for (typename std::map<size_t,NodeType* const>::const_iterator vIter=tempNodeMap.begin();vIter!=tempNodeMap.end();++vIter)
            {
                vIter->second->formNetworkComponent(psn);
            }
            
            for (typename std::map<std::pair<size_t,size_t>,LinkType* const>::const_iterator lIter=tempLinkMap.begin();lIter!=tempLinkMap.end();++lIter)
            {
                lIter->second->formNetworkComponent(psn);
            }
            
        }
        
        /**********************************************************************/
        void makeTopologyChange()
        {
            // Add this to NetworkLinkObserver
            loopNetwork->addLink(this->p_derived());

            // Add this to neighobors of source and sink
            source->addToNeighborhood(this->p_derived());
            sink  ->addToNeighborhood(this->p_derived());
            
            //! 2 - Joins source and sink NetworkComponents
            if (source->pSN()==sink->pSN()) // source and sink are already in the same NetworkComponent
            {
                psn=source->pSN();				// redirect psn to the source psn
                psn->add(this->p_derived());	// add this to the existing NetworkComponent
            }
            else // source and sink are in different NetworkComponents
            {
                // find the size of the source and sink
                size_t sourceSNsize(source->pSN()->nodeOrder());
                size_t   sinkSNsize(  sink->pSN()->nodeOrder());
                if (sourceSNsize>=sinkSNsize)
                {
                    psn=source->pSN();					   // redirect psn to the source psn
                    psn->add(this->p_derived());		   // add this to the source NetworkComponent
                    changeSN(*(sink->psn.get()));
                }
                else
                {
                    psn=sink->pSN();				       // redirect psn to the sink psn
                    psn->add(this->p_derived());	       // add this to the source NetworkComponent
                    changeSN(*(source->psn.get()));
                }
            }
        }
        
    public:
        
        const std::shared_ptr<NodeType> source;
        const std::shared_ptr<NodeType> sink;
        const std::pair<size_t,size_t> nodeIDPair;
        
        /**********************************************************************/
        NetworkLink(const std::shared_ptr<NodeType>& nI,
                    const std::shared_ptr<NodeType>& nJ) :
        /* init */ loopNetwork(nI->loopNetwork),
        /* init */ source(nI->sID<nJ->sID? nI : nJ),
        /* init */ sink(nI->sID<nJ->sID? nJ : nI),
        /* init */ nodeIDPair(std::make_pair(source->sID,sink->sID))
        {
            assert(nI->loopNetwork==nJ->loopNetwork && "source and sink in different networks");
            makeTopologyChange();
        }
        
        /**********************************************************************/
        ~NetworkLink()
        {
            loopNetwork->removeLink(this->p_derived());
            
            source->removeFromNeighborhood(this->p_derived());
            sink->removeFromNeighborhood(this->p_derived());
            
            this->psn->remove(this->p_derived());
            
            //! 3- If Now Source and Sink are disconnected then reset the NetworkComponent in the sink
            const bool sourceCanReachSink(source->depthFirstSearch(sink->sID));
            
            if (!sourceCanReachSink)
            {
                sink -> resetPSN();
            }

            
        }
        
        /**********************************************************************/
        NetworkLink(const NetworkLink&) =delete;
        const NetworkLink& operator=(const NetworkLink&) =delete;
        
        
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
            assert(success && "Could not insert LoopLink in NetworkLink");
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const size_t erased=loopLinks().erase(pL);
            assert(erased==1 && "Could not erase LoopLink from NetworkLink");
        }
        
        /**********************************************************************/
        std::string tag() const
        {
            return std::to_string(source->sID)+"->"+std::to_string(sink->sID);
        }
        
    };
    
    
}
#endif
