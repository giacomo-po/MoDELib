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
#include <model/Utilities/CRTP.h>
#include <model/LoopNetwork/LoopLink.h>
#include <model/LoopNetwork/NetworkLinkObserver.h>
#include <model/LoopNetwork/NetworkComponent.h>


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
        typedef typename TypeTraits<Derived>::FlowType FlowType;
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
//            std::cout<<"1"<<std::endl;
            // Add this to NetworkLinkObserver
            loopNetwork->addLink(this->p_derived());

//                        std::cout<<"2"<<std::endl;
            // Add this to neighobors of source and sink
            source->addToNeighborhood(this->p_derived());
            sink  ->addToNeighborhood(this->p_derived());
            
            //! 2 - Joins source and sink NetworkComponents
            if (source->pSN()==sink->pSN()) // source and sink are already in the same NetworkComponent
            {
//                            std::cout<<"3"<<std::endl;
                psn=source->pSN();				// redirect psn to the source psn
                psn->add(this->p_derived());	// add this to the existing NetworkComponent
            }
            else // source and sink are in different NetworkComponents
            {
//                            std::cout<<"4"<<std::endl;
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
//                     std::cout<<"end"<<std::endl;
        }
        
    public:
        
        const std::shared_ptr<NodeType> source;
        const std::shared_ptr<NodeType> sink;
        const std::pair<size_t,size_t> nodeIDPair;
        
        
        /**********************************************************************/
        NetworkLink(const std::shared_ptr<NodeType>& nI,
                    const std::shared_ptr<NodeType>& nJ) :
//        /* init */ _flow(TypeTraits<Derived>::zeroFlow),
        /* init */ loopNetwork(nI->loopNetwork),
        /* init */ source(nI->sID<nJ->sID? nI : nJ),
        /* init */ sink(nI->sID<nJ->sID? nJ : nI),
        /* init */ nodeIDPair(std::make_pair(source->sID,sink->sID))
        {
            //                        std::cout<<"Constructing NetworkLink ("<<source->sID<<","<<sink->sID<<")"<<std::endl;
            //            const bool sourceInserted=source->insert(this->p_derived()).second;
            //            assert(sourceInserted);
            //            const bool sinkInserted=sink->insert(this->p_derived()).second;
            //            assert(sinkInserted);
            assert(nI->loopNetwork==nJ->loopNetwork && "source and sink in different networks");
            makeTopologyChange();
        }
        
        /**********************************************************************/
        ~NetworkLink()
        {
            //                        std::cout<<"Destroying NetworkLink "<<source->sID<<" "<<sink->sID<<std::endl;
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
        void addLink(LoopLinkType* const pL)
        {
            const bool success=loopLinks().insert(pL).second;
            assert(success && "Could not insert LoopLink in NetworkLink");
            
//            if(pL->source()->sID==source->sID)
//            {
//                _flow+=pL->flow();
//            }
//            else
//            {
//                _flow-=pL->flow();
//            }
            
            
        }
        
        /**********************************************************************/
        void removeLink(LoopLinkType* const pL)
        {
            const size_t erased=loopLinks().erase(pL);
            assert(erased==1 && "Could not erase LoopLink from NetworkLink");
            
//            if(pL->source()->sID==source->sID)
//            {
//                _flow-=pL->flow();
//            }
//            else
//            {
//                _flow+=pL->flow();
//            }
        }
        
        /**********************************************************************/
        std::string tag() const
        {
            return std::to_string(source->sID)+"->"+std::to_string(sink->sID);
        }
        
//        const FlowType& flow() const
//        {
//            return _flow;
//        }
        
    };
    
    
}
#endif
