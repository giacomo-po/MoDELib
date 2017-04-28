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


namespace model
{
    template<typename Derived>
    class NetworkLink : public CRTP<Derived>,
    /*               */ private std::set<const LoopLink<Derived>*>
    {
        
        public:
        
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef typename TypeTraits<Derived>::LoopType LoopType;
        typedef LoopLink<Derived> LoopLinkType;
        typedef std::set<const LoopLinkType*> LoopLinkContainerType;
        
        
//        const NodeType* const source;
//        const NodeType* const sink;
        std::shared_ptr<NodeType> source;
        std::shared_ptr<NodeType> sink;

        
        /**********************************************************************/
//        NetworkLink(const NodeType* const nI,
//                    const NodeType* const nJ) :
//        /* init */ source(nI->sID<nJ->sID? nI : nJ),
//        /* init */ sink(nI->sID<nJ->sID? nJ : nI)
        NetworkLink(const std::shared_ptr<NodeType>& nI,
                    const std::shared_ptr<NodeType>& nJ) :
        /* init */ source(nI->sID<nJ->sID? nI : nJ),
        /* init */ sink(nI->sID<nJ->sID? nJ : nI)
        {
//                        std::cout<<"Constructing NetworkLink ("<<source->sID<<","<<sink->sID<<")"<<std::endl;
            NetworkLinkObserver<LinkType>::addLink(this->p_derived());
         
//            const bool sourceInserted=source->insert(this->p_derived()).second;
//            assert(sourceInserted);
//            const bool sinkInserted=sink->insert(this->p_derived()).second;
//            assert(sinkInserted);

        }
        
        /**********************************************************************/
        ~NetworkLink()
        {
//                        std::cout<<"Destroying NetworkLink "<<source->sID<<" "<<sink->sID<<std::endl;
            NetworkLinkObserver<LinkType>::removeLink(this->p_derived());
            
//            const int sourceErased=source->erase(this->p_derived());
//            assert(sourceErased==1);
//            const bool sinkErased=sink->erase(this->p_derived());
//            assert(sinkErased==1);

        }
        
        /**********************************************************************/
        NetworkLink(const NetworkLink&) =delete;
        const NetworkLink& operator=(const NetworkLink&) =delete;
        
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
        void addLink(const LoopLinkType* const pL)
        {
            const bool success=loopLinks().insert(pL).second;
            assert(success && "Could not insert LoopLink in NetworkLink");
        }
        
        /**********************************************************************/
        void removeLink(const LoopLinkType* const pL)
        {
            const size_t erased=loopLinks().erase(pL);
            assert(erased==1 && "Could not erase LoopLink from NetworkLink");
        }
        
    };
    
    
}
#endif
