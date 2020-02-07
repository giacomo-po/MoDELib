/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Loop_H_
#define model_Loop_H_

#include <iostream>
#include <list>
#include <deque>
#include <map>
#include <memory>
#include <iterator>

#include <StaticID.h>
#include <CRTP.h>
#include <LoopObserver.h>


namespace model
{
    template<typename Derived>
    class Loop : public CRTP<Derived>,
    /*        */ public StaticID<Derived>,
    /*        */ private std::map<std::pair<size_t,size_t>,LoopLink<typename TypeTraits<Derived>::LinkType>* const>
    {
        
    public:

        
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::map<std::pair<size_t,size_t>,LoopLinkType* const> LoopLinkContainerType;
        typedef std::deque<const LoopLinkType*> LoopLinkSequenceType;
        typedef std::deque<std::pair<std::shared_ptr<NodeType>,std::shared_ptr<NodeType>>> LoopNodeSequenceType;
        typedef LoopObserver<Derived> LoopObserverType;
        typedef typename TypeTraits<LinkType>::FlowType FlowType;
        
    private:
        
        LoopLinkSequenceType linkSeq;
        LoopNetworkType* const _loopNetwork;
        FlowType _flow;

    public:
        

//        const LoopNetworkType& loopNetwork;

        /**********************************************************************/
        Loop(LoopNetworkType* const loopNetwork_in,
             const FlowType& f) :
        /* init */ _loopNetwork(loopNetwork_in),
        /* init */ _flow(f)
//        /* init */ loopNetwork(*_loopNetwork)
        {
//            std::cout<<"Constructing Loop "<<this->sID<<std::endl;
            _loopNetwork->addLoop(this->p_derived());
//            LoopObserverType::addLoop(this->p_derived());
        }
        
//        Loop(const Loop<Derived>& other) =delete;
        
        Loop& operator=(const Loop<Derived>& other) =delete;
        
        /**********************************************************************/
        Loop(const Loop<Derived>& other) :
        /* init */ StaticID<Derived>()
        /* init */,LoopLinkContainerType()
        /* init */,_loopNetwork(other._loopNetwork)
        /* init */,_flow(other._flow)
//        /* init */ loopNetwork(*_loopNetwork)
//        /* init */ loopNetwork(*_loopNetwork)
        {
//            std::cout<<"Copying Loop: "<<this->sID<<std::endl;

//            std::cout<<"Constructing Loop "<<this->sID<<std::endl;
            _loopNetwork->addLoop(this->p_derived());
            //            LoopObserverType::addLoop(this->p_derived());
        }
        
//        SOME CONSTRUCTOR HERE IS NOT DOING addLoop
        
        /**********************************************************************/
        ~Loop()
        {
            
            assert(links().size()==0 && "DESTROYING NON-EMPTY LOOP.");
            
            _loopNetwork->removeLoop(this->p_derived());
//            LoopObserverType::removeLoop(this->p_derived());
        }
        
        /**********************************************************************/
        const LoopNetworkType& network() const
        {
            return *_loopNetwork;
        }
        
        LoopNetworkType& network()
        {
            return *_loopNetwork;
        }
        
//        /**********************************************************************/
//        std::shared_ptr<Derived> clone() const
//        {/* Returns a copy of this loop. The new loop, however, does not contain
//          * any LoopLink.
//          */
//            std::shared_ptr<Derived> temp(new Derived(this->derived()));
//            assert(temp->links().empty());
////            temp->links().clear();
//            
//            return temp;
//        }
        
//        /**********************************************************************/
//        void flip()
//        {
//            _flow*=-1;
//            for(auto link : links())
//            {
//                link.second->flip();
//            }
//        }
        
        /**********************************************************************/
        const FlowType& flow() const
        {
            return _flow;
        }
        
        /**********************************************************************/
        LoopLinkContainerType& links()
        {
            return *this;
        }
        
        /**********************************************************************/
        const LoopLinkContainerType& links() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {
            const bool success=links().insert(std::make_pair(LoopLinkType::networkLinkKey(pL->source()->sID,pL->sink()->sID),pL)).second;
            if(!success)
            {
                std::cout<<"DislocationLoop "<<this->sID<<" cannot add LoopLink "<<pL->tag()<<std::endl;
                exit(EXIT_FAILURE);
            }
//            assert(success && "Could not insert in linkMap");
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {
            const size_t erased=links().erase(LoopLinkType::networkLinkKey(pL->source()->sID,pL->sink()->sID));
            assert(erased==1 && "Could not erase from linkMap");
        }

        
        /**********************************************************************/
        LoopLinkSequenceType linkSequence() const
        {
            //RECODE THIS USING prev/next
            //typename LoopLinkContainerType::const_iterator iter;
            LoopLinkSequenceType temp;
            if(links().size())
            {
                const LoopLinkType* pL=links().begin()->second;
                for(size_t k=0;k<links().size();++k)
                {
                    if(pL)
                    {
                        if(pL->next)
                        {
                            temp.push_back(pL);
                            pL=pL->next;
                        }
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        LoopNodeSequenceType nodeSequence() const
        {
            LoopNodeSequenceType temp;
            for(const auto& link : linkSequence())
            {
                temp.emplace_back(link->source(),link->sink());
            }
            return temp;
        }
        
        /**********************************************************************/
        std::deque<size_t> nodeIDSequence() const
        {
            std::deque<size_t> temp;
            for(const auto& link : linkSequence())
            {
                temp.emplace_back(link->source()->sID);
            }
            return temp;
        }
        
        /**********************************************************************/
        std::pair<bool,LoopLinkType*> linkStartingAt(const size_t& i) const
        {
            
            std::pair<bool,LoopLinkType*> temp=std::make_pair(false,static_cast<LoopLinkType*>(nullptr));
            for(const auto& link : links())
            {
                if(link.second->source()->sID==i)
                {
                    temp=std::make_pair(true,link.second);
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isLoop() const
        {
            const LoopLinkSequenceType linkSeq(linkSequence());
            bool temp(linkSeq.size()>=3);
            for(typename LoopLinkSequenceType::const_iterator iter=linkSeq.begin();iter!=linkSeq.end();++iter)
            {
                auto next=std::next(iter,1);
                if(next==linkSeq.end())
                {
                    next=linkSeq.begin();
                }
                
                temp*=((*iter)->sink()->sID==(*next)->source()->sID);
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }

        /**********************************************************************/
        bool isIsolated() const
        {
            bool temp(true);
            for(const auto& link : links())
            {
                temp*=(link.second->pLink->loopLinks().size()==1);
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        void printLoop() const
        {
            std::cout<<"Loop "<<this->sID<<std::endl;
            const LoopLinkSequenceType linkSeq(linkSequence());
            for(typename LoopLinkSequenceType::const_iterator iter=linkSeq.begin();iter!=linkSeq.end();++iter)
            {
                std::cout<<"    "<<(*iter)->source()->sID<<"->"<<(*iter)->sink()->sID
                <<" (prev "<<(*iter)->prev->source()->sID<<"->"<<(*iter)->prev->sink()->sID<<")"
                <<" (next "<<(*iter)->next->source()->sID<<"->"<<(*iter)->next->sink()->sID<<")"<<std::endl;
            }
        }
        
    };
    
    
}
#endif
