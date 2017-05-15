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

#include <model/Utilities/StaticID.h>
#include <model/Utilities/CRTP.h>
#include <model/LoopNetwork/LoopObserver.h>


namespace model
{
    template<typename Derived>
    class Loop : public CRTP<Derived>,
    /*        */ public StaticID<Derived>,
    /*        */ private std::map<std::pair<size_t,size_t>,const LoopLink<typename TypeTraits<Derived>::LinkType>* const>
    {
        
    public:

        
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef typename TypeTraits<Derived>::NodeType NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
        typedef LoopLink<LinkType> LoopLinkType;
        typedef std::map<std::pair<size_t,size_t>,const LoopLinkType* const> LoopLinkContainerType;
        typedef std::deque<const LoopLinkType*> LoopLinkSequenceType;
        typedef LoopObserver<Derived> LoopObserverType;
        typedef typename TypeTraits<LinkType>::FlowType FlowType;
        
    private:
        
        LoopLinkSequenceType linkSeq;
        
    public:
        
        const LoopNetworkType& loopNetwork;
        const FlowType flow;

        /**********************************************************************/
        Loop(const LoopNetworkType& loopNetwork_in,
             const FlowType& f) :
        /* init */ loopNetwork(loopNetwork_in),
        /* init */ flow(f)
        {
            std::cout<<"Constructing Loop "<<this->sID<<std::endl;
            LoopObserverType::addLoop(this->p_derived());
        }
        
        /**********************************************************************/
        ~Loop()
        {
            std::cout<<"Destorying Loop "<<this->sID<<std::endl;
            LoopObserverType::removeLoop(this->p_derived());
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
        void addLink(const LoopLinkType* const pL)
        {
            const bool success=links().insert(std::make_pair(std::make_pair(pL->source->sID,pL->sink->sID),pL)).second;
            assert(success && "Could not insert in linkMap");
            
//            if(linkSeq.isempty())
//            {
//                linkSeq.insert(linkSeq.end(),pL);
//            }
            
        }
        
        /**********************************************************************/
        void removeLink(const LoopLinkType* const pL)
        {
            const size_t erased=links().erase(std::make_pair(pL->source->sID,pL->sink->sID));
            assert(erased==1 && "Could not erase from linkMap");
        }

//        /**********************************************************************/
//        const LoopLinkSequenceType& linkSequence() const
//        {
//            return linkSeq;
//        }
        
        /**********************************************************************/
        LoopLinkSequenceType linkSequence() const
        {
            //RECODE THIS USING prev/next
            //typename LoopLinkContainerType::const_iterator iter;
            LoopLinkSequenceType temp;
            const LoopLinkType* pL=links().begin()->second;
            for(size_t k=0;k<links().size();++k)
            {
//                typename LoopLinkSequenceType::const_iterator iter;
//                for(iter=temp.begin();
//                    iter!=temp.end();
//                    ++iter)
//                {
//                    if((*iter)->source->sID==link.second->sink->sID)
//                    {
//                        break;
//                    }
//                }
                temp.push_back(pL);
                pL=pL->next;
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isLoop() const
        {
            std::cout<<"Loop "<<this->sID<<std::endl;
            bool temp=true;
            const LoopLinkSequenceType linkSeq(linkSequence());
            for(typename LoopLinkSequenceType::const_iterator iter=linkSeq.begin();iter!=linkSeq.end();++iter)
            {
                std::cout<<(*iter)->source->sID<<"->"<<(*iter)->sink->sID<<std::endl;
                
                auto next=std::next(iter,1);
                if(next==linkSeq.end())
                {
                    next=linkSeq.begin();
                }
                
                temp*=((*iter)->sink->sID==(*next)->source->sID);
                
            }
            return temp;
        }
        
    };
    
    
}
#endif
