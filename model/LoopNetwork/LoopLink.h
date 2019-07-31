/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopLink_H_
#define model_LoopLink_H_

#include <algorithm>
#include <memory>
#include <string>
#include <iterator>
#include <TypeTraits.h>
#include <MPIcout.h>
#include <NodeObserver.h>


#define VerboseLoopLink(N,x) if(verboseLevel>=N){model::cout<<x;}

namespace model
{
    template<typename LinkType>
    class LoopLink
    {
        
        typedef typename TypeTraits<LinkType>::NodeType LoopNodeType;
        typedef typename TypeTraits<LinkType>::LoopType LoopType;
        typedef typename TypeTraits<LinkType>::FlowType FlowType;
        
        
        const std::shared_ptr<LoopNodeType> _source;
        const std::shared_ptr<LoopNodeType> _sink;
        std::shared_ptr<LoopType> pLoop; // THIS SHOULD BE CONST, THAT WAY WE COULD USE A MAP WITH 3 IDs TO STORE LOOPLINKS
        
    public:
        
        
        typedef std::pair<size_t,size_t> KeyType;
        
        
        /**********************************************************************/
        static KeyType networkLinkKey(const std::shared_ptr<LoopNodeType>& Ni,const std::shared_ptr<LoopNodeType>& Nj)
        {
            return networkLinkKey(Ni->sID,Nj->sID);
        }
        
        /**********************************************************************/
        static KeyType networkLinkKey(const size_t& i,const size_t& j)
        {
            assert(i!=j && "i and j cannot be the same");
            return KeyType(std::min(i,j),std::max(i,j));
        }
        
        static int verboseLevel;
        
        
        //        const LoopNodeType* const source;
        //        const LoopNodeType* const sink;
        
        
        
        const std::shared_ptr<LinkType> pLink;
        
        LoopLink* prev;
        LoopLink* next;
        
        /**********************************************************************/
        LoopLink(const std::shared_ptr<LoopNodeType>& so,
                 const std::shared_ptr<LoopNodeType>& si,
                 const std::shared_ptr<LoopType>& pL) :
        /* init */ _source(so),
        /* init */ _sink(si),
        /* init */ pLoop(pL),
        /* init */ pLink(pLoop->network().sharedLink(_source,_sink)),
        /* init */ prev(nullptr),
        /* init */ next(nullptr)
        {
            VerboseLoopLink(1,"Constructing LoopLink "<<tag()<<" (loop "<<pLoop->sID<<")"<<std::endl);
            pLink->addLoopLink(this);
            _source->addLoopLink(this);
            _sink->addLoopLink(this);
            pLoop->addLoopLink(this);
        }

        /**********************************************************************/
        LoopLink(const LoopLink&) =delete;

        
        /**********************************************************************/
        ~LoopLink()
        {// call removeLoopLink in inverse order compared to addLoopLink
            VerboseLoopLink(1,"Destroying LoopLink "<<tag()<<" (loop "<<pLoop->sID<<")"<<std::endl);
            pLoop->removeLoopLink(this);
            _sink->removeLoopLink(this);
            _source->removeLoopLink(this);
            pLink->removeLoopLink(this);
        }
        
        /**********************************************************************/
        const std::shared_ptr<LoopNodeType>& source() const
        {
            return _source;
        }
        
        /**********************************************************************/
        const std::shared_ptr<LoopNodeType>& sink() const
        {
            return _sink;
        }
        
        /**********************************************************************/
        const std::shared_ptr<LoopType>& loop() const
        {
            return pLoop;
        }
        
        /**********************************************************************/
        void resetLoop(const std::shared_ptr<LoopType>& pL)
        {
            if(pL.get()!=pLoop.get())
            {
                VerboseLoopLink(1,"LoopLink "<<tag()<<", resetting loop: old loop="<<pLoop->sID<<std::flush);
                
                pLoop->removeLoopLink(this);
                pLoop=pL;
                pLoop->addLoopLink(this);
                VerboseLoopLink(1,", new loop="<<pLoop->sID<<std::endl);
                
                if(next!=nullptr)
                {
                    next->resetLoop(pLoop);
                }
                
                if(prev!=nullptr)
                {
                    prev->resetLoop(pLoop);
                }
            }
        }
        
        /**********************************************************************/
        void resetLoop(const std::shared_ptr<LoopType>& pL,
                       const size_t& startID,
                       const size_t& endID)
        {
            if(pL.get()!=pLoop.get())
            {
                
                if(source()->sID==startID)
                {
                    VerboseLoopLink(1,"LoopLink "<<tag()<<", resetting loop: old loop="<<pLoop->sID<<std::flush);
                    pLoop->removeLoopLink(this);
                    pLoop=pL;
                    pLoop->addLoopLink(this);
                    VerboseLoopLink(1,", new loop="<<pLoop->sID<<std::endl);
                    
                    if(sink()->sID!=endID)
                    {
                        next->resetLoop(pLoop,sink()->sID,endID);
                    }
                }
                else
                {
                    if(next!=nullptr)
                    {
                        next->resetLoop(pLoop,startID,endID);
                    }
                    else
                    {
                        assert(0 && "next cannot be nullptr for this function");
                    }
                }
            }
            
        }
        
//        /**********************************************************************/
//        void flip()
//        {/*!Swaps source-sink, and prev-next
//          */
//            _source.swap(_sink);
//            std::swap(prev,next);
//        }
        
        /**********************************************************************/
        const FlowType& flow() const
        {/*!\returns a reference to the loop flow
          */
            return pLoop->flow();
        }
        
        /**********************************************************************/
        std::string tag() const
        {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
          */
            return std::to_string(source()->sID) + "->" + std::to_string(sink()->sID) + " ("+std::to_string(loop()->sID)+")";
        }
        
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LoopLink<LinkType>& ll)
        {
            os  << ll.loop()->sID<<"\t"
            /**/<< ll.source()->sID<<"\t"
            /**/<< ll.sink()->sID<<"\t";
            
            return os;
        }
        
    };
    
    template<typename LinkType>
    int LoopLink<LinkType>::verboseLevel=0;
    
}
#endif
