/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopLink_H_
#define model_LoopLink_H_

#include <memory>
#include <iterator>
#include <model/MPI/MPIcout.h>

#define VerboseLoopLink(N,x) if(verboseLevel>=N){model::cout<<x;}

namespace model
{
    template<typename LinkType>
    class LoopLink
    {
        
        typedef typename TypeTraits<LinkType>::NodeType LoopNodeType;
        typedef typename TypeTraits<LinkType>::LoopType LoopType;
        typedef typename TypeTraits<LinkType>::FlowType FlowType;

        
    public:
        
        static int verboseLevel;
        
        LoopLink(const LoopLink&) =delete;
        
        //        const LoopNodeType* const source;
        //        const LoopNodeType* const sink;
        const std::shared_ptr<LoopNodeType> source;
        const std::shared_ptr<LoopNodeType> sink;
        
        
        std::shared_ptr<LoopType> pLoop;
        const std::shared_ptr<LinkType> pLink;
        
        LoopLink* prev;
        LoopLink* next;
        
        /**********************************************************************/
        LoopLink(const std::shared_ptr<LoopNodeType>& so,
                 const std::shared_ptr<LoopNodeType>& si,
                 const std::shared_ptr<LoopType>& pL) :
        /* init */ source(so),
        /* init */ sink(si),
        /* init */ pLoop(pL),
        /* init */ pLink(pLoop->loopNetwork.pLink(source,sink)),
        /* init */ prev(nullptr),
        /* init */ next(nullptr)
        {
            VerboseLoopLink(1,"Constructing LoopLink "<<source->sID<<"->"<<sink->sID<<" (loop "<<pLoop->sID<<")"<<std::endl);
            pLoop->addLink(this);
            pLink->addLink(this);
            
            source->addLoopLink(this);
            sink->addLoopLink(this);
        }
        
        /**********************************************************************/
        ~LoopLink()
        {
            VerboseLoopLink(1,"Destroying LoopLink "<<source->sID<<" "<<sink->sID<<" (loop "<<pLoop->sID<<")"<<std::endl);
            pLoop->removeLink(this);
            pLink->removeLink(this);
            
            source->removeLoopLink(this);
            sink->removeLoopLink(this);
            
        }
        
        /**********************************************************************/
        void resetLoop(const std::shared_ptr<LoopType>& pL)
        {
            if(pL.get()!=pLoop.get())
            {
                VerboseLoopLink(1,"LoopLink "<<source->sID<<"->"<<sink->sID<<", resetting loop: old loop="<<pLoop->sID<<std::flush);
                
                pLoop->removeLink(this);
                pLoop=pL;
                pLoop->addLink(this);
                VerboseLoopLink(1,", new loop="<<pLoop->sID<<std::endl);

                if(next!=nullptr)
                {
                    next->resetLoop(pLoop);
                }
                
                if(prev!=nullptr)
                {
                    prev->resetLoop(pLoop);
                }
//                else
//                {
//  
//                }

            }
            
        }
        
        /**********************************************************************/
        const FlowType& flow() const
        {
            return pLoop->flow();
        }
        
        //        /**********************************************************************/
        //        const LoopLink* const previous() const
        //        {
        //            const auto iter=pLoop->links().find(std::make_pair(source->sID,source->sID));
        //            assert(iter!=pLoop->links().end() && "LoopLink not found in Loop.");
        //            return (iter==pLoop->links().begin())? pLoop->links().rbegin()->second : std::prev(iter)->second;
        //        }
        //
        //        /**********************************************************************/
        //        const LoopLink* const next() const
        //        {
        //            const auto iter=pLoop->links().find(std::make_pair(source->sID,source->sID));
        //            assert(iter!=pLoop->links().end() && "LoopLink not found in Loop.");
        //            return (iter==pLoop->links().rbegin())? pLoop->links().begin()->second : std::next(iter)->second;
        //        }
        
    };
    
    template<typename LinkType>
    int LoopLink<LinkType>::verboseLevel=1;
    
}
#endif
