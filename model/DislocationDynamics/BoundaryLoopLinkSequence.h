/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury <ypachaur@purdue.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundaryLoopLinkSequence_H_
#define model_BoundaryLoopLinkSequence_H_

#include <memory>
#include <deque>
#include <TypeTraits.h>
#include <PeriodicDislocationLoopPair.h>

namespace model
{
    template<typename LoopType>
    struct BoundaryLoopLinkSequence : public std::deque<std::deque<const typename TypeTraits<LoopType>::LoopLinkType*>>
    {
        
        const LoopType* loop;
        const std::set<size_t> faceIDs;
        const std::shared_ptr<PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>> loopPair;
        
        /**********************************************************************/
        BoundaryLoopLinkSequence(const LoopType* loop_in,
                                 const std::set<size_t>& faceIDs_in) :
        /* init */ loop(loop_in)
        /* init */,faceIDs(faceIDs_in)
        /* init */,loopPair(new PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>())
        {
            loopPair->insertLoop(this);
        }
        
        /**********************************************************************/
        BoundaryLoopLinkSequence(const LoopType* loop_in,
                                 const std::set<size_t>& faceIDs_in,
                                 const std::shared_ptr<PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>>& loopPair_in) :
        /* init */ loop(loop_in)
        /* init */,faceIDs(faceIDs_in)
        /* init */,loopPair(loopPair_in)
        {
            loopPair->insertLoop(this);
        }
        
        /**********************************************************************/
        ~BoundaryLoopLinkSequence()
         {
             loopPair->eraseLoop(this);
         }
        
    };
}
#endif
