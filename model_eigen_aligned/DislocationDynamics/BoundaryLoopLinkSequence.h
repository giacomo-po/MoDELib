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
//#include <PeriodicDislocationLoopPair.h>

namespace model
{
    template<typename LoopType>
    struct BoundaryLoopLinkSequence : public std::deque<std::deque<const typename TypeTraits<LoopType>::LoopLinkType*>>
    {
        
        const LoopType* loop;
        const std::set<size_t> faceIDs;
//        const std::shared_ptr<PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>> loopPair;
        
        /**********************************************************************/
        BoundaryLoopLinkSequence(const LoopType* loop_in,
                                 const std::set<size_t>& faceIDs_in) :
        /* init */ loop(loop_in)
        /* init */,faceIDs(faceIDs_in)
//        /* init */,loopPair(new PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>())
        {
 //           loopPair->insertLoop(this);
        }
        
//        /**********************************************************************/
//        BoundaryLoopLinkSequence(const LoopType* loop_in,
//                                 const std::set<size_t>& faceIDs_in,
//                                 const std::shared_ptr<PeriodicDislocationLoopPair<BoundaryLoopLinkSequence<LoopType>>>& loopPair_in) :
//        /* init */ loop(loop_in)
//        /* init */,faceIDs(faceIDs_in)
//        /* init */,loopPair(loopPair_in)
//        {
//            assert(loopPair->size()==1);
//            loopPair->insertLoop(this);
//        }
        
//        /**********************************************************************/
//        ~BoundaryLoopLinkSequence()
//         {
//             std::cout<<"Destroying BoundaryLoopLinkSequence of loop "<<loop->sID<<std::endl;
//             std::cout<<"current loopPair->size()="<<loopPair->size()<<std::endl;
//             loopPair->eraseLoop(this);
//             std::cout<<"new loopPair->size()="<<loopPair->size()<<std::endl;
//
//             print();
//
//         }
        
        
//        void imageLinkAtSequenceStart(const std::deque<const typename TypeTraits<LoopType>::LoopLinkType*>& otherSeq,
//                                 const std::set<size_t>& otherFaceIDs) const
//        {
//            
//            if(otherSeq.size())
//            {
//                if(!otherSeq.front()->masterNode)
//                {
//                    for(const auto& otherLink : otherSeq)
//                    {
//                        const auto& otherSource(otherLink->source());
//                        //                if(!otherSource->masterNode)
//                        //                {// start of otherSeq is not an image
//                        const auto& otherSourceImage(otherSource->sharedImage(otherFaceIDs));
//                        bool found(false);
//
//                        if(otherSourceImage->meshFaceIDs()==faceIDs)
//                        {
//                            
//                            for(const auto& linkDeq : *this)
//                            {
//                                for(const auto& link : linkDeq)
//                                {
//                                    found=(link->source()->sID==otherSourceImage->sID);
//                                    if(found)
//                                    {
//                                        break;
//                                    }
//                                }
//                                if(found)
//                                {
//                                    break;
//                                }
//                            }
//                            
//                            
//                        }
//                        
//                        
//                        FINISH HERE
//                        
//                        //                }
//                    }
//                }
//            }
//            
//
//            
//        }
        
        void print() const
        {
            std::cout<<"face "<<std::flush;
            for(const auto& val : faceIDs)
            {
                std::cout<<val<<" "<<std::endl;
            }
            
            for(const auto& deq : *this)
            {
                for(const auto& link : deq)
                {
                    std::cout<<link->tag()<<std::endl;
                }
                std::cout<<"---------"<<std::endl;
            }
        }
        
    };
}
#endif
