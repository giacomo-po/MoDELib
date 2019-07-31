/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury <ypachaur@purdue.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDislocationSuperLoop_H_
#define model_PeriodicDislocationSuperLoop_H_

#include <string>
#include <TextFileParser.h>
#include <IDreader.h>
#include <MPIcout.h>
#include <PeriodicDislocationLoopPair.h>
#include <utility>

namespace model
{
    template<typename DislocationNetworkType>
    struct PeriodicDislocationSuperLoop
    {
        static constexpr int dim=DislocationNetworkType::dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
        typedef typename TypeTraits<DislocationNetworkType>::NodeType NodeType;
        typedef typename TypeTraits<DislocationNetworkType>::LoopLinkType LoopLinkType;
        //typedef typename TypeTraits<DislocationNetworkType>::NodeType NodeType;
        
        typedef BoundaryLoopLinkSequence<LoopType> BoundaryLoopLinkSequenceType;
        typedef PeriodicDislocationLoopPair<BoundaryLoopLinkSequenceType> PeriodicDislocationLoopPairType;
        
        DislocationNetworkType& DN;
        const PeriodicDislocationLoopPairType& loopPair;
        
        /**********************************************************************/
        PeriodicDislocationSuperLoop(DislocationNetworkType& DN_in,
                                     const PeriodicDislocationLoopPairType& loopPair_in) :
        /* */ DN(DN_in)
        /* */,loopPair(loopPair_in)
        {
            
        }
        
        /**********************************************************************/
        LoopLinkType* outLinkInLoop(const std::shared_ptr<NodeType>& node,
                                          const LoopType* const loop) const
        {
            const auto linksByLoopID(node->linksByLoopID());
            const auto iter(linksByLoopID.find(loop->sID));
            if(iter!=linksByLoopID.end())
            {
                assert(iter->second.size()==2);
                return ((*iter->second.begin())->source()->sID==node->sID? *iter->second.begin() : *iter->second.rbegin());
            }
            else
            {
                return nullptr;
            }
        }
        
        /**********************************************************************/
        std::pair<std::shared_ptr<NodeType>,LoopLinkType *const> isFaultedSuperLoop()
        {

            auto currentLinkSeq(*loopPair.begin());
            auto otherLinkSeq(*loopPair.rbegin());
            auto currentLoop(currentLinkSeq->loop);
            auto otherLoop(otherLinkSeq->loop);
            
            LoopLinkType* currentLoopLink(nullptr);

            for(const auto& loopLink : currentLoop->links())
            {
                if(!loopLink.second->source()->isOnMeshFaces(currentLinkSeq->faceIDs))
                {
                    currentLoopLink=loopLink.second;
                    break;
                }
            }
            
            assert(currentLoopLink && "INTERNAL NODE NOT FOUND");
            
            const auto startSource(currentLoopLink->source());
            while(currentLoopLink->sink()->sID!=startSource->sID)
            {
                
                if(currentLoopLink->sink()->isOnMeshFaces(currentLinkSeq->faceIDs))
                {// sink of currentLoopLink on mesh face
                    const std::shared_ptr<NodeType> sinkImage(currentLoopLink->sink()->sharedImage(currentLinkSeq->faceIDs));
                    
                    const auto outLinkInOtherLoop(outLinkInLoop(sinkImage,otherLoop));
                    if(!outLinkInOtherLoop)
                    {// sinkImage is no in otherLoop, super loop is faulted
                        std::shared_ptr<NodeType> lastImageInOtherLoop(nullptr);
                        //bool found(false);
                        for(size_t k=0;k<currentLoop->links().size();k++)
                        {
                            currentLoopLink=currentLoopLink->prev;
                            if(currentLoopLink->sink()->isOnMeshFaces(currentLinkSeq->faceIDs))
                            {
                                std::shared_ptr<NodeType> testImage(currentLoopLink->sink()->sharedImage(currentLinkSeq->faceIDs));
                                const auto outLink(outLinkInLoop(testImage,otherLoop));
                                if(outLink)
                                {
                                    return std::make_pair(sinkImage,outLink); // returning the link to be expanded, and the new node to be inserted
                                }
                            }
                        }
                        assert(false && "NO NODE IMAGES IN OTHER LOOP");
                    }

                    else
                    {// sinkImage found in otherLoop, continue superLoop on other side
                        //assert(iter->second.size()==2);
                        currentLoopLink=outLinkInOtherLoop;
                        std::swap(currentLinkSeq,otherLinkSeq);
                        currentLoop=currentLinkSeq->loop;
                        otherLoop=otherLinkSeq->loop;
                    }
                    
                }
                else
                {// sink of currentLoopLink not on mesh face, keep mooving along current loop
                    currentLoopLink=currentLoopLink->next;
                }
                
            }
            return std::make_pair(std::shared_ptr<NodeType>(nullptr),nullptr); // super loop not faulted
        }
        
        /**********************************************************************/
        void repair()
        {
            assert(DN.mesh.regions().size()==1);
            const auto& region(*DN.mesh.regions().begin()->second);
 
            switch (loopPair.size())
            {
                case 1:
                {// nucleate the new loop
                    std::cout<<"LoopPair size="<<loopPair.size()<<std::endl;
                    
                    std::vector<std::shared_ptr<NodeType>> imageNodes;
                    const BoundaryLoopLinkSequence<LoopType>& loopLinkSequence(**loopPair.begin());
                    const VectorDim bndNormal(region.outNormal(loopLinkSequence.faceIDs));
                    const VectorDim glideNormal(loopLinkSequence.loop->glidePlane->unitNormal);
                    const VectorDim glideDir(bndNormal-bndNormal.dot(glideNormal)*glideNormal);
                    
                    for(const auto& linkSequence : loopLinkSequence)
                    {
                        const auto startImage(linkSequence.front()->source()->sharedImage(loopLinkSequence.faceIDs));
                        const auto   endImage(linkSequence.back()   ->sink()->sharedImage(loopLinkSequence.faceIDs));
                        
                        const VectorDim centerNodePos(0.5*(startImage->get_P()+endImage->get_P())+glideDir.normalized()*10.0);
                        std::shared_ptr<NodeType> centerNode(new NodeType(&DN,centerNodePos,VectorDim::Zero(),1.0));
                        
                        imageNodes.push_back(startImage);
                        imageNodes.push_back(centerNode);
                        imageNodes.push_back(endImage);
                    }
                    
                    DN.insertLoop(imageNodes,loopLinkSequence,glideNormal,imageNodes[0]->get_P());
                    break;
                }
                    
                case 2:
                {// update loops
                    std::cout<<"LoopPair size="<<loopPair.size()<<std::endl;
                    
                    const auto faultedSuperLoop(isFaultedSuperLoop());
                    if(faultedSuperLoop.first!=nullptr && faultedSuperLoop.second!=nullptr)
                    {// superLoop is faulted, fix it by expanding
                        std::cout<<"faulted SuperLoop. Expanding "<<faultedSuperLoop.second->tag()<<", with Node "<<faultedSuperLoop.first->sID<<std::endl;

                        DN.expandLoopLink(faultedSuperLoop.second,faultedSuperLoop.first);
                        
                        PeriodicDislocationSuperLoop(DN,loopPair).repair(); // call recursively to fix other possible faults
                    }
                    else
                    {
                        std::cout<<"SuperLoop NOT faulted."<<std::endl;

                    }
                    break;
                }
                    
                default:
                {
                    std::cout<<"LoopPair size="<<loopPair.size()<<std::endl;
                    assert(false && "LoopPair must gave size 1 or 2");
                    break;
                }
            }
        }
        
    };
}
#endif
