/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNode_cpp_
#define model_DislocationLoopNode_cpp_

#include <DislocationLoopNode.h>

namespace model
{
    
//    template <int dim, short unsigned int corder>
//    DislocationLoopNode<dim,corder>::DislocationLoopNode(typename DislocationLoopNode<dim,corder>::LoopNetworkType* const net,
//                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::LoopType>& loop,
//                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::NetworkNodeType>& networkNode,
//                                                                           const size_t& edgeID) :
//    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
//    /* init */,_periodicPlanePatch(loop->periodicGlidePlane? loop->periodicGlidePlane->getPatch(this->get_P()) : nullptr)
//    /* init */,periodicPlaneEdge(_periodicPlanePatch? _periodicPlanePatch->edges()[edgeID] : nullptr)
//    {
//
//    }
    
    // template <int dim, short unsigned int corder>
    // DislocationLoopNode<dim,corder>::DislocationLoopNode(typename DislocationLoopNode<dim,corder>::LoopNetworkType* const net,
    //                                                                        const std::shared_ptr<typename DislocationLoopNode<dim,corder>::LoopType>& loop,
    //                                                                        const std::shared_ptr<typename DislocationLoopNode<dim,corder>::NetworkNodeType>& networkNode,
    //                                                                        const VectorDim& P,
    //                                                                        const std::shared_ptr<PeriodicPlanePatch<dim>>& patch_in,
    //                                                                        const std::pair<const std::shared_ptr<PeriodicPlaneEdge<dim>>,const std::shared_ptr<PeriodicPlaneEdge<dim>>>& edge_in) :    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    // /* init */,SplineNodeType(P)
    // /* init */,_periodicPlanePatch(patch_in)
    // /* init */,periodicPlaneEdge(edge_in)
    // {
    //     VerboseDislocationLoopNode(1,"Creating LoopNode "<<this->tag()<<std::endl;);
        
    // }

    template <int dim, short unsigned int corder>
    DislocationLoopNode<dim,corder>::DislocationLoopNode(typename DislocationLoopNode<dim,corder>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::NetworkNodeType>& networkNode,
                                                                           const VectorDim& P,
                                                                           const std::shared_ptr<PeriodicPlanePatch<dim>>& patch_in,
                                                                           const std::pair<const std::shared_ptr<PeriodicPlaneEdge<dim>>,const std::shared_ptr<PeriodicPlaneEdge<dim>>>& edge_in
//                                                                           const std::set<std::shared_ptr<PeriodicPlanePatch<dim>>>& auxPatch_in
                                                                           ) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(P)
    /* init */,_periodicPlanePatch(patch_in)
    /* init */,periodicPlaneEdge(edge_in)
//    /* init */,_auxperiodicPlanePatch(auxPatch_in)
    {
        VerboseDislocationLoopNode(1,"Creating LoopNode "<<this->tag()<<std::endl;);
        if (periodicPlaneEdge.second!=nullptr)
        {
            if (periodicPlaneEdge.first == nullptr)
            {
                assert(false && "Inconsisten definition of periodic plane edges for the dislocation loop node");
            }
        }
        
    }

    template <int dim, short unsigned int corder>
    DislocationLoopNode<dim,corder>::DislocationLoopNode(typename DislocationLoopNode<dim,corder>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder>::NetworkNodeType>& networkNode,
                                                                           const LoopLinkType* const loopLink) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(loopLink->periodicPlanePatch()? networkNode->get_P()-loopLink->periodicPlanePatch()->shift : networkNode->get_P())
    /* init */,_periodicPlanePatch(loopLink->periodicPlanePatch())
    /* init */,periodicPlaneEdge(std::make_pair(nullptr,nullptr))
//    /* init */,_auxperiodicPlanePatch(std::set<std::shared_ptr<PeriodicPlanePatch<dim>>>{nullptr}
    {
        VerboseDislocationLoopNode(1,"Creating LoopNode without PeriodicPlaneEdge "<<this->tag()<<std::endl;);
        
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationLoopNode<dim,corder>::isContractableTo(const LoopNodeType* const other) const
    {
        return LoopNode<LoopNodeType>::isContractableTo(other) && _periodicPlanePatch==other->periodicPlanePatch();
    }

    
    template <int dim, short unsigned int corder>
    std::shared_ptr<DislocationLoopNode<dim,corder>> DislocationLoopNode<dim,corder>::clone(const std::shared_ptr<typename DislocationLoopNode<dim,corder>::LoopType>& otherLoop,
                                                                                                                                const std::shared_ptr<typename DislocationLoopNode<dim,corder>::NetworkNodeType>& otherNetworkNode) const
    {
        return std::shared_ptr<DislocationLoopNode<dim,corder>>(new DislocationLoopNode(this->p_network(),otherLoop,otherNetworkNode,this->get_P(),_periodicPlanePatch,periodicPlaneEdge));
    }
    
    template <int dim, short unsigned int corder>
    void DislocationLoopNode<dim,corder>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopNode=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopNode",true);
    }

    // template <int dim, short unsigned int corder>
    // const DislocationLoopNode<dim, corder> *DislocationLoopNode<dim, corder>::periodicPrev() const
    // {
    //     auto currentPrev(this->prev.first);
    //     while (currentPrev->periodicPlaneEdge)
    //     {
    //         currentPrev = currentPrev->prev.first;
    //     }
    //     return currentPrev;
    // }

    // template <int dim, short unsigned int corder>
    // const DislocationLoopNode<dim, corder> *DislocationLoopNode<dim, corder>::periodicNext() const
    // {
    //     auto currentNext(this->next.first);
    //     while (currentNext->periodicPlaneEdge)
    //     {
    //         currentNext = currentNext->next.first;
    //     }
    //     return currentNext;
    // }

    template <int dim, short unsigned int corder>
    const DislocationLoopNode<dim,corder>* DislocationLoopNode<dim,corder>::periodicPrev() const
    {
        if (this->prev.first)
        {
            auto currentPrev(this->prev.first);
            while (currentPrev->periodicPlaneEdge.first) //First is sufficient
            {
                if (currentPrev == this)   //This if statement can return the current node as periodicPrev. Discuss this with Dr. Po
                {
                    return nullptr;
                }
                currentPrev = currentPrev->prev.first;
            }
            return currentPrev;
        }
        return nullptr;

    }
    
    template <int dim, short unsigned int corder>
    const DislocationLoopNode<dim,corder>* DislocationLoopNode<dim,corder>::periodicNext() const
    {
         auto currentNext(this->next.first);
        if (currentNext)
        {
            while (currentNext->periodicPlaneEdge.first)
            {
                if (currentNext == this)  //This if statement can return the current node as periodicNext. Discuss this with Dr. Po
                {
                    return nullptr;
                }
                currentNext = currentNext->next.first;

            }
            return currentNext;
        }
        return nullptr;

    }
    
    template <int dim, short unsigned int corder>
    std::vector<DislocationLoopNode<dim,corder>*> DislocationLoopNode<dim,corder>::boundaryPrev() const
    {
        assert(!this->periodicPlaneEdge.first);
        std::vector<DislocationLoopNode<dim,corder>*> temp;
        if (this->prev.first)
        {
            auto currentPrev(this->prev.first);
            while (currentPrev->periodicPlaneEdge.first)
            {
                temp.push_back(currentPrev);
                currentPrev = currentPrev->prev.first;
            }
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder>
    std::vector<DislocationLoopNode<dim,corder>*> DislocationLoopNode<dim,corder>::boundaryNext() const
    {
                assert(!this->periodicPlaneEdge.first);
        std::vector<DislocationLoopNode<dim,corder>*> temp;
        if (this->next.first)
        {
            auto currentNext(this->next.first);
            while (currentNext->periodicPlaneEdge.first)
            {
                temp.push_back(currentNext);
                currentNext = currentNext->next.first;
            }
        }

        return temp;
    }

//     /**********************************************************************/
//     template <int dim, short unsigned int corder>
//     void DislocationLoopNode<dim,corder>::addLoopLink(LoopLinkType* const pL)
//     {/*@param[in] pL LoopLink pointer
//       *
//       * This functin overrides LoopNode::addLoopLink
//       */
        
//         VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->sID<<" addLoopLink "<<pL->tag()<<std::endl;);
//         LoopNode<LoopNodeType>::addLoopLink(pL); // forward to base class
//         if(periodicPlanePatch())
//         {
//             this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
//         }
//         else
//         {
//             this->networkNode->addGlidePlane(pL->loop->glidePlane.get());
//         }

// //        updatePeriodicNodes(pL);
//     }
    
    
//     /**********************************************************************/
//     template <int dim, short unsigned int corder>
//     void DislocationLoopNode<dim,corder>::removeLoopLink(LoopLinkType* const pL)
//     {/*@param[in] pL LoopLink pointer
//       * This functin overrides LoopNode::removeLoopLink
//       */
// //LoopNode        VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
//         VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
//         LoopNode<LoopNodeType>::removeLoopLink(pL); // forward to base class
//         VerboseDislocationLoopNode(3,"networkNode->glidePlanes().size()= "<<this->networkNode->glidePlanes().size()<<std::endl;);
//         this->networkNode->confinedObject().clear();
//         //        VerboseDislocationLoopNode(3,"networkNode->glidePlanes().size()= "<<this->networkNode->glidePlanes().size()<<std::endl;);
//         //Add the glidePlanes due to other loopnodes in the networknode //Added by Yash

//         for (const auto& loopNode : this->networkNode->loopNodes())
//         {
//             if (loopNode != this)
//             {
//                 if (loopNode->periodicPlanePatch())
//                 {
//                     this->networkNode->addGlidePlane(loopNode->periodicPlanePatch()->glidePlane.get());
//                 }
//             }
//         }

//         //Giacomo Version
//         if(this->prev.second)
//         {
//             if(periodicPlanePatch())
//             {
//                 VerboseDislocationLoopNode(3,"prev="<<this->prev.second->tag()<<std::endl;);
//                 this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
//             }
//             else
//             {
//                 this->networkNode->addGlidePlane(this->prev.second->loop->glidePlane.get());
//             }
//         }
//         if(this->next.second)
//         {
//             if(periodicPlanePatch())
//             {
//                 VerboseDislocationLoopNode(3,"next="<<this->next.second->tag()<<std::endl;);
//                 this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
//             }
//             else
//             {
//                 this->networkNode->addGlidePlane(this->next.second->loop->glidePlane.get());
//             }
//         }
//         // //Yash Version (Do this for all the loop nodes of this network node)
//         // for (const auto& loopNode : this->networkNode->loopNodes())
//         // {
//         //     if (loopNode->prev.second)
//         //     {
//         //         if (loopNode->periodicPlanePatch())
//         //         {
//         //             VerboseDislocationLoopNode(3, "prev=" << loopNode->prev.second->tag() << std::endl;);
//         //             loopNode->networkNode->addGlidePlane(loopNode->periodicPlanePatch()->glidePlane.get());
//         //         }
//         //         else
//         //         {
//         //             loopNode->networkNode->addGlidePlane(loopNode->prev.second->loop->glidePlane.get());
//         //         }
//         //     }
//         //     if (loopNode->next.second)
//         //     {
//         //         if (loopNode->periodicPlanePatch())
//         //         {
//         //             VerboseDislocationLoopNode(3, "next=" << loopNode->next.second->tag() << std::endl;);
//         //             loopNode->networkNode->addGlidePlane(loopNode->periodicPlanePatch()->glidePlane.get());
//         //         }
//         //         else
//         //         {
//         //             loopNode->networkNode->addGlidePlane(loopNode->next.second->loop->glidePlane.get());
//         //         }
//         //     }
//         // }
        

//         VerboseDislocationLoopNode(3,"DislocationLoopNode "<<this->sID<<" removeLoopLink DONE"<<std::endl;);

// //        for(const auto& loopLink : this->loopLinks())
// //        {
// //            this->networkNode->addGlidePlane(loopLink->loop()->glidePlane.get());
// //        }
        
// //        VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" finished removeLoopLink "<<pL->tag()<<std::endl;);
//     }

    template <int dim, short unsigned int corder>
    void DislocationLoopNode<dim,corder>::set_P(const typename DislocationLoopNode<dim,corder>::VectorLowerDim& newP)
    {
       VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P (lowerDim)"<<std::endl;);
        assert(periodicPlaneEdge.first);

        VectorDim globalPosition(this->loop()->periodicGlidePlane->referencePlane->globalPosition(newP));
        VectorDim localPosition(VectorDim::Zero());

        if (periodicPlaneEdge.second==nullptr)
        {
            localPosition=periodicPlaneEdge.first->meshIntersection->snap(globalPosition + periodicPlaneEdge.first->patch->shift);
            // if (this->networkNode->meshFaces().size()==2)
            // {
            //     //This means that the node was created at a mesh face and not at the edge
            //     //Clear the confinement so that it can move away from the edge
            //     this->networkNode->meshFaces().clear();
            // }
        }
        else
        {
            SegmentSegmentDistance<3> ssd(periodicPlaneEdge.first->meshIntersection->P0,periodicPlaneEdge.first->meshIntersection->P1,
            periodicPlaneEdge.second->meshIntersection->P0,periodicPlaneEdge.second->meshIntersection->P1);
            assert(ssd.dMin<FLT_EPSILON &&  "The two periodic plane edges must intersect");
            assert((periodicPlaneEdge.first->patch->shift-periodicPlaneEdge.second->patch->shift).norm()<FLT_EPSILON && "Two patch shifts must match ");
            localPosition=0.5 * (ssd.x0 + ssd.x1);
        }
        globalPosition = localPosition - periodicPlaneEdge.first->patch->shift;

        // This global local is necessary so as to have the precision for the boundary node and eliminate any floating point error

        SplineNodeType::set_P(globalPosition);
        this->networkNode->set_P(localPosition);
        // this->networkNode->updateConfinement();
    }

//Yash's Version temporary

    // template <int dim, short unsigned int corder>
    // void DislocationLoopNode<dim,corder>::set_P(const typename DislocationLoopNode<dim,corder>::VectorLowerDim& newP)
    // {
    //     VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P (lowerDim)"<<std::endl;);
    //     assert(periodicPlaneEdge);

    //     //If we can get rid of the boundary junctions we do not need the snap operation
    //     const VectorDim globalRVEPosition(this->networkNode->snapToGlidePlanes(this->loop()->periodicGlidePlane->referencePlane->globalPosition(newP)+periodicPlaneEdge->patch->shift));
    //     SplineNodeType::set_P(globalRVEPosition-periodicPlaneEdge->patch->shift);
    //     // std::cout<<"this->networkNode "<<this->networkNode->sID<<this->networkNode->glidePlanes().size()<<std::endl;
    //     this->networkNode->set_P(globalRVEPosition);
    // }

    
    template <int dim, short unsigned int corder>
    void DislocationLoopNode<dim,corder>::set_P(const typename DislocationLoopNode<dim,corder>::VectorDim& newP)
    {
       VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P "<<std::endl;);

        if(this->loop()->glidePlane)
        {
            if (!this->loop()->glidePlane->contains(newP))
            {
                std::cout <<this->tag() << "loopNode is " << this->get_P().transpose() << std::endl;
                std::cout << "Position different is " << (newP - (this->loop()->glidePlane->snapToPlane(newP))).squaredNorm()  << std::endl;
            }
            assert(this->loop()->glidePlane->contains(newP));
        }
        
        if(this->network().simulationParameters.isPeriodicSimulation())
        {
            if(!periodicPlaneEdge.first)
            {// only move loop nodes not on patch boundaries
                VerboseDislocationLoopNode(3,"DislocationLoopNode "<<this->tag()<<" not on periodicPlaneEdge"<<std::endl;);
                VerboseDislocationLoopNode(3,"DislocationLoopNode "<<this->tag()<<" @ "<<this->get_P().transpose()<<std::endl;);
                SplineNodeType::set_P(newP);
                VerboseDislocationLoopNode(3,"DislocationLoopNode now "<<this->tag()<<" @ "<<this->get_P().transpose()<<std::endl;);
                //assert(this->loop()->glidePlane->contains(this->get_P()));
                const auto oldPatch(_periodicPlanePatch);
                const auto pLocalNew(this->loop()->periodicGlidePlane->referencePlane->localPosition(this->get_P()));
                _periodicPlanePatch=this->loop()->periodicGlidePlane->getPatch(this->loop()->periodicGlidePlane->findPatch(pLocalNew,oldPatch->shift));
                VerboseDislocationLoopNode(4,"old patch= "<<oldPatch->shift.transpose()<<std::endl;);
                VerboseDislocationLoopNode(4,"new patch= "<<_periodicPlanePatch->shift.transpose()<<std::endl;);
                // if(oldPatch!=_periodicPlanePatch)
                // {
                //     this->networkNode->confinedObject().clear();
                //     for(auto& neighbor : this->networkNode->neighbors())
                //     {
                //         std::get<1>(neighbor.second)->confinedObject().clear();
                //     }
                // }
                //Added for the nodes which are not classified as boundary by present at the boundary (Confined object classifies those nodes as onExternalBoundary())
                // if (this->networkNode->isOnExternalBoundary())
                // {
                //     //The node is an internal node present at the boundary(Clear the mesh faces so that it can move within the patches)
                //     this->networkNode->meshFaces().clear();
                //     for(auto& neighbor : this->networkNode->neighbors())
                //     {
                //         std::get<1>(neighbor.second)->meshFaces().clear();
                //     }
                // }
                this->networkNode->set_P(_periodicPlanePatch->glidePlane->snapToPlane(this->get_P()+_periodicPlanePatch->shift));
                // this->networkNode->set_P(this->get_P()+_periodicPlanePatch->shift);
                // if(oldPatch!=_periodicPlanePatch)
                // {
                //     //Add the glide plane for the networkNode for the new patch
                //     this->networkNode->addGlidePlane(_periodicPlanePatch->glidePlane.get());
                // }
                // if(this->prev.second->hasNetworkLink() && this->prev.second->loop->glidePlane)
                // {
                //     if(this->prev.second->source->periodicPlanePatch()==this->prev.second->sink->periodicPlanePatch())
                //     {
                //         this->prev.second->networkLink()->addGlidePlane(this->prev.second->source->periodicPlanePatch()->glidePlane.get());
                //     }
                // }
                
                // if(this->next.second->hasNetworkLink() && this->next.second->loop->glidePlane)
                // {
                //     if(this->next.second->source->periodicPlanePatch()==this->next.second->sink->periodicPlanePatch())
                //     {
                //         this->next.second->networkLink()->addGlidePlane(this->next.second->source->periodicPlanePatch()->glidePlane.get());
                //     }
                // }
                
                
                for(auto& bndNode : boundaryPrev())
                {
                    VerboseDislocationLoopNode(3,"Setting P for  PrevBndNode "<<bndNode->tag()<<std::endl;);
                    const auto pPrev(bndNode->periodicPrev());
                    const auto pNext(bndNode->periodicNext());
                    if (pPrev && pNext)
                    {
                        assert(pNext == this);
                        VerboseDislocationLoopNode(4, "pPrev= " << pPrev->tag() << " @ " << pPrev->get_P().transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "pNext= " << pNext->tag() << " @ " << pNext->get_P().transpose() << std::endl;);

                        const auto pPrevLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                        VerboseDislocationLoopNode(4, "pPrevLocal= " << pPrevLocal.transpose() << std::endl;);
                        const auto pNextLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                        VerboseDislocationLoopNode(4, "pNextLocal= " << pNextLocal.transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "periodicPlaneEdge.first->source= " << bndNode->periodicPlaneEdge.first->source->transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "bndNode->periodicPlaneEdge.first->sink= " << bndNode->periodicPlaneEdge.first->sink->transpose() << std::endl;);
                        if (periodicPlaneEdge.second)
                        {
                            VerboseDislocationLoopNode(4, "periodicPlaneEdge.second->source= " << bndNode->periodicPlaneEdge.second->source->transpose() << std::endl;);
                            VerboseDislocationLoopNode(4, "bndNode->periodicPlaneEdge.second->sink= " << bndNode->periodicPlaneEdge.second->sink->transpose() << std::endl;);

                        }
                        if (periodicPlaneEdge.second)
                        {
                            SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                            SegmentSegmentDistance<dim - 1> ssd2(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.second->source, *bndNode->periodicPlaneEdge.second->sink);
                            if (ssd1.dMin < FLT_EPSILON)
                            {
                                assert(ssd2.dMin<FLT_EPSILON && "Other intersection must exist");
                                const VectorLowerDim ssd1Pos(0.5 * (ssd1.x0 + ssd1.x1));
                                const VectorLowerDim ssd2Pos(0.5 * (ssd2.x0 + ssd2.x1));
                                assert((ssd1Pos-ssd2Pos).norm()<FLT_EPSILON && "The two intersection point must be same ");
                                VerboseDislocationLoopNode(3, "dMin= " << ssd1.dMin << std::endl;);
                                bndNode->set_P(ssd1Pos);
                            }
                            else
                            {
                                bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                                this->network().danglingBoundaryLoopNodes.insert(bndNode);
                            }
                        }
                        else
                        {
                            SegmentSegmentDistance<dim - 1> ssd(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                            if (ssd.dMin < FLT_EPSILON)
                            {
                                VerboseDislocationLoopNode(3, "dMin= " << ssd.dMin << std::endl;);
                                bndNode->set_P(VectorLowerDim(0.5 * (ssd.x0 + ssd.x1)));
                            }
                            else
                            {
                                bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                                this->network().danglingBoundaryLoopNodes.insert(bndNode);
                            }
                        }

                    }
                }
                for(auto& bndNode : boundaryNext())
                {
                                        VerboseDislocationLoopNode(3,"Setting P for  NextBndNode "<<bndNode->tag()<<std::endl;);
                    const auto pPrev(bndNode->periodicPrev());
                    const auto pNext(bndNode->periodicNext());
                    if (pPrev && pNext)
                    {
                        assert(pPrev == this);
                        VerboseDislocationLoopNode(4, "pPrev= " << pPrev->tag() << " @ " << pPrev->get_P().transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "pNext= " << pNext->tag() << " @ " << pNext->get_P().transpose() << std::endl;);

                        const auto pPrevLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                        VerboseDislocationLoopNode(4, "pPrevLocal= " << pPrevLocal.transpose() << std::endl;);
                        const auto pNextLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                        VerboseDislocationLoopNode(4, "pNextLocal= " << pNextLocal.transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "periodicPlaneEdge->source= " << bndNode->periodicPlaneEdge.first->source->transpose() << std::endl;);
                        VerboseDislocationLoopNode(4, "bndNode->periodicPlaneEdge->sink= " << bndNode->periodicPlaneEdge.first->sink->transpose() << std::endl;);
                        if (periodicPlaneEdge.second)
                        {
                            VerboseDislocationLoopNode(4, "periodicPlaneEdge.second->source= " << bndNode->periodicPlaneEdge.second->source->transpose() << std::endl;);
                            VerboseDislocationLoopNode(4, "bndNode->periodicPlaneEdge.second->sink= " << bndNode->periodicPlaneEdge.second->sink->transpose() << std::endl;);
                        }
                        if (periodicPlaneEdge.second)
                        {
                            SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                            SegmentSegmentDistance<dim - 1> ssd2(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.second->source, *bndNode->periodicPlaneEdge.second->sink);
                            if (ssd1.dMin < FLT_EPSILON)
                            {
                                assert(ssd2.dMin < FLT_EPSILON && "Other intersection must exist");
                                const VectorLowerDim ssd1Pos(0.5 * (ssd1.x0 + ssd1.x1));
                                const VectorLowerDim ssd2Pos(0.5 * (ssd2.x0 + ssd2.x1));
                                assert((ssd1Pos - ssd2Pos).norm() < FLT_EPSILON && "The two intersection point must be same ");
                                VerboseDislocationLoopNode(3, "dMin= " << ssd1.dMin << std::endl;);
                                bndNode->set_P(ssd1Pos);
                            }
                            else
                            {
                                bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                                this->network().danglingBoundaryLoopNodes.insert(bndNode);
                            }
                        }
                        else
                        {
                            SegmentSegmentDistance<dim - 1> ssd(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                            if (ssd.dMin < FLT_EPSILON)
                            {
                                VerboseDislocationLoopNode(3, "dMin= " << ssd.dMin << std::endl;);
                                bndNode->set_P(VectorLowerDim(0.5 * (ssd.x0 + ssd.x1)));
                            }
                            else
                            {
                                bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                                this->network().danglingBoundaryLoopNodes.insert(bndNode);
                            }
                        }
                    }
                    
                }
            }
        }
        else
        {
            SplineNodeType::set_P(newP);
//            assert(false && "FINISH HERE");
        }
    }
    
    template <int dim, short unsigned int corder>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopNode<dim,corder>::periodicPlanePatch() const
    {
        return _periodicPlanePatch;
    }

    // template <int dim, short unsigned int corder>
    // void DislocationLoopNode<dim,corder>::updateConfinedGeometry()
    // {
    //     //This will update the confinement of the loop nodes
    //     // for (const auto &loopNode : this->networkNode->loopNodes())
    //     // {
    //     //     if (loopNode->periodicPlanePatch())
    //     //     {
    //     //         this->networkNode->addGlidePlane(loopNode->periodicPlanePatch()->glidePlane.get());
    //     //     }
    //     // }
    //     if (this->periodicPlanePatch())
    //     {
    //         this->networkNode->addGlidePlane(this->periodicPlanePatch()->glidePlane.get());
    //     }
    // }

    template <int dim, short unsigned int corder>
    void DislocationLoopNode<dim,corder>::updateGeometry()
    {

    }

    template <int dim, short unsigned int corder>
    bool DislocationLoopNode<dim,corder>::isGeometricallyRemovable(const double &Lmin,const double &Lmax,const double& relAreaTh)
    {
        const auto pPrev(periodicPrev());
        const auto pNext(periodicNext());
        if((pNext->get_P()-pPrev->get_P()).norm()<0.95*Lmax)
        {
            const double deltaArea(0.5 * (this->get_P() - pPrev->get_P()).cross(pNext->get_P() - this->get_P()).norm());

            // VectorDim vectorArea(VectorDim::Zero());
            // if (this->loop()->loopLinks().size())
            // {
            //     const VectorDim P0((*this->loop()->loopLinks().begin())->source->get_P());
            //     for (const auto &loopLink : this->loop()->loopLinks())
            //     {
            //         vectorArea += 0.5 * (loopLink->source->get_P() - P0).cross(loopLink->sink->get_P() - loopLink->source->get_P());
            //     }
            // }

            // const double loopArea(vectorArea.norm());
            const double loopArea(this->loop()->slippedArea());
            const double loopAreaRate(this->loop()->slippedAreaRate());
            // std::cout<<" For loop node "<<this->tag()<<" with "<<this->loop()->loopLinks().size()<<" loopLinks. loop area is "<<std::setprecision(15)<<loopArea<<std::endl;
            if (loopArea > std::pow(Lmin, 2))
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " loopArea= " << loopArea << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " deltaArea= " << deltaArea << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " relAreaTh= " << relAreaTh << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " removable? " << (deltaArea / loopArea < relAreaTh) << std::endl;);

                return deltaArea / loopArea < relAreaTh;
            }
            else if (loopArea > FLT_EPSILON && loopAreaRate<-FLT_EPSILON)
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " removable (small loop)" << std::endl;);
                return true;
            }
            else
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (zero area loop)" << std::endl;);
                return false;
            }
        }
        else
        {
            VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (short segment)" << std::endl;);
            return false;
        }
    }
    
//Updated Version for allowing for the removal of twinned loopNode only (This version has worked very well (Reached 1% strain with both FR source and circular loops))
//Boundary node accumulation for annihilation at the boundary    
    // template <int dim, short unsigned int corder>
    // std::pair<bool,size_t> DislocationLoopNode<dim, corder>::isRemovable(const double &Lmin, const double &relAreaTh)
    // {
    //     //The pair will be populated as follows:
    //     //If the node is removable and it does not have a twin the pair will contain the same node
    //     //Otherwise it will contrain the twin (in the removal then we have to remove both the current node and twin)

    //     VerboseDislocationLoopNode(2, " Checking if DislocationLoopNode " << this->tag() << " is removable " << std::endl;);
    //     if (periodicPlaneEdge)
    //     { // a boundary node
    //         VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " onPlaneEdge" << std::endl;);
    //         // const auto pPrev(periodicPrev());
    //         // const auto pNext(periodicNext());
    //         // if (pPrev && pNext)
    //         // {
    //         //     VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (both periodicPrev and periodicNext exist)" << std::endl;);
    //         //     return std::make_pair(false,0);
    //         // }
    //         // else
    //         // {
    //         //     VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << "  removable (periodicPrev and periodicNext do not exist)" << std::endl;);
    //         //     return std::make_pair(true,this->sID);
    //         // }
    //             VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (on Edge)" << std::endl;);
    //             return std::make_pair(false,0);

    //     }
    //     else
    //     { // not a boundary node
    //         VerboseDislocationLoopNode(2, " Dislocation LoopNode Not on Boundary " << std::endl;);

    //         const auto pPrev(periodicPrev());
    //         const auto pNext(periodicNext());
    //         // if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax)
    //         // {
    //         if (pPrev && pNext)
    //         {
    //             if (this->networkNode->loopNodes().size() == 1)
    //             {
    //                 if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax &&
    //                     isGeometricallyRemovable(Lmin, relAreaTh))
    //                 {
    //                     return std::make_pair(true,this->sID);
    //                 }
    //                 else
    //                 {
    //                     return std::make_pair(false,0);
    //                 }
    //             }
    //             else if (this->networkNode->loopNodes().size() >= 2)
    //             {
    //                 const auto prevTwin(this->prev.second->twinnedLink());
    //                 const auto nextTwin(this->next.second->twinnedLink());

    //                 if (prevTwin && nextTwin)
    //                 {

    //                     if (prevTwin->loop == nextTwin->loop && (prevTwin->periodicPlanePatch() == nextTwin->periodicPlanePatch()))
    //                    {
    //                         //Check based on neighbors of the networkNode
    //                         //Here the other node will be inserted
    //                         size_t prevtwinID((prevTwin->source->networkNode==this->networkNode)? prevTwin->source->sID : prevTwin->sink->sID);
    //                         size_t nexttwinID((nextTwin->source->networkNode==this->networkNode)? nextTwin->source->sID : nextTwin->sink->sID);
    //                         // if (prevtwinID!=nexttwinID)
    //                         // {
    //                         //     // prevTwin->loop->printLoop();
    //                         //     std::cout<<"Prev twin->tag() "<<prevTwin->tag()<<std::endl;
    //                         //     std::cout<<"Next twin->tag() "<<nextTwin->tag()<<std::endl;

    //                         //     std::cout<<" Current node id "<<this->tag()<<std::endl;
    //                         //     std::cout<<" PrevTwin tag "<<((prevTwin->source->networkNode==this->networkNode)? prevTwin->source->tag() : prevTwin->sink->tag())<<std::endl;
    //                         //     std::cout<<" NextTwin tag "<<((nextTwin->source->networkNode==this->networkNode)? nextTwin->source->tag() : nextTwin->sink->tag())<<std::endl;

    //                         //     std::cout<<" PrevTwin periodicPrev "<<((prevTwin->source->networkNode==this->networkNode)? prevTwin->source->periodicPrev()->tag() : prevTwin->sink->periodicPrev()->tag())<<std::endl;
    //                         //     std::cout<<" NextTwin periodicNext "<<((nextTwin->source->networkNode==this->networkNode)? nextTwin->source->periodicNext()->tag() : nextTwin->sink->periodicNext()->tag())<<std::endl;

    //                         // }
    //                         // assert((prevtwinID==nexttwinID) && "LoopNode mismatch in the removal of the nodes 0B nodes");
    //                         // return std::make_pair(true,prevtwinID);

    //                         if (prevtwinID==nexttwinID)
    //                         {
    //                             return std::make_pair(true,prevtwinID);
    //                         }
    //                         else
    //                         {
    //                             // prevTwin->loop->printLoop();
    //                             std::cout << "Prev twin->tag() " << prevTwin->tag() << std::endl;
    //                             std::cout << "Next twin->tag() " << nextTwin->tag() << std::endl;

    //                             std::cout << " Current node id " << this->tag() << std::endl;
    //                             std::cout << " PrevTwin tag " << ((prevTwin->source->networkNode == this->networkNode) ? prevTwin->source->tag() : prevTwin->sink->tag()) << std::endl;
    //                             std::cout << " NextTwin tag " << ((nextTwin->source->networkNode == this->networkNode) ? nextTwin->source->tag() : nextTwin->sink->tag()) << std::endl;

    //                             std::cout << " PrevTwin periodicPrev " << ((prevTwin->source->networkNode == this->networkNode) ? prevTwin->source->periodicPrev()->tag() : prevTwin->sink->periodicPrev()->tag()) << std::endl;
    //                             std::cout << " NextTwin periodicNext " << ((nextTwin->source->networkNode == this->networkNode) ? nextTwin->source->periodicNext()->tag() : nextTwin->sink->periodicNext()->tag()) << std::endl;

    //                             assert(false && "LoopNode mismatch in the removal of the nodes 0B nodes");
    //                             return std::make_pair(false,0);
    //                         }
                            
    //                     }
    //                     // else
    //                     // {
    //                         //This is not valid for the case of junctions with different loops
    //                         //In one case it is becoming removable but for the other node it is not removable
    //                     //     if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax)
    //                     //     {
    //                     //         //Geometrically Removable condition should be checked for all the loop nodes corresponding to this network node
    //                     //         bool tempRemovable(true);
    //                     //         for (const auto &loopN : this->networkNode->loopNodes())
    //                     //         {
    //                     //             const auto LoopNpPrev(periodicPrev());
    //                     //             const auto LoopNpNext(periodicNext());

    //                     //             if (LoopNpPrev && LoopNpNext)
    //                     //             {
    //                     //                 tempRemovable *= (loopN->isGeometricallyRemovable(Lmin, FLT_EPSILON) &&
    //                     //                                   (LoopNpPrev->get_P() - LoopNpNext->get_P()).norm() < this->network().networkRemesher.Lmax);
    //                     //             }
    //                     //         }
    //                     //         return std::make_pair(tempRemovable,this->sID);
    //                     //     }
    //                     //     else
    //                     //     {
    //                     //         return std::make_pair(false,0);
    //                     //     }
    //                     // }
    //                     else
    //                     {
    //                         return std::make_pair(false, 0);
    //                     }
    //                 }
    //                 else
    //                 {
    //                     return std::make_pair(false,0);
    //                 }
    //             }
    //             else
    //             {
    //                 return std::make_pair(false,0);
    //             }
    //         }
    //         else
    //         {
    //             return std::make_pair(false,0);
    //         }

    //         // }
    //         // else
    //         // {
    //         //     return false;
    //         // }
    //     }
    // }

    //Version which can classiy the boundary nodes as well...
    template <int dim, short unsigned int corder>
    std::pair<bool,size_t> DislocationLoopNode<dim, corder>::isRemovable(const double &Lmin,const double &Lmax, const double &relAreaTh)
    {
       //The pair will be populated as follows:
        //If the node is removable and it does not have a twin the pair will contain the same node
        //Otherwise it will contrain the twin (in the removal then we have to remove both the current node and twin)

        VerboseDislocationLoopNode(2, " Checking if DislocationLoopNode " << this->tag() << " is removable " << std::endl;);
        if (periodicPlaneEdge.first)
        { // a boundary node
            VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " onPlaneEdge" << std::endl;);
            const auto pPrev(periodicPrev());
            const auto pNext(periodicNext());
            if (pPrev && pNext)
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (both periodicPrev and periodicNext exist)" << std::endl;);
                return std::make_pair(false,0);
            }
            else
            {
                const double loopArea(this->loop()->slippedArea());
                if (loopArea < std::pow(Lmin, 2))
                {
                    VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << "  removable (on Edge via Geometric condition)" << std::endl;);
                    return std::make_pair(true,this->sID);
                }
                else
                {
                    VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (on Edge)" << std::endl;);
                    return std::make_pair(false, 0);
                }
                
            }

        }
        else
        { // not a boundary node
            VerboseDislocationLoopNode(2, " Dislocation LoopNode Not on Boundary " << std::endl;);

            const auto pPrev(periodicPrev());
            const auto pNext(periodicNext());
            // if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax)
            // {
            if (pPrev && pNext)
            {
                VerboseDislocationLoopNode(2, "  Has pPrev and pNext" << std::endl;);

                if (this->networkNode->loopNodes().size() == 1)
                {
                    VerboseDislocationLoopNode(2, "  Has 1 loop nodes" << std::endl;);
                    if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax
//                    if (((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax || pPrev->periodicPlanePatch()!=pNext->periodicPlanePatch())
                        && isGeometricallyRemovable(Lmin, Lmax,relAreaTh))
                    {
                        return std::make_pair(true,this->sID);
                    }
                    else if (this->networkNode->isOnExternalBoundary())
                    {
                        //This means that an internal node is exactly on the the boundary based on the geometrically removable condition, the node can be removed
                        if (isGeometricallyRemovable(Lmin, Lmax, relAreaTh))
                        {
                            return std::make_pair(true,this->sID);
                        }
                        else
                        {
                            return std::make_pair(false,0);
                        }
                        
                    }
                    else
                    {
                        return std::make_pair(false,0);
                    }
                }
                else if (this->networkNode->loopNodes().size() >= 2)
                {
                    VerboseDislocationLoopNode(2, "  Has 2 loop nodes" << std::endl;);
                    const auto prevTwin(this->prev.second->twinnedLink());
                    const auto nextTwin(this->next.second->twinnedLink());

                    if (prevTwin && nextTwin)
                    {
                        VerboseDislocationLoopNode(2, "  Has prevTwin and nextTwin" << std::endl;);

                        if (prevTwin->loop == nextTwin->loop && (prevTwin->periodicPlanePatch() == nextTwin->periodicPlanePatch()))
                       {
                           VerboseDislocationLoopNode(2, "  Same patches" << std::endl;);
                            //Check based on neighbors of the networkNode
                            //Here the other node will be inserted
                            size_t prevtwinID((prevTwin->source->networkNode==this->networkNode)? prevTwin->source->sID : prevTwin->sink->sID);
                            size_t nexttwinID((nextTwin->source->networkNode==this->networkNode)? nextTwin->source->sID : nextTwin->sink->sID);
                            

                            if (prevtwinID==nexttwinID)
                            {
                                VerboseDislocationLoopNode(2, "  Same IDs" << std::endl;);
                                return std::make_pair(true,prevtwinID);
                            }
                            else
                            {
                                // prevTwin->loop->printLoop();
                                std::cout << "Prev twin->tag() " << prevTwin->tag() << std::endl;
                                std::cout << "Next twin->tag() " << nextTwin->tag() << std::endl;

                                std::cout << " Current node id " << this->tag() << std::endl;
                                std::cout << " PrevTwin tag " << ((prevTwin->source->networkNode == this->networkNode) ? prevTwin->source->tag() : prevTwin->sink->tag()) << std::endl;
                                std::cout << " NextTwin tag " << ((nextTwin->source->networkNode == this->networkNode) ? nextTwin->source->tag() : nextTwin->sink->tag()) << std::endl;

                                std::cout << " PrevTwin periodicPrev " << ((prevTwin->source->networkNode == this->networkNode) ? prevTwin->source->periodicPrev()->tag() : prevTwin->sink->periodicPrev()->tag()) << std::endl;
                                std::cout << " NextTwin periodicNext " << ((nextTwin->source->networkNode == this->networkNode) ? nextTwin->source->periodicNext()->tag() : nextTwin->sink->periodicNext()->tag()) << std::endl;

                                throw std::runtime_error("LoopNode mismatch in the removal of the nodes 0B nodes");
                                
//                                assert(false && "LoopNode mismatch in the removal of the nodes 0B nodes");
                                return std::make_pair(false,0);
                            }
                            
                        }
                        else
                        {
                            VerboseDislocationLoopNode(2, "  Different patches" << std::endl;);
                            return std::make_pair(false, 0);
                        }
                    }
                    else
                    {
                        VerboseDislocationLoopNode(2, "  No prevTwin or nextTwin" << std::endl;);
                        return std::make_pair(false,0);
                    }
                }
                else
                {
                    VerboseDislocationLoopNode(2, "this->networkNode->loopNodes().size()="<<this->networkNode->loopNodes().size()<< std::endl;);
                    return std::make_pair(false,0);
                }
            }
            else
            {
                VerboseDislocationLoopNode(2, "  No pPrev or pNext" << std::endl;);
                return std::make_pair(false,0);
            }

           
        }
    }

    template <int dim, short unsigned int corder>
    bool DislocationLoopNode<dim,corder>::isMovableTo(const VectorDim& X) const
    {
       VerboseDislocationLoopNode(4,"checking if DislocationLoopNode "<<this->sID<< " isMovable:"<<std::endl;);

        bool isMovable=true;
        
        if(periodicPlaneEdge.first)
        {
            isMovable=(isMovable && ((X-this->get_P()).squaredNorm()<FLT_EPSILON));
        }
        
        if(isMovable && this->loop()->glidePlane)
        {
            isMovable=(isMovable && this->loop()->glidePlane->contains(X));
        }
        else
        {
            isMovable=false;
        }
        
        return isMovable;

        
//        if(isMovable && this->loop()->isSessile)
//        {
//
//        }
//
//        if(this->network().simulationParameters.isPeriodicSimulation() && this->isBoundaryNode())
//        {// cannot move boundary nodes under periodic simulations
//            return (X-this->get_P()).squaredNorm()<FLT_EPSILON;
//        }
//        else
//        {
//            bool isMovable=true;
//
//            VerboseDislocationNode(4,"checking if PlanarDislocationNode "<<this->sID<< " isMovable:"<<std::endl;);
//
//            for(const auto& gp : this->glidePlanes())
//            {// X must be contained by all glidePlanes
//                isMovable*=gp->contains(X);
//            }
//            VerboseDislocationNode(4,"  meshPlanes contains X? "<<isMovable<<std::endl;);
//
//            if(isMovable)
//            {
//                for(const auto& pair : this->neighbors())
//                {
//                    if(std::get<1>(pair.second)->isSessile())
//                    {// sessile segments cannot change direction if this node is moved
//                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
//                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
//                        if(currentNorm>FLT_EPSILON && newNorm>FLT_EPSILON)
//                        {
//                            const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
//                            VerboseDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
//                            isMovable=(isMovable&&sessileNeighborMovable);
//                            //                                isMovable*=sessileNeighborMovable;
//                            if(!isMovable)
//                            {
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//
//            if(isMovable && this->isOnBoundary())
//            {
//                isMovable*=this->boundingBoxSegments().contains(X);
//            }
//        }
        
        
    }
    
    template <int dim, short unsigned int corder>
    int DislocationLoopNode<dim,corder>::verboseDislocationLoopNode=0;

    template class DislocationLoopNode<3,0>;
}
#endif
