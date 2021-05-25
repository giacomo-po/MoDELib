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
#include <SegmentSegmentDistance.h>

namespace model
{
    
//    template <int dim, short unsigned int corder, typename InterpolationType>
//    DislocationLoopNode<dim,corder,InterpolationType>::DislocationLoopNode(typename DislocationLoopNode<dim,corder,InterpolationType>::LoopNetworkType* const net,
//                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& loop,
//                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& networkNode,
//                                                                           const size_t& edgeID) :
//    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
//    /* init */,_periodicPlanePatch(loop->periodicGlidePlane? loop->periodicGlidePlane->getPatch(this->get_P()) : nullptr)
//    /* init */,periodicPlaneEdge(_periodicPlanePatch? _periodicPlanePatch->edges()[edgeID] : nullptr)
//    {
//
//    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationLoopNode<dim,corder,InterpolationType>::DislocationLoopNode(typename DislocationLoopNode<dim,corder,InterpolationType>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& networkNode,
                                                                           const VectorDim& P,
                                                                           const std::shared_ptr<PeriodicPlanePatch<dim>>& patch_in,
                                                                           const std::shared_ptr<PeriodicPlaneEdge<dim>>& edge_in) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(P)
    /* init */,_periodicPlanePatch(patch_in)
    /* init */,periodicPlaneEdge(edge_in)
    {
        
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationLoopNode<dim,corder,InterpolationType>::DislocationLoopNode(typename DislocationLoopNode<dim,corder,InterpolationType>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& networkNode,
                                                                           const LoopLinkType* const loopLink) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(loopLink->periodicPlanePatch()? networkNode->get_P()-loopLink->periodicPlanePatch()->shift : networkNode->get_P())
    /* init */,_periodicPlanePatch(loopLink->periodicPlanePatch())
    /* init */,periodicPlaneEdge(nullptr)
    {
        
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationLoopNode<dim,corder,InterpolationType>::isContractableTo(const LoopNodeType* const other) const
    {
        return LoopNode<LoopNodeType>::isContractableTo(other) && _periodicPlanePatch==other->periodicPlanePatch();
    }

    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<DislocationLoopNode<dim,corder,InterpolationType>> DislocationLoopNode<dim,corder,InterpolationType>::clone(const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& otherLoop,
                                                                                                                                const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& otherNetworkNode) const
    {
        return std::shared_ptr<DislocationLoopNode<dim,corder,InterpolationType>>(new DislocationLoopNode(this->p_network(),otherLoop,otherNetworkNode,this->get_P(),_periodicPlanePatch,periodicPlaneEdge));
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopNode=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopNode",true);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const DislocationLoopNode<dim,corder,InterpolationType>* DislocationLoopNode<dim,corder,InterpolationType>::periodicPrev() const
    {
        auto currentPrev(this->prev.first);
        while(currentPrev->periodicPlaneEdge)
        {
            currentPrev=currentPrev->prev.first;
        }
        return currentPrev;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const DislocationLoopNode<dim,corder,InterpolationType>* DislocationLoopNode<dim,corder,InterpolationType>::periodicNext() const
    {
        auto currentNext(this->next.first);
        while(currentNext->periodicPlaneEdge)
        {
            currentNext=currentNext->next.first;
        }
        return currentNext;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::vector<DislocationLoopNode<dim,corder,InterpolationType>*> DislocationLoopNode<dim,corder,InterpolationType>::boundaryPrev() const
    {
        assert(!this->periodicPlaneEdge);
        std::vector<DislocationLoopNode<dim,corder,InterpolationType>*> temp;
        auto currentPrev(this->prev.first);
        while(currentPrev->periodicPlaneEdge)
        {
            temp.push_back(currentPrev);
            currentPrev=currentPrev->prev.first;
        }
        return temp;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::vector<DislocationLoopNode<dim,corder,InterpolationType>*> DislocationLoopNode<dim,corder,InterpolationType>::boundaryNext() const
    {
                assert(!this->periodicPlaneEdge);
        std::vector<DislocationLoopNode<dim,corder,InterpolationType>*> temp;
        auto currentNext(this->next.first);
        while(currentNext->periodicPlaneEdge)
        {
            temp.push_back(currentNext);
            currentNext=currentNext->next.first;
        }
        return temp;
    }

    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::addLoopLink(LoopLinkType* const pL)
    {/*@param[in] pL LoopLink pointer
      *
      * This functin overrides LoopNode::addLoopLink
      */
        
        VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->sID<<" addLoopLink "<<pL->tag()<<std::endl;);
        LoopNode<LoopNodeType>::addLoopLink(pL); // forward to base class
        if(periodicPlanePatch())
        {
            this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
        }
        else
        {
            this->networkNode->addGlidePlane(pL->loop->glidePlane.get());
        }

//        updatePeriodicNodes(pL);
    }
    
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::removeLoopLink(LoopLinkType* const pL)
    {/*@param[in] pL LoopLink pointer
      * This functin overrides LoopNode::removeLoopLink
      */
//LoopNode        VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
        VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
        LoopNode<LoopNodeType>::removeLoopLink(pL); // forward to base class
        VerboseDislocationLoopNode(3,"networkNode->glidePlanes().size()= "<<this->networkNode->glidePlanes().size()<<std::endl;);
        this->networkNode->confinedObject().clear();
//        VerboseDislocationLoopNode(3,"networkNode->glidePlanes().size()= "<<this->networkNode->glidePlanes().size()<<std::endl;);
        if(this->prev.second)
        {
            if(periodicPlanePatch())
            {
                VerboseDislocationLoopNode(3,"prev="<<this->prev.second->tag()<<std::endl;);
                this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
            }
            else
            {
                this->networkNode->addGlidePlane(this->prev.second->loop->glidePlane.get());
            }
        }
        if(this->next.second)
        {
            if(periodicPlanePatch())
            {
                VerboseDislocationLoopNode(3,"next="<<this->next.second->tag()<<std::endl;);
                this->networkNode->addGlidePlane(periodicPlanePatch()->glidePlane.get());
            }
            else
            {
                this->networkNode->addGlidePlane(this->next.second->loop->glidePlane.get());
            }
        }
        VerboseDislocationLoopNode(3,"DislocationLoopNode "<<this->sID<<" removeLoopLink DONE"<<std::endl;);

//        for(const auto& loopLink : this->loopLinks())
//        {
//            this->networkNode->addGlidePlane(loopLink->loop()->glidePlane.get());
//        }
        
//        VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" finished removeLoopLink "<<pL->tag()<<std::endl;);
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::set_P(const typename DislocationLoopNode<dim,corder,InterpolationType>::VectorLowerDim& newP)
    {
        VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P (lowerDim)"<<std::endl;);
        assert(periodicPlaneEdge);
        SplineNodeType::set_P(this->loop()->periodicGlidePlane->referencePlane->globalPosition(newP));
        this->networkNode->set_P(this->get_P()+periodicPlaneEdge->patch->shift);
    }

    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::set_P(const typename DislocationLoopNode<dim,corder,InterpolationType>::VectorDim& newP)
    {
        VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P "<<std::endl;);

        if(this->loop()->glidePlane)
        {
            assert(this->loop()->glidePlane->contains(newP));
        }
        
        if(this->network().simulationParameters.isPeriodicSimulation())
        {
            if(!periodicPlaneEdge)
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
                if(oldPatch!=_periodicPlanePatch)
                {
                    this->networkNode->confinedObject().clear();
                    for(auto& neighbor : this->networkNode->neighbors())
                    {
                        std::get<1>(neighbor.second)->confinedObject().clear();
                    }
                }
                this->networkNode->set_P(this->get_P()+_periodicPlanePatch->shift);
                
                if(this->prev.second->hasNetworkLink() && this->prev.second->loop->glidePlane)
                {
                    if(this->prev.second->source->periodicPlanePatch()==this->prev.second->sink->periodicPlanePatch())
                    {
                        this->prev.second->networkLink()->addGlidePlane(this->prev.second->source->periodicPlanePatch()->glidePlane.get());
                    }
                }
                
                if(this->next.second->hasNetworkLink() && this->next.second->loop->glidePlane)
                {
                    if(this->next.second->source->periodicPlanePatch()==this->next.second->sink->periodicPlanePatch())
                    {
                        this->next.second->networkLink()->addGlidePlane(this->next.second->source->periodicPlanePatch()->glidePlane.get());
                    }
                }
                
                
                for(auto& bndNode : boundaryPrev())
                {
                    VerboseDislocationLoopNode(3,"Setting P for  PrevBndNode "<<bndNode->tag()<<std::endl;);
                    const auto pPrev(bndNode->periodicPrev());
                    const auto pNext(bndNode->periodicNext());
                    assert(pNext==this);
                    VerboseDislocationLoopNode(4,"pPrev= "<<pPrev->tag()<<" @ "<<pPrev->get_P().transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"pNext= "<<pNext->tag()<<" @ "<<pNext->get_P().transpose()<<std::endl;);

                    const auto pPrevLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                    VerboseDislocationLoopNode(4,"pPrevLocal= "<<pPrevLocal.transpose()<<std::endl;);
                    const auto pNextLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                    VerboseDislocationLoopNode(4,"pNextLocal= "<<pNextLocal.transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"periodicPlaneEdge->source= "<<bndNode->periodicPlaneEdge->source->transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"bndNode->periodicPlaneEdge->sink= "<<bndNode->periodicPlaneEdge->sink->transpose()<<std::endl;);

//                    SegmentSegmentDistance<dim> ssd3(pPrev->get_P(),pNext->get_P(),bndNode->periodicPlaneEdge->meshIntersection->P0,bndNode->periodicPlaneEdge->meshIntersection->P1);
//                    VerboseDislocationLoopNode(4,"ssd3.dMin= "<<ssd3.dMin<<std::endl;);

                    
                    SegmentSegmentDistance<dim-1> ssd(pPrevLocal,pNextLocal,*bndNode->periodicPlaneEdge->source,*bndNode->periodicPlaneEdge->sink);
                    if(ssd.dMin<FLT_EPSILON)
                    {
                        VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                        bndNode->set_P(VectorLowerDim(0.5*(ssd.x0+ssd.x1)));
                    }
                    else
                    {
                        bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                        this->network().danglingBoundaryLoopNodes.insert(bndNode);
                    }
                    
                }
                for(auto& bndNode : boundaryNext())
                {
                                        VerboseDislocationLoopNode(3,"Setting P for  NextBndNode "<<bndNode->tag()<<std::endl;);
                    const auto pPrev(bndNode->periodicPrev());
                    const auto pNext(bndNode->periodicNext());
                    assert(pPrev==this);
                    VerboseDislocationLoopNode(4,"pPrev= "<<pPrev->tag()<<" @ "<<pPrev->get_P().transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"pNext= "<<pNext->tag()<<" @ "<<pNext->get_P().transpose()<<std::endl;);

                    
                    const auto pPrevLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                    VerboseDislocationLoopNode(4,"pPrevLocal= "<<pPrevLocal.transpose()<<std::endl;);
                    const auto pNextLocal(this->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                    VerboseDislocationLoopNode(4,"pNextLocal= "<<pNextLocal.transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"periodicPlaneEdge->source= "<<bndNode->periodicPlaneEdge->source->transpose()<<std::endl;);
                    VerboseDislocationLoopNode(4,"bndNode->periodicPlaneEdge->sink= "<<bndNode->periodicPlaneEdge->sink->transpose()<<std::endl;);

//                    SegmentSegmentDistance<dim> ssd3(pPrev->get_P(),pNext->get_P(),bndNode->periodicPlaneEdge->meshIntersection->P0,bndNode->periodicPlaneEdge->meshIntersection->P1);
//                    VerboseDislocationLoopNode(4,"ssd3.dMin= "<<ssd3.dMin<<std::endl;);

                    
//                    std::cout<<"edges"<<std::endl;
//                    for(const auto& edge : bndNode->periodicPlaneEdge->patch->edges())
//                    {
//                        VerboseDislocationLoopNode(4,"periodicPlaneEdge->source= "<<edge->edgeID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<std::endl;);
//                    }
                    
                    SegmentSegmentDistance<dim-1> ssd(pPrevLocal,pNextLocal,*bndNode->periodicPlaneEdge->source,*bndNode->periodicPlaneEdge->sink);
                    if(ssd.dMin<FLT_EPSILON)
                    {
                        VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                        bndNode->set_P(VectorLowerDim(0.5*(ssd.x0+ssd.x1)));
                    }
                    else
                    {
                        bndNode->set_P(this->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                        this->network().danglingBoundaryLoopNodes.insert(bndNode);
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopNode<dim,corder,InterpolationType>::periodicPlanePatch() const
    {
        return _periodicPlanePatch;
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::updateGeometry()
    {

    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationLoopNode<dim,corder,InterpolationType>::isRemovable(const double& Lmin,const double& relAreaTh)
    {
        
        if(periodicPlaneEdge)
        {// a boundary node
            VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " NOT removable (onPlaneEdge)"<<std::endl;);
            return false;
        }
        else
        {// not a boundary node
            
            const auto pPrev(periodicPrev());
            const auto pNext(periodicNext());

            if((pPrev->get_P()-pNext->get_P()).norm()<this->network().networkRemesher.Lmax)
            {
                const double loopArea(this->loop()->slippedArea());
                if(loopArea>std::pow(Lmin,2))
                {
                    //                const auto& periodicPrevPos(periodicPrev()->get_P());
                    //                const auto& periodicNextPos(periodicNext()->get_P());
                    const double deltaArea(0.5*(this->get_P()-pPrev->get_P()).cross(pNext->get_P()-this->get_P()).norm());
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " loopArea= "<<loopArea<<std::endl;);
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " deltaArea= "<<deltaArea<<std::endl;);
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " relAreaTh= "<<relAreaTh<<std::endl;);
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " removable? "<<(deltaArea/loopArea<relAreaTh)<<std::endl;);
                    return deltaArea/loopArea<relAreaTh;
                }
                else if(loopArea>FLT_EPSILON)
                {
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " removable (small loop)"<<std::endl;);
                    return true;
                }
                else
                {
                    VerboseDislocationLoopNode(2,"  DislocationLoopNode "<<this->tag()<< " NOT removable (zero area loop)"<<std::endl;);
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationLoopNode<dim,corder,InterpolationType>::isMovableTo(const VectorDim& X) const
    {
        VerboseDislocationLoopNode(4,"checking if DislocationLoopNode "<<this->sID<< " isMovable:"<<std::endl;);

        bool isMovable=true;
        
        if(periodicPlaneEdge)
        {
            isMovable*=((X-this->get_P()).squaredNorm()<FLT_EPSILON);
        }
        
        if(isMovable && this->loop()->glidePlane)
        {
            isMovable*=this->loop()->glidePlane->contains(X);
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
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationLoopNode<dim,corder,InterpolationType>::verboseDislocationLoopNode=0;
}
#endif
