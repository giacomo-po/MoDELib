/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_cpp_
#define model_PeriodicGlidePlane_cpp_

#include <numbers>
#include <memory>
#include <math.h>

#include <PlanarMeshFace.h>
#include <MeshBoundarySegment.h>
#include <PeriodicGlidePlane.h>
//#include <PeriodicGlidePlane.hpp>
#include <PeriodicGlidePlaneFactory.h>

namespace model
{

    template<int dim>
    PeriodicPlaneEdge<dim>::PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                                              const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in,
                                              const short int& edgeID_in) :
    /* init */ patch(patch_in)
    /* init */,source(source_in)
    /* init */,sink(sink_in)
    /* init */,meshIntersection(meshIntersection_in)
    /* init */,edgeID(edgeID_in)
    //    /* init */,deltaShift(getDeltaShift(patch,meshIntersection))
    /* init */,deltaShift(meshIntersection->periodicShift())
    /* init */,next(nullptr)
    /* init */,prev(nullptr)
    /* init */,twin(nullptr)
    {
        //    std::cout<<"Creating PeriodicPlaneEdge "<<this->tag()<<std::endl;
        source->addLink(this);
        sink->addLink(this);
        if(twin)
        {
            patch->patchBoundary->removeUntwinnedEdge(twin);
        }
        else
        {
            patch->patchBoundary->addUntwinnedEdge(this);
        }
    }

    template<int dim>
    PeriodicPlaneEdge<dim>::~PeriodicPlaneEdge()
    {
        //        std::cout<<"Destroying PeriodicPlaneEdge "<<this->tag()<<std::endl;
        
        source->removeLink(this);
        sink->removeLink(this);
        if(next)
        {
            next->prev=nullptr;
        }
        if(prev)
        {
            prev->next=nullptr;
        }
        if(twin)
        {
            twin->twin=nullptr;
            patch->patchBoundary->addUntwinnedEdge(twin);
        }
        else
        {
            patch->patchBoundary->removeUntwinnedEdge(this);
        }
    }

    //    template<int dim>
    //    typename PeriodicPlaneEdge<dim>::VectorDim PeriodicPlaneEdge<dim>::getDeltaShift(const PeriodicPlanePatch<dim>* const patch,const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection)
    //    {
    //
    //        // VectorDim shift(VectorDim::Zero());
    //        // for(const PlanarMeshFace<dim>* const face : meshIntersection->faces)
    //        // {
    //        //     const auto parallelFaceID(patch->glidePlane->grain.region.parallelFaces().at(face->sID));
    //        //     const auto parallelFace(patch->glidePlane->grain.region.faces().at(parallelFaceID));
    //        //     shift+=(parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal();
    //
    //        // }
    //        // return shift;
    //        // std::cout<<"Starting to calculate delta shift "<<std::endl;
    //        VectorDim shift(VectorDim::Zero());
    //        for(const PlanarMeshFace<dim>* const face : meshIntersection->faces)
    //        {
    //            // std::cout<<"Normal corresponding to this face "<<face->sID<<face->outNormal().transpose()<<std::endl;
    //            const auto parallelFaceID(patch->glidePlane->grain.region.parallelFaces().at(face->sID));
    //            const auto parallelFace(patch->glidePlane->grain.region.faces().at(parallelFaceID));
    //            const double parallelFD(((parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal()).norm());
    //
    //            std::vector<VectorDim> vectorofNormals;
    //
    //
    //            for (const auto& pface : patch->glidePlane->grain.region.parallelFaces())
    //            {
    //                if (pface.first!= face->sID && pface.second!= face->sID)
    //                {
    //                    const auto otherFaceOutNormal(patch->glidePlane->grain.region.faces().at(pface.second)->outNormal());
    //                    vectorofNormals.push_back(otherFaceOutNormal);
    //                    // const auto shiftDirection((face->outNormal()-(face->outNormal().dot(otherFaceOutNormal)*otherFaceOutNormal)).normalized());
    //                    // std::cout<<"Shift direction "<<shiftDirection.transpose()<<std::endl;
    //                    // VectorDim faceShift(shiftDirection*parallelFD/face->outNormal().dot(shiftDirection));
    //                    // shift+=0.5*faceShift;
    //                }
    //            }
    //            GramSchmidt::orthoNormalize(vectorofNormals);
    //            VectorDim faceShift(face->outNormal());
    //            for (const auto& normal : vectorofNormals)
    //            {
    //                faceShift-=faceShift.dot(normal)*normal;
    //            }
    //            faceShift.normalize();
    //            shift-=faceShift*parallelFD/faceShift.dot(face->outNormal());
    //        }
    //        // std::cout<<"Finished to calculate delta shift "<<shift.transpose()<<std::endl;
    //
    //        return shift;
    //    }

    template<int dim>
    std::string PeriodicPlaneEdge<dim>::tag() const
    {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
      */
        return std::to_string(source->sID) + " " + std::to_string(sink->sID)+" "+std::to_string(patch->sID)+" ("+std::to_string(edgeID)+")";
    }


    /**********************************************************************/
    template<int dim>
    PeriodicPlaneNode<dim>::PeriodicPlaneNode(PeriodicPatchBoundary<dim>* const,const VectorLowerDim& pos) :
    /* init*/ VectorLowerDim(pos)
    {
        
    }

    template<int dim>
    void PeriodicPlaneNode<dim>::addLink(PeriodicPlaneEdgeType* const link)
    {
        
        if(link->sink->sID==this->sID)
        {// edge ends at this node, so link is an inLink
            
            // Update _patchConnectivities
            NodalConnectivityType& patchConnectivity(patchConnectivities()[link->patch->sID]);
            assert(patchConnectivity.inEdge==nullptr);
            patchConnectivity.inEdge=link;
            if(patchConnectivity.outEdge)
            {// and outEdge of the same patch is connected to this node
                assert(patchConnectivity.outEdge->prev==nullptr);
                assert(link->next==nullptr);
                patchConnectivity.outEdge->prev=link;
                link->next=patchConnectivity.outEdge;
            }
            
            // Update _neighborConnectivities
            NodalConnectivityType& neighborConnectivity(neighborConnectivities()[link->source->sID]);
            assert(neighborConnectivity.inEdge==nullptr || neighborConnectivity.inEdge==link);
            neighborConnectivity.inEdge=link;
            if(neighborConnectivity.outEdge)
            {// and outEdge of the same patch is connected to this node
                assert(neighborConnectivity.outEdge->twin==nullptr || neighborConnectivity.outEdge->twin==link);
                assert(link->twin==nullptr || link->twin==neighborConnectivity.outEdge);
                neighborConnectivity.outEdge->twin=link;
                link->twin=neighborConnectivity.outEdge;
            }
            
        }
        else if(link->source->sID==this->sID)
        {// edge starts at this node, so link is an outLink
            
            // Update _patchConnectivities
            NodalConnectivityType& patchConnectivity(patchConnectivities()[link->patch->sID]);
            assert(patchConnectivity.outEdge==nullptr);
            patchConnectivity.outEdge=link;
            if(patchConnectivity.inEdge)
            {
                assert(patchConnectivity.inEdge->next==nullptr);
                assert(link->prev==nullptr);
                patchConnectivity.inEdge->next=link;
                link->prev=patchConnectivity.inEdge;
            }
            
            // Update _neighborConnectivities
            NodalConnectivityType& neighborConnectivity(neighborConnectivities()[link->sink->sID]);
            assert(neighborConnectivity.outEdge==nullptr || neighborConnectivity.outEdge==link);
            neighborConnectivity.outEdge=link;
            if(neighborConnectivity.inEdge)
            {// and outEdge of the same patch is connected to this node
                assert(neighborConnectivity.inEdge->twin==nullptr || neighborConnectivity.inEdge->twin==link);
                assert(link->twin==nullptr || link->twin==neighborConnectivity.inEdge);
                neighborConnectivity.inEdge->twin=link;
                link->twin=neighborConnectivity.inEdge;
            }
        }
        else
        {
            assert(false && "CONNECTING LINK TO NON-INCIDENT NODE");
        }
        
    }

    template<int dim>
    void PeriodicPlaneNode<dim>::removeLink(PeriodicPlaneEdgeType* const link)
    {
        
        if(link->sink->sID==this->sID)
        {// edge ends at this node, so link is an inLink
            
            // Update _patchConnectivities
            auto patchIter(patchConnectivities().find(link->patch->sID));
            assert(patchIter!=patchConnectivities().end() && "PATCH NOT FOUND IN patchConnectivities");
            NodalConnectivityType& patchConnectivity(patchIter->second);
            assert(patchConnectivity.inEdge==link);
            patchConnectivity.inEdge=nullptr;
            if(patchConnectivity.inEdge==nullptr && patchConnectivity.outEdge==nullptr)
            {
                patchConnectivities().erase(patchIter);
            }
            
            // Update _neighborConnectivities
            auto neighborIter(neighborConnectivities().find(link->source->sID));
            assert(neighborIter!=neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
            NodalConnectivityType& neighborConnectivity(neighborIter->second);
            assert(neighborConnectivity.inEdge==link);
            neighborConnectivity.inEdge=nullptr;
            if(neighborConnectivity.inEdge==nullptr && neighborConnectivity.outEdge==nullptr)
            {
                neighborConnectivities().erase(neighborIter);
            }
            
        }
        else if(link->source->sID==this->sID)
        {// edge starts at this node, so link is an outLink
            
            // Update _patchConnectivities
            auto patchIter(patchConnectivities().find(link->patch->sID));
            assert(patchIter!=patchConnectivities().end() && "PATCH NOT FOUND IN patchConnectivities");
            NodalConnectivityType& patchConnectivity(patchIter->second);
            assert(patchConnectivity.outEdge==link);
            patchConnectivity.outEdge=nullptr;
            if(patchConnectivity.inEdge==nullptr && patchConnectivity.outEdge==nullptr)
            {
                patchConnectivities().erase(patchIter);
            }
            
            // Update _neighborConnectivities
            auto neighborIter(neighborConnectivities().find(link->sink->sID));
            assert(neighborIter!=neighborConnectivities().end() && "NEIGHBOR NOT FOUND IN neighborConnectivity");
            NodalConnectivityType& neighborConnectivity(neighborIter->second);
            assert(neighborConnectivity.outEdge==link);
            neighborConnectivity.outEdge=nullptr;
            if(neighborConnectivity.inEdge==nullptr && neighborConnectivity.outEdge==nullptr)
            {
                neighborConnectivities().erase(neighborIter);
            }
        }
        else
        {
            assert(false && "DISCONNECTING LINK FROM NON-INCIDENT NODE");
        }
    }

    template<int dim>
    const typename PeriodicPlaneNode<dim>::NodalConnectivityContainerType& PeriodicPlaneNode<dim>::patchConnectivities() const
    {
        return _patchConnectivities;
    }

    template<int dim>
    typename PeriodicPlaneNode<dim>::NodalConnectivityContainerType& PeriodicPlaneNode<dim>::patchConnectivities()
    {
        return _patchConnectivities;
    }

    template<int dim>
    const typename PeriodicPlaneNode<dim>::NodalConnectivityContainerType& PeriodicPlaneNode<dim>::neighborConnectivities() const
    {
        return _neighborConnectivities;
    }

    template<int dim>
    typename PeriodicPlaneNode<dim>::NodalConnectivityContainerType& PeriodicPlaneNode<dim>::neighborConnectivities()
    {
        return _neighborConnectivities;
    }

    template<int dim>
    typename PeriodicPlaneNode<dim>::InOutEdgeContainerType PeriodicPlaneNode<dim>::inEdges() const
    {
        InOutEdgeContainerType temp;
        for(const auto& connectivity : neighborConnectivities())
        {
            if(connectivity.second.inEdge && !connectivity.second.outEdge)
            {
                assert(connectivity.second.inEdge->twin==nullptr);
                temp.insert(connectivity.second.inEdge);
            }
        }
        return temp;
    }

    template<int dim>
    typename PeriodicPlaneNode<dim>::InOutEdgeContainerType PeriodicPlaneNode<dim>::outEdges() const
    {
        InOutEdgeContainerType temp;
        for(const auto& connectivity : neighborConnectivities())
        {
            if(connectivity.second.outEdge && !connectivity.second.inEdge)
            {
                assert(connectivity.second.outEdge->twin==nullptr);
                temp.insert(connectivity.second.outEdge);
            }
        }
        return temp;
    }

    /**********************************************************************/
    template<int dim>
    PeriodicPlanePatch<dim>::PeriodicPlanePatch(PeriodicPatchBoundary<dim>* const periodicPatchBoundary_in,
                                                const VectorDim& shift_in) :
    /* init */ patchBoundary(periodicPatchBoundary_in)
    /* init */,shift(shift_in)
    /* init */,glidePlane(patchBoundary->getGlidePlane(shift))
    {
        //    std::cout<<cyanColor<<"Creating patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
        
        if(isRightHandedBoundary(glidePlane->meshIntersections,*patchBoundary->referencePlane))
        {
            addMeshIntersections(glidePlane->meshIntersections);
        }
        else
        {
            BoundingMeshSegments<dim> flippedTemp;
            for (auto it = glidePlane->meshIntersections.rbegin(); it != glidePlane->meshIntersections.rend(); ++it)
            {
                flippedTemp.emplace_back(new MeshBoundarySegment<dim>((*it)->P1,(*it)->P0,(*it)->faces));
            }
            addMeshIntersections(flippedTemp);
        }
        
        patchBoundary->updateBoundaries();
        
    }

    template<int dim>
    PeriodicPlanePatch<dim>::~PeriodicPlanePatch()
    {
        //                    std::cout<<cyanColor<<"Destroying patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
        this->clear(); // call before updateBoundaries
        patchBoundary->updateBoundaries();
    }

    template<int dim>
    bool PeriodicPlanePatch<dim>::isRightHandedBoundary(const BoundingMeshSegments<dim>& bnd,const Plane<dim>& plane)
    {
        if(bnd.size())
        {
            VectorDim nA(VectorDim::Zero());
            const VectorDim P0(bnd.front()->P0);
            for(const auto& seg : bnd)
            {
                nA+= 0.5*(seg->P0-P0).cross(seg->P1-seg->P0);
            }
            return nA.dot(plane.unitNormal)>0.0;
        }
        else
        {
            return false;
        }
    }

    template<int dim>
    void PeriodicPlanePatch<dim>::addMeshIntersections(const BoundingMeshSegments<dim>& bms)
    {
        //
        //        std::vector<std::tuple<std::shared_ptr<PeriodicPlaneNode<dim>>,std::shared_ptr<PeriodicPlaneNode<dim>>,std::shared_ptr<MeshBoundarySegment<dim>>,double>> tempEdges;
        //        for(size_t k=0;k<bms.size();++k)
        //        {// compute 2d points on PeriodicPlane, and coneect them
        //            std::shared_ptr<PeriodicPlaneNode<dim>> source(patchBoundary->getSharedNode(bms[k]->P0,shift));
        //            std::shared_ptr<PeriodicPlaneNode<dim>>   sink(patchBoundary->getSharedNode(bms[k]->P1,shift));
        //            tempEdges.emplace_back(source,sink,bms[k],0.0);
        //        }
        //
        //        VectorLowerDim C(VectorLowerDim::Zero());
        //        for(const auto& tup : tempEdges)
        //        {
        //            C+=*std::get<0>(tup);
        //        }
        //        C/=tempEdges.size(); // position of center of patch
        //
        //        std::pair<double,size_t> minAngle(std::numbers::pi,tempEdges.size());
        //        for(size_t k=0;k<tempEdges.size();++k)
        //        {
        //            const auto x(*std::get<0>(tempEdges[k])-C);
        //            const double angle(std::fabs(atan2(x(1),x(0))));
        //            if(angle<minAngle.first)
        //            {
        //                minAngle.first=angle;
        //                minAngle.second=k;
        //            }
        //            std::get<3>(tempEdges[k])=angle;
        //        }
        //
        //        // store edges starting from source of minimum angle
        //        for(size_t s=0;s<tempEdges.size();++s)
        //        {
        //            const short int k((minAngle.second+s)%tempEdges.size());
        //            this->emplace_back(new PeriodicPlaneEdge<dim>(this,std::get<0>(tempEdges[k]),std::get<1>(tempEdges[k]),std::get<2>(tempEdges[k]),s));
        //        }
        
        //
        for(size_t k=0;k<bms.size();++k)
        {// compute 2d points on PeriodicPlane, and coneect them
            std::shared_ptr<PeriodicPlaneNode<dim>> source(patchBoundary->getSharedNode(bms[k]->P0,shift));
            std::shared_ptr<PeriodicPlaneNode<dim>>   sink(patchBoundary->getSharedNode(bms[k]->P1,shift));
            //this->emplace_back(new PeriodicPlaneEdge<dim>(this,source,sink,bms[k]));
            this->emplace_back(new PeriodicPlaneEdge<dim>(this,source,sink,bms[k],k));
        }
    }

    template<int dim>
    const typename PeriodicPlanePatch<dim>::PeriodicPlaneEdgeContainerType& PeriodicPlanePatch<dim>::edges() const
    {
        return *this;
    }

    template<int dim>
    std::shared_ptr<PeriodicPlaneEdge<dim>> PeriodicPlanePatch<dim>::edgeOnFaces(const std::set<const PlanarMeshFace<dim>*>& faces) const
    {
        for(const auto& edge : edges())
        {
            if(edge->meshIntersection->faces==faces)
            {
                return edge;
            }
        }
        return std::shared_ptr<PeriodicPlaneEdge<dim>>(nullptr);
    }

    template<int dim>
    int PeriodicPlanePatch<dim>::contains(const VectorLowerDim& test) const
    {
        //                const auto& outerBoundary(outerBoundaries()[0]);
        size_t i, j, c = 0;
        for (i = 0, j = edges().size()-1; i < edges().size(); j = i++)
        {
            const VectorLowerDim& Pi(*edges()[i]->source);
            const VectorLowerDim& Pj(*edges()[j]->source);
            
            if ( ((Pi(1)>test(1)) != (Pj(1)>test(1))) &&
                (test(0) < (Pj(0)-Pi(0)) * (test(1)-Pi(1)) / (Pj(1)-Pi(1)) + Pi(0)) )
            {
                c = !c;
            }
        }
        return c;
    }

    template<int dim>
    PeriodicPatchBoundary<dim>::PeriodicPatchBoundary(GlidePlaneFactory<dim>& glidePlaneFactory_in, const std::shared_ptr<GlidePlane<dim>>& referencePlane_in) :
    /* init */ glidePlaneFactory(glidePlaneFactory_in)
    /* init */,referencePlane(referencePlane_in)
    //    /* init */,L2G(getL2G(referencePlane->unitNormal))
    {
        
        assert(referencePlane->meshIntersections.size()>=3);
    }


    template<int dim>
    typename PeriodicPatchBoundary<dim>::BoundariesContainerType& PeriodicPatchBoundary<dim>::outerBoundaries()
    {
        return _outerBoundaries;
    }

    template<int dim>
    const typename PeriodicPatchBoundary<dim>::BoundariesContainerType& PeriodicPatchBoundary<dim>::outerBoundaries() const
    {
        return _outerBoundaries;
    }

    template<int dim>
    typename PeriodicPatchBoundary<dim>::BoundariesContainerType& PeriodicPatchBoundary<dim>::innerBoundaries()
    {
        return _innerBoundaries;
    }

    template<int dim>
    const typename PeriodicPatchBoundary<dim>::BoundariesContainerType& PeriodicPatchBoundary<dim>::innerBoundaries() const
    {
        return _innerBoundaries;
    }

    template<int dim>
    const typename PeriodicPatchBoundary<dim>::UntwinnedEdgeContainerType& PeriodicPatchBoundary<dim>::untwinnedEdges() const
    {
        return *this;
    }

    template<int dim>
    typename PeriodicPatchBoundary<dim>::UntwinnedEdgeContainerType& PeriodicPatchBoundary<dim>::untwinnedEdges()
    {
        return *this;
    }

    template<int dim>
    void PeriodicPatchBoundary<dim>::addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
    {
        untwinnedEdges().insert(link);
    }

    template<int dim>
    void PeriodicPatchBoundary<dim>::removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
    {
        const size_t erased(untwinnedEdges().erase(link));
        if(erased!=1)
        {
            throw std::runtime_error("COULD NOT ERASE LINK FROM BOUNDARYLINKS");
        }
//        assert(erased==1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
    }

    template<int dim>
    typename PeriodicPatchBoundary<dim>::NodeCointainerType& PeriodicPatchBoundary<dim>::nodes()
    {
        return *this;
    }

    template<int dim>
    const typename PeriodicPatchBoundary<dim>::NodeCointainerType& PeriodicPatchBoundary<dim>::nodes() const
    {
        return *this;
    }

    template<int dim>
    int PeriodicPatchBoundary<dim>::isInsideOuterBoundary(const VectorLowerDim& test)
    {
        if(isCompact())
        {
            const auto& outerBoundary(outerBoundaries()[0]);
            size_t i, j, c = 0;
            for (i = 0, j = outerBoundary.size()-1; i < outerBoundary.size(); j = i++)
            {
                const VectorLowerDim& Pi(*outerBoundary[i]->source);
                const VectorLowerDim& Pj(*outerBoundary[j]->source);
                
                if ( ((Pi(1)>test(1)) != (Pj(1)>test(1))) &&
                    (test(0) < (Pj(0)-Pi(0)) * (test(1)-Pi(1)) / (Pj(1)-Pi(1)) + Pi(0)) )
                {
                    c = !c;
                }
            }
            return c;
        }
        else
        {
            return 0;
        }
    }

    template<int dim>
    bool PeriodicPatchBoundary<dim>::isCompact() const
    {
        return outerBoundaries().size()==1 && innerBoundaries().size()==0;
    }

    template<int dim>
    typename PeriodicPatchBoundary<dim>::VectorDim PeriodicPatchBoundary<dim>::rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
    {// TO DO: THIS SHOULD USE THE 2D POINTS
        if(bnd.size()>=3)
        {
            VectorDim nA(VectorDim::Zero());
            const VectorDim P0(bnd.front()->meshIntersection->P0-(bnd.front())->patch->shift);
            for(const auto& seg : bnd)
            {
                VectorDim temp_shift=seg->patch->shift;
                VectorDim P0_temp=seg->meshIntersection->P0-temp_shift;
                VectorDim P1_temp=seg->meshIntersection->P1-temp_shift;
                nA+= 0.5*(P0_temp-P0).cross(P1_temp-P0_temp);
            }
            return nA;
        }
        else
        {
            return VectorDim::Zero();
        }
    }

    template<int dim>
    double PeriodicPatchBoundary<dim>::rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
    {
        return rightHandedNormal(bnd).dot(referencePlane->unitNormal);
    }

    template<int dim>
    bool PeriodicPatchBoundary<dim>::isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
    {
        return rightHandedArea(bnd)>=0.0;
    }

    template<int dim>
    std::shared_ptr<PeriodicPlaneNode<dim> > PeriodicPatchBoundary<dim>::getSharedNode(const VectorDim& pointDim,const VectorDim& shift)
    {
        //        const VectorLowerDim point(referencePlane->localPosition(pointDim-shift));
        return nodes().getFromKey(referencePlane->localPosition(pointDim-shift));
    }

    template<int dim>
    void PeriodicPatchBoundary<dim>::createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy)
    {
        BoundaryContainerType temp;
        temp.reserve(untwinnedCopy.size());
        while(true)
        {
            temp.push_back(currentEdge);
            const size_t erased(untwinnedCopy.erase(currentEdge));
            if(erased!=1)
            {
                
                std::cout<<"Trying to erase "<<currentEdge->tag()<<std::endl;
                std::cout<<"untwinnedCopy is"<<std::endl;
                for(const auto& edgePtr : untwinnedCopy)
                {
                    std::cout<<"    "<<edgePtr->tag()<<std::endl;
                }
                
                std::ofstream pointsFile("points.txt");
                std::cout<<"points"<<std::endl;
                for(const auto& node : this->nodes())
                {
                    if(!node.second.expired())
                    {
                        pointsFile<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                        std::cout<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                    }
                }
                
                
                //                std::ofstream edgesFile("edges.txt");
                //                for(const auto& patch : patches())
                //                {
                //                    for(const auto& edge : patch.second->edges())
                //                    {
                //                        edgesFile<<edge->tag()<<std::endl;
                //                    }
                //                }
                
                throw std::runtime_error("could not find link in untwinnedEdges 2");
                
//                assert(erased==1 && "could not find link in untwinnedEdges 2");
            }
            if(temp.back()->sink->sID==temp.front()->source->sID)
            {
                break;
            }
            currentEdge=currentEdge->next;
            while(currentEdge->twin)
            {
                currentEdge=currentEdge->twin->next;
            }
            
        }
        
        if(temp.size())
        {// if temp is not empty a loop was closed.
            if(this->isRightHandedBoundary(temp))
            {// assign temp to _outerBoundaries
                this->_outerBoundaries.push_back(temp);
            }
            else
            {// assign temp to _innerBoundaries
                this->_innerBoundaries.push_back(temp);
                
            }
        }
    }

    //Other version (for the case where "must start from node with only one edge in and one edge out" is not possible)
    template<int dim>
    void PeriodicPatchBoundary<dim>::updateBoundaries()
    {
        //    std::cout<<"Updating OuterBoundary"<<std::endl;
        this->outerBoundaries().clear();
        this->innerBoundaries().clear();
        UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
        while(untwinnedCopy.size())
        {
            size_t untwinnedSize(untwinnedCopy.size());
            for(const auto& edge : untwinnedCopy)
            {
                if(edge->source->outEdges().size()==1 && edge->source->inEdges().size()==1) //Discuss this with Dr. Po(Yash)
                {// must start from node with only one edge in and one edge out
                    createNewBoundary(edge,untwinnedCopy);
                    break;
                }
            }
            if (untwinnedCopy.size() && untwinnedCopy.size()==untwinnedSize) //No node exists with only one edge in and one edge out
            {
                createNewBoundary(*untwinnedCopy.begin(),untwinnedCopy);
                
            }
        }
    }

    template<int dim>
    std::shared_ptr<PeriodicPlanePatch<dim>> PeriodicPatchBoundary<dim>::getPatch(const VectorDim& shift)
    {
        return patches().getFromKey(shift);
    }

    template<int dim>
    const typename PeriodicPatchBoundary<dim>::PatchContainerType& PeriodicPatchBoundary<dim>::patches() const
    {
        return *this;
    }

    template<int dim>
    typename PeriodicPatchBoundary<dim>::PatchContainerType& PeriodicPatchBoundary<dim>::patches()
    {
        return *this;
    }

    template<int dim>
    GlidePlaneKey<dim> PeriodicPatchBoundary<dim>::getGlidePlaneKey(const VectorDim& shift)
    {
        return GlidePlaneKey<dim>(referencePlane->P+shift,referencePlane->n);
    }

    template<int dim>
    std::shared_ptr<GlidePlane<dim>> PeriodicPatchBoundary<dim>::getGlidePlane(const VectorDim& shift)
    {
        return glidePlaneFactory.getFromKey(getGlidePlaneKey(shift));
    }

    template<int dim>
    PeriodicGlidePlane<dim>::PeriodicGlidePlane(PeriodicGlidePlaneFactoryType* const pgpf,
                                                const KeyType& key_in) :
    /* init */ PeriodicPatchBoundary<dim>(pgpf->glidePlaneFactory,pgpf->glidePlaneFactory.getFromKey(key_in))
    //    /* init */ PeriodicPatchBoundary<dim>(pgpf->glidePlaneFactory.getFromKey(key_in))
    /* init */,glidePlaneFactory(pgpf->glidePlaneFactory)
    /* init */,periodicGlidePlaneFactory(*pgpf)
    /* init */,referencePlane(pgpf->glidePlaneFactory.getFromKey(key_in))
    {
        
    }

    template<int dim>
    typename PeriodicGlidePlane<dim>::VectorDim PeriodicGlidePlane<dim>::findPatch(const typename PeriodicGlidePlane<dim>::VectorLowerDim& P0,const VectorDim& guessShift)
    {
        
        PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory,referencePlane);
        std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> currentPatches; // keep patches alive during search
        
        
        // Set up a starting point
        std::shared_ptr<PeriodicPlanePatch<dim>> lastPatch(patchBoundary.getPatch(guessShift));
        currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
        const std::shared_ptr<GlidePlane<dim>>& centralPlane(lastPatch->glidePlane);
        VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
        for(const auto& seg : centralPlane->meshIntersections) // problem is here for referencePlane not cutting mesh
        {
            insideReferencePoint+=referencePlane->localPosition(seg->P0-lastPatch->shift);
        }
        insideReferencePoint/=centralPlane->meshIntersections.size();
        assert(lastPatch->contains(insideReferencePoint));
        
        // Find initial point
        //        const VectorLowerDim P0(polyPoints[0]);
        while(!lastPatch->contains(P0))
        {
            std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
            for(const auto& bndEdge : patchBoundary.untwinnedEdges())
            {// loop over outer boundaries and holes
                SegmentSegmentDistance<dim-1> ssd(P0,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
                if(ssd.dMin<FLT_EPSILON)
                {// intersection with current boundary found
                    //                    std::cout<<ssd.x0.transpose()<<std::endl;
                    crossdEdges.insert(bndEdge);
                }
            }
            
            switch (crossdEdges.size())
            {
                case 1:
                {
                    assert(lastPatch.get()==(*crossdEdges.begin())->patch);
                    const VectorDim newShift((*crossdEdges.begin())->deltaShift+(*crossdEdges.begin())->patch->shift);
                    //                    const auto shiftIter(usedPatches.find(newShift));
                    lastPatch=patchBoundary.getPatch(newShift);
                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                    //                    if(shiftIter==usedPatches.end())
                    //                    {
                    //                        lastPatch=usedPatches.emplace(newShift,new PeriodicPlanePatch<dim>(*this,newShift)).first->second;
                    //                    }
                    //                    else
                    //                    {
                    //                        lastPatch=shiftIter->second;
                    //                    }
                    //                    usedPatches.emplace(new PeriodicPlanePatch<dim>(*this,));
                    //getPatch(getShift(**crossdEdges.begin())+(*crossdEdges.begin())->patch->shift);
                    break;
                }
                    
                case 2:
                {
                    assert(lastPatch.get()==(*crossdEdges.begin())->patch);
                    assert((*crossdEdges.begin())->patch==(*crossdEdges.rbegin())->patch);
                    const VectorDim newShift((*crossdEdges.begin())->deltaShift+(*crossdEdges.rbegin())->deltaShift+(*crossdEdges.begin())->patch->shift);
                    lastPatch=patchBoundary.getPatch(newShift);
                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                    //                    const auto shiftIter(usedPatches.find(newShift));
                    //                    if(shiftIter==usedPatches.end())
                    //                    {
                    //                        lastPatch=usedPatches.emplace(newShift,new PeriodicPlanePatch<dim>(*this,newShift)).first->second;
                    //                    }
                    //                    else
                    //                    {
                    //                        lastPatch=shiftIter->second;
                    //                    }
                    //
                    //
                    //                    usedPatches.emplace(new PeriodicPlanePatch<dim>(*this,getShift(**crossdEdges.begin())+getShift(**crossdEdges.rbegin())+(*crossdEdges.begin())->patch->shift));
                    //                    getPatch(getShift(**crossdEdges.begin())+getShift(**crossdEdges.rbegin())+(*crossdEdges.begin())->patch->shift);
                    break;
                }
                    
                default:
                {
                    std::cout<<"untwinnedEdges().size()="<<patchBoundary.untwinnedEdges().size()<<std::endl;
                    std::cout<<"crossdEdges.size()="<<crossdEdges.size()<<std::endl;
                    assert(false && "1 or 2 edges must be crossed");
                    break;
                }
            }
        }
        assert(lastPatch->contains(P0));
        return lastPatch->shift;
    }

    template<int dim>
    void PeriodicGlidePlane<dim>::print()
    {
        std::cout<<"PeriodiPlane nodes:"<<std::endl;
        for(const auto& node : this->nodes())
        {
            if(!node.second.expired())
            {
                std::cout<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
            }
        }
        
    }

template<int dim>
std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> PeriodicGlidePlane<dim>::filledPatches(const std::vector<VectorDim>& patchShifts)
{/*\returns a vector of patch-shifts (including startShifts) which closes the "holes" on pgp
  */
    // Collect patches of the loop. We create a temporary PeriodicPatchBoundary for this, having only the patches of this loop
    PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory,referencePlane);
    std::set<std::shared_ptr<PeriodicPlanePatch<dim>>> tempPatches;
    for(const auto& ps : patchShifts)
    {
        tempPatches.insert(patchBoundary.getPatch(ps));
    }
    // There may be "holes" in the patches. Close them
    while(patchBoundary.innerBoundaries().size())
    {
        const PeriodicPlaneEdge<dim>* const holeEdge(*patchBoundary.innerBoundaries().front().begin());
        tempPatches.insert(patchBoundary.getPatch(holeEdge->deltaShift+holeEdge->patch->shift));
    }

    // Now grab corresponding patches from this->patchBoundary
    std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> temp;
    for(const auto& patch : tempPatches)
    {
        temp.push_back(this->getPatch(patch->shift));
    }
    return temp;
}

    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO() :
    /* init */ glidePlaneID(0)
    /* init */,patchBoundaryID(0)
    /* init */,shift(VectorDim::Zero())
    {
        
    }

    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch) :
    /* init */ glidePlaneID(patch.glidePlane->sID)
    /* init */,patchBoundaryID(patch.patchBoundary->sID)
    /* init */,shift(patch.shift)
    {
        
    }

    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO(std::stringstream& ss) :
    /* init */ glidePlaneID(0)
    /* init */,patchBoundaryID(0)
    /* init */,shift(VectorDim::Zero())
    {
        ss>>glidePlaneID;
        ss>>patchBoundaryID;
        for(int d=0;d<dim;++d)
        {
            ss>>shift(d);
        }
    }

    template struct PeriodicPlaneEdge<3>;
    template struct PeriodicPatchBoundary<3>;
    template class PeriodicPlaneNode<3>;
    template struct PeriodicPlanePatch<3>;
    template class PeriodicGlidePlane<3>;
    template struct PeriodicPlanePatchIO<3>;


}
#endif
