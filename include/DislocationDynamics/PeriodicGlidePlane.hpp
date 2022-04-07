/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_cpp_
#define model_PeriodicGlidePlane_cpp_


#include <PeriodicGlidePlane.h>

#include <math.h>


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
    /* init */,deltaShift(getDeltaShift(patch,meshIntersection))
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
    
    template<int dim>
    typename PeriodicPlaneEdge<dim>::VectorDim PeriodicPlaneEdge<dim>::getDeltaShift(const PeriodicPlanePatch<dim>* const patch,const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection)
    {
        
        // VectorDim shift(VectorDim::Zero());
        // for(const PlanarMeshFace<dim>* const face : meshIntersection->faces)
        // {
        //     const auto parallelFaceID(patch->glidePlane->grain.region.parallelFaces().at(face->sID));
        //     const auto parallelFace(patch->glidePlane->grain.region.faces().at(parallelFaceID));
        //     shift+=(parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal();
            
        // }
        // return shift;
        // std::cout<<"Starting to calculate delta shift "<<std::endl;
        VectorDim shift(VectorDim::Zero());
        for(const PlanarMeshFace<dim>* const face : meshIntersection->faces)
        {
            // std::cout<<"Normal corresponding to this face "<<face->sID<<face->outNormal().transpose()<<std::endl;
            const auto parallelFaceID(patch->glidePlane->grain.region.parallelFaces().at(face->sID));
            const auto parallelFace(patch->glidePlane->grain.region.faces().at(parallelFaceID));
            const double parallelFD(((parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal()).norm());

            std::vector<VectorDim> vectorofNormals;
            

            for (const auto& pface : patch->glidePlane->grain.region.parallelFaces())
            {
                if (pface.first!= face->sID && pface.second!= face->sID)
                {
                    const auto otherFaceOutNormal(patch->glidePlane->grain.region.faces().at(pface.second)->outNormal());
                    vectorofNormals.push_back(otherFaceOutNormal);
                    // const auto shiftDirection((face->outNormal()-(face->outNormal().dot(otherFaceOutNormal)*otherFaceOutNormal)).normalized());
                    // std::cout<<"Shift direction "<<shiftDirection.transpose()<<std::endl;
                    // VectorDim faceShift(shiftDirection*parallelFD/face->outNormal().dot(shiftDirection));
                    // shift+=0.5*faceShift;
                }
            }
            GramSchmidt::orthoNormalize(vectorofNormals);
            VectorDim faceShift(face->outNormal());
            for (const auto& normal : vectorofNormals)
            {
                faceShift-=faceShift.dot(normal)*normal;
            }
            faceShift.normalize();
            shift-=faceShift*parallelFD/faceShift.dot(face->outNormal());
        }
        // std::cout<<"Finished to calculate delta shift "<<shift.transpose()<<std::endl;

        return shift;
    }
    
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
//        std::pair<double,size_t> minAngle(M_PI,tempEdges.size());
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
    int PeriodicPlanePatch<dim>::contains(const VectorLowerDim& test)
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
    
//    /**********************************************************************/
//    template<int dim>
//    typename PeriodicPatchBoundary<dim>::MatrixDim PeriodicPatchBoundary<dim>::getL2G(VectorDim z)
//    {
//        //            const double xNorm(x.norm());
//        const double zNorm(z.norm());
//        assert(zNorm>FLT_EPSILON);
//        z/=zNorm;
//
//        VectorDim x(VectorDim::UnitX().cross(z));
//        double xNorm(x.norm());
//        if(xNorm>FLT_EPSILON)
//        {
//            x=x/xNorm;
//        }
//        else
//        {
//            x=VectorDim::UnitY().cross(z);
//            xNorm=x.norm();
//            if(xNorm>FLT_EPSILON)
//            {
//                x=x/xNorm;
//            }
//            else
//            {
//                x=VectorDim::UnitZ().cross(z);
//                xNorm=x.norm();
//                if(xNorm>FLT_EPSILON)
//                {
//                    x=x/xNorm;
//                }
//                else
//                {
//                    assert(false && "CANNOT FIND VECTOR ORTHOGONAL TO z");
//                }
//            }
//        }
//
//        assert(std::fabs(x.norm()-1.0)<FLT_EPSILON);
//        assert(std::fabs(z.norm()-1.0)<FLT_EPSILON);
//        assert(fabs(x.dot(z)<FLT_EPSILON));
//        MatrixDim temp(Eigen::Matrix3d::Identity());
//        temp.col(2)=z;
//        temp.col(0)=x;
//        temp.col(1)=temp.col(2).cross(temp.col(0));
//        return temp;
//    }
    
    template<int dim>
    PeriodicPatchBoundary<dim>::PeriodicPatchBoundary(GlidePlaneFactory<dim>& glidePlaneFactory_in, const std::shared_ptr<GlidePlane<dim>>& referencePlane_in) :
    /* init */ glidePlaneFactory(glidePlaneFactory_in)
    /* init */,referencePlane(referencePlane_in)
//    /* init */,L2G(getL2G(referencePlane->unitNormal))
    {
        
        assert(referencePlane->meshIntersections.size()>=3);
    }
    
    //        ~PeriodicPatchBoundary()
    //        {
    //            std::cout<<"DESTROYING PeriodicPatchBoundary"<<std::endl;
    //        }
    
    
    
    
    
    

    
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
        assert(erased==1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
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
    
//    template<int dim>
//    typename PeriodicPatchBoundary<dim>::VectorLowerDim PeriodicPatchBoundary<dim>::getLocalPosition(const VectorDim& point,const VectorDim& shift) const
//    {
//        const VectorDim pointLocal(L2G.transpose()*(point-shift-referencePlane->P));
//        if(fabs(pointLocal(2))>FLT_EPSILON)
//        {
//            std::cout<<"point"<<point.transpose()<<std::endl;
//            std::cout<<"shift"<<shift.transpose()<<std::endl;
//            std::cout<<"referencePlane->P"<<referencePlane->P.transpose()<<std::endl;
//            std::cout<<"L2G=\n"<<L2G<<std::endl;
//            std::cout<<"pointLocal"<<pointLocal.transpose()<<std::endl;
//            assert(false && "local point has non-zero z-coordinate");
//        }
//        return pointLocal.template segment<2>(0);
//    }
//
//    template<int dim>
//    typename PeriodicPatchBoundary<dim>::VectorDim PeriodicPatchBoundary<dim>::getGlobalPosition(const VectorLowerDim& point) const
//    {// terurns the position on the plane in global goordinates
//        return L2G.template block<dim,dim-1>(0,0)*point+referencePlane->P;
//    }
    
//    template<int dim>
//    std::shared_ptr<PeriodicPlaneNode<dim> > PeriodicPatchBoundary<dim>::getSharedNode(const VectorDim& pointDim,const VectorDim& shift)
//    {
//        const VectorLowerDim point(referencePlane->localPosition(pointDim-shift));
//        const auto iter(nodes().find(point));
//        if(iter==nodes().end())
//        {// point does not exist
//            //                std::shared_ptr<PeriodicPlaneNode<dim>> newNode(new PeriodicPlaneNode<dim>(point));
//            //                nodes().emplace(point,newNode.get());
//            return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
//            //                return newNode;
//        }
//        else
//        {// point exists
//            if(iter->second.expired())
//            {// node deleted elsewhere
//                nodes().erase(iter);
//                return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
//            }
//            else
//            {
//                return iter->second.lock();
//            }
//        }
//    }
    
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
                
                
                assert(erased==1 && "could not find link in untwinnedEdges 2");
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
        //            }
    }

//Original Version    
    // template<int dim>
    // void PeriodicPatchBoundary<dim>::updateBoundaries()
    // {
    //             //    std::cout<<"Updating OuterBoundary"<<std::endl;
    //     this->outerBoundaries().clear();
    //     this->innerBoundaries().clear();
    //     UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
    //     while(untwinnedCopy.size())
    //     {
    //         for(const auto& edge : untwinnedCopy)
    //         {
    //             if(edge->source->outEdges().size()==1 && edge->source->inEdges().size()==1) //Discuss this with Dr. Po(Yash)
    //             {// must start from node with only one edge in and one edge out
    //                 createNewBoundary(edge,untwinnedCopy);
    //                 break;
    //             }
    //         }
    //     }
    //             //    std::cout<<"Updated OuterBoundary"<<std::endl;

    // }

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
                //    std::cout<<"Updated OuterBoundary"<<std::endl;

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

    
   
    
    /**********************************************************************/
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
    template<typename T>
     std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints)
    {
        std::vector<std::pair<VectorLowerDim,const T* const>> lowerPolyPoints;
        for(const auto& point : polyPoints)
        {
            assert(this->referencePlane->contains(point.first) && "reference plane does not cointain point");
            lowerPolyPoints.push_back(std::make_pair(referencePlane->localPosition(point.first-VectorDim::Zero()),point.second));
        }
        return polygonPatchIntersection(lowerPolyPoints);
    }


//     //Function updated to account for the insertion of the nodes on the periodically opposite faces
// template<int dim>
//     template<typename T>
//     std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints)
//     {
//         std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short int>,const T* const>> ppi; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
 
//         if (polyPoints.size())
//         {
//             PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory, referencePlane);
//             auto lastPatch(patchBoundary.getPatch(findPatch(polyPoints[0].first, VectorDim::Zero())));
//             // const size_t loopID((*polyPoints.begin()).second->loop()->sID);
//             // const size_t runID((*polyPoints.begin()).second->network().simulationParameters.runID);

//             //                std::ofstream polyFile("poly.txt");
//             //                for(const auto& node : polyPoints)
//             //                {
//             //                    polyFile<<"    "<<node.transpose()<<std::endl;
//             //                }

//             //        int removeMe=0;

//             //        ppi.emplace_back(polyPoints[0],lastPatch,nullptr);

//             // Find other points
//             for (size_t k = 0; k < polyPoints.size(); ++k)
//             {
//                 std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> currentPatches;
//                 //            std::cout<<"k point ="<<k<<std::endl;
//                 currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
//                 const auto startPoint(polyPoints[k]);
//                 const auto endPoint(k == polyPoints.size() - 1 ? polyPoints[0] : polyPoints[k + 1]);
//                 while (true)
//                 {
//                     std::map<const PeriodicPlanePatch<dim> *, std::vector<std::pair<const VectorLowerDim, const PeriodicPlaneEdge<dim> *const>>> crossdEdges;

//                     for (const auto &bndEdge : patchBoundary.untwinnedEdges())
//                     { // loop over outer boundaries and holes
//                         SegmentSegmentDistance<dim - 1> ssd(startPoint.first, endPoint.first, *bndEdge->source, *bndEdge->sink);
//                         // std::cout << " Determining intersection between " << startPoint.second->tag() << " and " << endPoint.second->tag()<<" dMin "<<ssd.dMin <<" den is "<< std::flush;
//                         // const double alignemntAngle(((endPoint.first - startPoint.first).normalized()).dot((*bndEdge->sink - *bndEdge->source).normalized()));
//                         // const double maxAngleAlignment(cos(1.0e-4 * M_PI / 180.0)); // Alignment angle is fixed to 10 degrees

//                         if (ssd.isIntersecting(FLT_EPSILON)) 
//                         { // intersection with current boundary found
//                             // std::cout << std::setprecision(15) << " Intersection found " << ssd.t << " " << ssd.u << std::endl;
//                             crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
//                         }
//                         //Changed to allow for the intersection of the near parallel segment with the boundary
//                         // else
//                         // {
//                         //     if (startPoint.second && endPoint.second)
//                         //     {
//                         //         // Check the alignment of the links with the bndEdge Source
//                         //         //                    	std::cout<<"Coming here to check for 3D intersection"<<std::endl;
//                         //         const double maxAngleAlignment(cos(5.0 * M_PI / 180.0)); //Alignment angle is fixed to 10 degrees
//                         //         // std::cout<<" Alignemnt Angle =>maxAngleAlignment"<<alignemntAngle<<" "<<maxAngleAlignment<<" for nodes "<<startPoint.second->tag()
//                         //         // <<"=>"<<endPoint.second->tag()<<" ssdDmin "<<ssd.dMin<<std::endl;

//                         //         if (fabs(alignemntAngle) > maxAngleAlignment && ssd.dMin < 1000*FLT_EPSILON) //If very close to the boundary and aligned only then check in 3D
//                         //         {
//                         //             //The value of 1000 is determined by trial and error
//                         //             //Determine the patches of startPoint and endPoint
//                         //             if (startPoint.second->periodicPlanePatch() != endPoint.second->periodicPlanePatch())
//                         //             {
//                         //                 //The two loop nodes belong to different patches, a node is needed to be inserted on the edge
//                         //                 //This is only valid with the if for the alignment angle and the ssd.dMin
//                         //                 crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
//                         //             }
//                         //         }
//                         //     }
//                         // }
                        
//                     }

//                     if (crossdEdges.size() == 0)
//                     { // No untwinned edges have been crossed, break and move to next pair of polyPoints
//                         ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
//                         //                    printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
//                         break;
//                     }
//                     else
//                     {
//                         for (const auto &pair : crossdEdges)
//                         { // loop over crossed patches
//                             //VectorDim shift(pair.first->shift);
//                             //                        std::cout<<"crossed patches="<<crossdEdges.size()<<std::endl;
//                             assert((pair.second.size()>0 && pair.second.size()<=2) && "At max two periodic plane edges can be intersected (2D geometry constraint)");

//                             // if (pair.second.size()==2)
//                             // {
//                             //     std::cout<<" Two intersection case "<<std::endl;
//                             //     std::cout<<" First "<<pair.second.begin()->first.transpose()<<std::endl;
//                             //     std::cout<<" second "<<pair.second.rbegin()->first.transpose()<<std::endl;
//                             //     std::cout<<" difference "<<(pair.second.begin()->first-pair.second.rbegin()->first).norm()<<std::endl;
//                             //     std::cout<<" difference "<<((pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)<<std::endl;
//                             // }
                            
//                             if (pair.second.size()==2 && (pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)
//                             {
//                                 //The two intersection points are the same... Insert at diagonally opposite end

//                                 // std::cout<<"Coming in two intersection case "<<std::endl;

//                                 ppi.emplace_back(pair.second.begin()->first, pair.first->shift, std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID), nullptr);

//                                 const VectorDim localShift(pair.first->shift + pair.second.begin()->second->deltaShift+pair.second.rbegin()->second->deltaShift);
//                                 //here three patches are needed to be inserted... two at periodically opposite faces and one at diagonally opposite faces
//                                 lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.begin()->second->deltaShift);
//                                 currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
//                                 lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.rbegin()->second->deltaShift);
//                                 currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
//                                 lastPatch = patchBoundary.patches().getFromKey(localShift);
//                                 currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
//                                 assert(pair.second.begin()->second->twin); //Now the twin for both the edges must exist
//                                 assert(pair.second.rbegin()->second->twin);
                                
//                                 //Insert a new node at the twinned edge of the current pair of edges
//                                 ppi.emplace_back(pair.second.begin()->first, localShift, std::make_pair(pair.second.begin()->second->twin->edgeID, pair.second.rbegin()->second->twin->edgeID), nullptr);
//                             }
//                             else
//                             {
//                                 for (const auto &edge : pair.second)
//                                 { // loop over crossed edges
//                                     ppi.emplace_back(edge.first, pair.first->shift, std::make_pair(edge.second->edgeID,-1), nullptr);
//                                     //  shift += edge.second->deltaShift;
//                                     const VectorDim localShift(pair.first->shift + edge.second->deltaShift);
//                                     // Insert patches corresponding to the individual intersections as well
//                                     lastPatch = patchBoundary.patches().getFromKey(localShift);
//                                     currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
//                                     assert(edge.second->twin);
//                                     ppi.emplace_back(edge.first, edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1), nullptr);
//                                 }
//                             }

//                         }
//                         if (lastPatch->contains(endPoint.first))
//                         {
//                             ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
//                             //                        printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
//                             break;
//                         }
//                     }
//                 }
//                 // if (true)
//                 // {
//                 //     std::ofstream edgesFile("Debug/edges_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt", std::ofstream::out | std::ofstream::app);
//                 //     std::ofstream pointsFile("Debug/points_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt",std::ofstream::out | std::ofstream::app);

//                 //     for (const auto &edge : lastPatch->edges())
//                 //     {
//                 //         edgesFile << edge->source->sID << " " << edge->sink->sID << " " << lastPatch->sID << std::endl;
//                 //         pointsFile << "    " << edge->source->sID << " " << (*edge->source.get()).transpose() << std::endl;
//                 //     }
//                 //     edgesFile.close();
//                 //     pointsFile.close();
//                 // }
//             }

//             // if (true)
//             // {

//             //     std::ofstream polyFile("Debug/poly_"+std::to_string(runID) + "_" + std::to_string(loopID) + ".txt");
//             //     //                    std::ofstream poly3DFile("poly3D.txt");
//             //     for (const auto &tup : ppi)
//             //     {
//             //         //                        if(!node.second.expired())
//             //         //                        {
//             //         polyFile << "    " << std::get<0>(tup).transpose() << std::endl;
//             //         //                        }
//             //     }
//             //     polyFile.close();
//             // }
//         }


//         return ppi;
//     }

    //Function updated on December 14_2021 to work with the new version of update boundary nodes
template<int dim>
    template<typename T>
    std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>,std::map<typename PeriodicGlidePlane<dim>::VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>>,size_t, const T* const>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim, const T* const>>& polyPoints)
    {
        std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short  int>, std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> ,size_t,const T* const>> ppiPPedges; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
 // 2dPos,shift,std::pair<edgeIDs>,std::set<std::pair<edgeIDs>>,size_t,LoopNodeType* 
 //This is done to take care of all the edges that the two nodes are intersecting, size_t indicates the remaining intersection
        if (polyPoints.size())
        {
            PeriodicPatchBoundary<dim> patchBoundary(glidePlaneFactory, referencePlane);
            auto lastPatch(patchBoundary.getPatch(findPatch(polyPoints[0].first, VectorDim::Zero())));
            // const size_t loopID((*polyPoints.begin()).second->loop()->sID);
            // const size_t runID((*polyPoints.begin()).second->network().simulationParameters.runID);

            //                std::ofstream polyFile("poly.txt");
            //                for(const auto& node : polyPoints)
            //                {
            //                    polyFile<<"    "<<node.transpose()<<std::endl;
            //                }

            //        int removeMe=0;

            //        ppi.emplace_back(polyPoints[0],lastPatch,nullptr);

            // Find other points
            for (size_t k = 0; k < polyPoints.size(); ++k)
            {
                std::vector<std::tuple<typename PeriodicGlidePlane<dim>::VectorLowerDim,typename PeriodicGlidePlane<dim>::VectorDim,std::pair<short int,short int>,const T* const>> ppi; // 2dPos,shift,std::pair<edgeIDs>,LoopNodeType*
                std::map<VectorDim,std::set<std::pair<short int,short  int>>,CompareVectorsByComponent<double,dim,float>> edgeMaps; //This will populate the set of edges that this link has crossed along with the patchShift
                std::vector<std::pair<short int,short  int>> edgeVector; //This will just give the number of edges intersected
                std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> currentPatches;
                //            std::cout<<"k point ="<<k<<std::endl;
                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                const auto startPoint(polyPoints[k]);
                const auto endPoint(k == polyPoints.size() - 1 ? polyPoints[0] : polyPoints[k + 1]);
                while (true)
                {
                    std::map<const PeriodicPlanePatch<dim> *, std::vector<std::pair<const VectorLowerDim, const PeriodicPlaneEdge<dim> *const>>> crossdEdges;

                    for (const auto &bndEdge : patchBoundary.untwinnedEdges())
                    { // loop over outer boundaries and holes
                        SegmentSegmentDistance<dim - 1> ssd(startPoint.first, endPoint.first, *bndEdge->source, *bndEdge->sink);
                        // std::cout << " Determining intersection between " << startPoint.second->tag() << " and " << endPoint.second->tag()<<" dMin "<<ssd.dMin <<" den is "<< std::flush;
                        // const double alignemntAngle(((endPoint.first - startPoint.first).normalized()).dot((*bndEdge->sink - *bndEdge->source).normalized()));
                        // const double maxAngleAlignment(cos(1.0e-4 * M_PI / 180.0)); // Alignment angle is fixed to 10 degrees

                        if (ssd.isIntersecting(FLT_EPSILON)) 
                        { // intersection with current boundary found
                            // std::cout << std::setprecision(15) << " Intersection found " << ssd.t << " " << ssd.u << std::endl;
                            crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
                        }
                        //Changed to allow for the intersection of the near parallel segment with the boundary
                        // else
                        // {
                        //     if (startPoint.second && endPoint.second)
                        //     {
                        //         // Check the alignment of the links with the bndEdge Source
                        //         //                    	std::cout<<"Coming here to check for 3D intersection"<<std::endl;
                        //         const double maxAngleAlignment(cos(5.0 * M_PI / 180.0)); //Alignment angle is fixed to 10 degrees
                        //         // std::cout<<" Alignemnt Angle =>maxAngleAlignment"<<alignemntAngle<<" "<<maxAngleAlignment<<" for nodes "<<startPoint.second->tag()
                        //         // <<"=>"<<endPoint.second->tag()<<" ssdDmin "<<ssd.dMin<<std::endl;

                        //         if (fabs(alignemntAngle) > maxAngleAlignment && ssd.dMin < 1000*FLT_EPSILON) //If very close to the boundary and aligned only then check in 3D
                        //         {
                        //             //The value of 1000 is determined by trial and error
                        //             //Determine the patches of startPoint and endPoint
                        //             if (startPoint.second->periodicPlanePatch() != endPoint.second->periodicPlanePatch())
                        //             {
                        //                 //The two loop nodes belong to different patches, a node is needed to be inserted on the edge
                        //                 //This is only valid with the if for the alignment angle and the ssd.dMin
                        //                 crossdEdges[bndEdge->patch].emplace_back(0.5 * (ssd.x0 + ssd.x1), bndEdge);
                        //             }
                        //         }
                        //     }
                        // }
                        
                    }

                    if (crossdEdges.size() == 0)
                    { // No untwinned edges have been crossed, break and move to next pair of polyPoints
                        ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                        // edgeMaps.emplace(VectorDim::Zero(),std::make_pair(-1,-1));
                        //                    printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
                        break;
                    }
                    else
                    {
                        for (const auto &pair : crossdEdges)
                        { // loop over crossed patches
                            //VectorDim shift(pair.first->shift);
                            //                        std::cout<<"crossed patches="<<crossdEdges.size()<<std::endl;
                            assert((pair.second.size()>0 && pair.second.size()<=2) && "At max two periodic plane edges can be intersected (2D geometry constraint)");

                            // if (pair.second.size()==2)
                            // {
                            //     std::cout<<" Two intersection case "<<std::endl;
                            //     std::cout<<" First "<<pair.second.begin()->first.transpose()<<std::endl;
                            //     std::cout<<" second "<<pair.second.rbegin()->first.transpose()<<std::endl;
                            //     std::cout<<" difference "<<(pair.second.begin()->first-pair.second.rbegin()->first).norm()<<std::endl;
                            //     std::cout<<" difference "<<((pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)<<std::endl;
                            // }
                            
                            if (pair.second.size()==2 && (pair.second.begin()->first-pair.second.rbegin()->first).norm()<FLT_EPSILON)
                            {
                                // pair.second.size() cannot be 4 even for diagonally opposite ends since we are looping over the untwineed edges for the intersection.
                                //Need to check that two subsequen links are at the same 2D position

                                //The two intersection points are the same... Insert at diagonally opposite end

                                // std::cout<<"Coming in two intersection case "<<std::endl;
                                assert((pair.second.begin()->second!=pair.second.rbegin()->second) && "Two intersection must be different");
                                ppi.emplace_back(pair.second.begin()->first, pair.first->shift, std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID), nullptr);
                                const auto edgeSetIter1 (edgeMaps.find(pair.first->shift));
                                edgeVector.push_back(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                
                                if (edgeSetIter1==edgeMaps.end())
                                {
                                    std::set<std::pair<short int,short  int>> edgeSet;
                                    edgeSet.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                    // edgeMaps.emplace(std::piecewise_construct,
                                    //                 std::forward_as_tuple(pair.first->shift),
                                    //                 std::forward_as_tuple(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                    edgeMaps.emplace(pair.first->shift,edgeSet);
                                }
                                else
                                {
                                    edgeSetIter1->second.insert(std::make_pair(pair.second.begin()->second->edgeID, pair.second.rbegin()->second->edgeID));
                                }

                                const VectorDim localShift(pair.first->shift + pair.second.begin()->second->deltaShift+pair.second.rbegin()->second->deltaShift);
                                //here three patches are needed to be inserted... two at periodically opposite faces and one at diagonally opposite faces
                                lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.begin()->second->deltaShift);
                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                lastPatch = patchBoundary.patches().getFromKey(pair.first->shift + pair.second.rbegin()->second->deltaShift);
                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                lastPatch = patchBoundary.patches().getFromKey(localShift);
                                currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                assert(pair.second.begin()->second->twin); //Now the twin for both the edges must exist
                                assert(pair.second.rbegin()->second->twin);
                                assert((pair.second.begin()->second->twin!=pair.second.rbegin()->second->twin) && "Two intersection must be different after twinned link");
                                //Insert a new node at the twinned edge of the current pair of edges
                                //Need to determine the correct edgeIDs
                                short int beginTwinEdgeID(-1);
                                short rbeginTwinEdgeID(-1);
                                //We need to get the edgeID corresponding to the diagonally opposite patch
                                if (pair.second.begin()->second==pair.second.rbegin()->second->prev)
                                {
                                    assert(pair.second.begin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");
                                    assert(pair.second.rbegin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                    beginTwinEdgeID=pair.second.begin()->second->twin->prev->twin->edgeID;
                                    rbeginTwinEdgeID=pair.second.rbegin()->second->twin->next->twin->edgeID;
                                }
                                else if (pair.second.begin()->second==pair.second.rbegin()->second->next)
                                {
                                    assert(pair.second.begin()->second->twin->next->twin && "A twin in the diagonally opposite patch must exist");
                                    assert(pair.second.rbegin()->second->twin->prev->twin && "A twin in the diagonally opposite patch must exist");

                                    beginTwinEdgeID=pair.second.begin()->second->twin->next->twin->edgeID;
                                    rbeginTwinEdgeID=pair.second.rbegin()->second->twin->prev->twin->edgeID;
                                }
                                else
                                {
                                    assert(false && "Edges must be consecutive");
                                }
                                
                                assert(beginTwinEdgeID!=-1 && "A twin ID must exist");
                                assert(rbeginTwinEdgeID!=-1 && "A twin ID must exist");

                                ppi.emplace_back(pair.second.begin()->first, localShift, std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID), nullptr);
                                edgeVector.push_back(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));

                                const auto edgeSetIter2 (edgeMaps.find(localShift));
                                if (edgeSetIter2==edgeMaps.end())
                                {
                                    std::set<std::pair<short int,short  int>> edgeSet;
                                    edgeSet.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                    // edgeMaps.emplace(std::piecewise_construct,
                                    //                 std::forward_as_tuple(localShift),
                                    //                 std::forward_as_tuple(pair.second.begin()->second->twin->edgeID, pair.second.rbegin()->second->twin->edgeID));
                                    edgeMaps.emplace(localShift,edgeSet);
                                }
                                else
                                {
                                    edgeSetIter2->second.insert(std::make_pair(beginTwinEdgeID, rbeginTwinEdgeID));
                                }
                            }
                            else
                            {
                                for (const auto &edge : pair.second)
                                { // loop over crossed edges
                                    ppi.emplace_back(edge.first, pair.first->shift, std::make_pair(edge.second->edgeID,-1), nullptr);
                                    edgeVector.push_back(std::make_pair(edge.second->edgeID,-1));

                                    const auto edgeSetIter1(edgeMaps.find(pair.first->shift));
                                    if (edgeSetIter1 == edgeMaps.end())
                                    {
                                        std::set<std::pair<short int,short  int>> edgeSet;
                                        edgeSet.insert(std::make_pair(edge.second->edgeID,-1));
                                        // edgeMaps.emplace(std::piecewise_construct,
                                        //             std::forward_as_tuple(pair.first->shift),
                                        //             std::forward_as_tuple(edge.second->edgeID,-1));
                                        edgeMaps.emplace(pair.first->shift, edgeSet);
                                    }
                                    else
                                    {
                                        edgeSetIter1->second.insert(std::make_pair(edge.second->edgeID,-1));
                                    }
                                    //  shift += edge.second->deltaShift;
                                    const VectorDim localShift(pair.first->shift + edge.second->deltaShift);
                                    // Insert patches corresponding to the individual intersections as well
                                    lastPatch = patchBoundary.patches().getFromKey(localShift);
                                    currentPatches.push_back(lastPatch); // keep patches alive durin current line segment search
                                    assert(edge.second->twin);
                                    ppi.emplace_back(edge.first, edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1), nullptr);
                                    edgeVector.push_back(std::make_pair(edge.second->twin->edgeID,-1));

                                    const auto edgeSetIter2(edgeMaps.find(edge.second->twin->patch->shift));
                                    if (edgeSetIter2 == edgeMaps.end())
                                    {
                                        std::set<std::pair<short int,short  int>> edgeSet;
                                        edgeSet.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                        // edgeMaps.emplace(std::piecewise_construct,
                                        //             std::forward_as_tuple(edge.second->twin->patch->shift),
                                        //             std::forward_as_tuple(edge.second->twin->edgeID,-1));
                                        // edgeMaps.emplace(edge.second->twin->patch->shift, std::make_pair(edge.second->twin->edgeID,-1));
                                        edgeMaps.emplace(edge.second->twin->patch->shift, edgeSet);
                                    }
                                    else
                                    {
                                        edgeSetIter2->second.insert(std::make_pair(edge.second->twin->edgeID,-1));
                                    }
                                }
                            }
                        }
                        if (lastPatch->contains(endPoint.first))
                        {
                            ppi.emplace_back(endPoint.first, lastPatch->shift, std::make_pair(-1,-1), endPoint.second);
                            // edgeMaps.emplace(VectorDim::Zero(), std::make_pair(-1, -1));
                            //                        printIntersections(ppi,currentPatches,lastPatch,removeMe,patchBoundary.untwinnedEdges());
                            break;
                        }
                    }
                }
                //Here populate the ppiPPedges
                size_t crossCount=0;
                for (const auto& ppiIter : ppi )
                {
                    // //LoopNode* is the last
                    // if (edgeVector.size()==0)
                    // {
                    //     ppiPPedges.emplace_back(std::get<0>(ppiIter),std::get<1>(ppiIter),std::get<2>(ppiIter),edgeMaps,0,std::get<3>(ppiIter));
                    // }
                    // else
                    // {
                    size_t edgesStillRemaining(0);
                    if (std::get<3>(ppiIter) == nullptr)
                    {
                        //That mean intersection with the edge is here
                        crossCount++;
                        edgesStillRemaining = edgeVector.size() - crossCount; 
                    }
                    ppiPPedges.emplace_back(std::get<0>(ppiIter), std::get<1>(ppiIter), std::get<2>(ppiIter), edgeMaps, edgesStillRemaining, std::get<3>(ppiIter));
                    // }
                }
                if (edgeVector.size() != 0)
                {
                    if ((edgeVector.size() - crossCount ) != 0)  
                    {
                        std::cout<<"Edge VectorSize "<<edgeVector.size()<<std::endl;
                        std::cout<<"Cross count "<<crossCount<<std::endl;
                    }
                    assert((edgeVector.size() - crossCount) == 0 && "No edges should be remaining after completion of polygonPatchIntersection");
                }
                // if (true)
                // {
                //     std::ofstream edgesFile("Debug/edges_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt", std::ofstream::out | std::ofstream::app);
                //     std::ofstream pointsFile("Debug/points_"+ std::to_string(runID) + "_" + std::to_string(loopID) + ".txt",std::ofstream::out | std::ofstream::app);

                //     for (const auto &edge : lastPatch->edges())
                //     {
                //         edgesFile << edge->source->sID << " " << edge->sink->sID << " " << lastPatch->sID << std::endl;
                //         pointsFile << "    " << edge->source->sID << " " << (*edge->source.get()).transpose() << std::endl;
                //     }
                //     edgesFile.close();
                //     pointsFile.close();
                // }
            }

            // if (true)
            // {

            //     std::ofstream polyFile("Debug/poly_"+std::to_string(runID) + "_" + std::to_string(loopID) + ".txt");
            //     //                    std::ofstream poly3DFile("poly3D.txt");
            //     for (const auto &tup : ppi)
            //     {
            //         //                        if(!node.second.expired())
            //         //                        {
            //         polyFile << "    " << std::get<0>(tup).transpose() << std::endl;
            //         //                        }
            //     }
            //     polyFile.close();
            // }
        }


        return ppiPPedges;
    }

    template<typename T1,typename T2, typename T3,typename T4>
    void printIntersections(const T1& ppi,const T2& usedPatches,const T3& lastPatch,int& removeMe, const T4& untwinnedEdges)
    {
        std::ofstream pointsFile("points_"+std::to_string(removeMe)+".txt");
        for(const auto& tup : ppi)
        {
            //                        if(!node.second.expired())
            //                        {
            pointsFile<<"    "<<std::get<0>(tup).transpose()<<std::endl;
            //                        }
        }
        
        std::ofstream edgesFile("edges_"+std::to_string(removeMe)+".txt");
        for(const auto& patch : usedPatches)
        {
            for(const auto& edge : patch->edges())
            {
                if(untwinnedEdges.find(edge.get())!=untwinnedEdges.end())
                {
                    edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 1"<<std::endl;
                }
                else
                {
                    edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 0"<<std::endl;
                }
                    
            }
        }
        
        std::ofstream lastPatchFile("lastPatch_"+std::to_string(removeMe)+".txt");
        for(const auto& edge : lastPatch->edges())
        {
            lastPatchFile<<lastPatch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<std::endl;
        }
        removeMe++;
    }
    
//    template<int dim>
//    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorLowerDim>& polyPoints)
//    {
//
//        if(polyPoints.size()>=3)
//        {
//
//            //                std::set<long int> pointsHeights;
//            //                for(const auto& point : polyPoints)
//            //                {
//            //                    const auto heightPair(LatticePlane::computeHeight(this->referencePlane->n,point));
//            //                    assert(heightPair.first && "Point not on a lattice plane");
//            //                    pointsHeights.insert(heightPair.second);
//            //                }
//            //                assert(pointsHeights.size()==1 && "polyPoints on different planes");
//            ////                const GlidePlaneKey<dim> pointsPlaneKey(this->referencePlane->grain.grainID,polyPoints[0],this->referencePlane->n);
//            //                const GlidePlaneKey<dim> pointsPlaneKey(polyPoints[0],this->referencePlane->n);
//
//            //                const auto pointsPlane(this->glidePlaneFactory.get(pointsPlaneKey));
//            //                const VectorDim pointsShift(pointsPlane->P-this->referencePlane->P);
//            getPatch(VectorDim::Zero());
//
//            VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
//            for(const auto& seg : this->referencePlane->meshIntersections) // problem is here for referencePlane not cutting mesh
//            {
//                insideReferencePoint+=this->getLocalPosition(seg->P0,VectorDim::Zero());
//            }
//            insideReferencePoint/=this->referencePlane->meshIntersections.size();
//            assert(this->isInsideOuterBoundary(insideReferencePoint));
//
//
//            //                std::deque<std::shared_ptr<PeriodicPlanePatch<dim>>> tempPatches;
//
//            const VectorLowerDim P0(polyPoints[0]);
//
//            while(!this->isInsideOuterBoundary(P0))
//            {
//                std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
//                for(const auto& bndEdge : this->untwinnedEdges())
//                {// loop over outer boundaries and holes
//                    SegmentSegmentDistance<dim-1> ssd(P0,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
//                    if(ssd.dMin<FLT_EPSILON)
//                    {// intersection with current boundary found
//                        crossdEdges.insert(bndEdge);
//                    }
//                }
//
//                switch (crossdEdges.size())
//                {
//                    case 1:
//                    {
//                        getPatch(getShift(**crossdEdges.begin())+(*crossdEdges.begin())->patch->shift);
//                        break;
//                    }
//
//                    case 2:
//                    {
//                        assert((*crossdEdges.begin())->patch==(*crossdEdges.rbegin())->patch);
//                        getPatch(getShift(**crossdEdges.begin())+getShift(**crossdEdges.rbegin())+(*crossdEdges.begin())->patch->shift);
//                        break;
//                    }
//
//                    default:
//                    {
//                        std::cout<<"crossdEdges.size()="<<crossdEdges.size()<<std::endl;
//                        assert(false && "1 or 2 edges must be crossed");
//                        break;
//                    }
//                }
//
//                assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
//            }
//
//            // erase patches needed to find P0
//            std::list<VectorDim> eraseKeys;
//            for(const auto& patch : patches())
//            {
//                if(!patch.second->contains(P0))
//                {
//                    eraseKeys.push_back(patch.first);
//                }
//            }
//            for(const auto& key : eraseKeys)
//            {
//                patches().erase(key);
//            }
//
//            for(size_t k=0;k<polyPoints.size();++k)
//            {
//                const VectorLowerDim startPoint(polyPoints[k]);
//                const VectorLowerDim endPoint(k==polyPoints.size()-1? polyPoints[0] : polyPoints[k+1]);
//                while(true)
//                {
//                    //                        std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
//                    std::map<const PeriodicPlanePatch<dim>*,std::set<const PeriodicPlaneEdge<dim>*>> crossdEdges;
//                    for(const auto& bndEdge : this->untwinnedEdges())
//                    {// loop over outer boundaries and holes
//                        SegmentSegmentDistance<dim-1> ssd(startPoint,endPoint,*bndEdge->source,*bndEdge->sink);
//                        if(ssd.dMin<FLT_EPSILON)
//                        {// intersection with current boundary found
//                            crossdEdges[bndEdge->patch].insert(bndEdge);
//                        }
//                    }
//
//                    if(crossdEdges.size()==0)
//                    {// No untwinned edges have been crossed, break and move to next pair of polyPoints
//                        break;
//                    }
//                    else
//                    {
//                        for(const auto& pair : crossdEdges)
//                        {
//                            VectorDim shift(pair.first->shift);
//                            for(const auto& edge : pair.second)
//                            {
//                                shift+=getShift(*edge);
//                            }
//                            getPatch(shift); // this will change untwinnedEdges. Remain in while loop
//
//                        }
//                    }
//                    assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
//                }
//            }
//
//            fillHoles();
//
//            if(true)
//            {
//
//
//                std::ofstream polyFile("poly.txt");
//                //                    std::ofstream poly3DFile("poly3D.txt");
//
//                polyFile<<insideReferencePoint.transpose()<<std::endl;
//                for(const auto& node : polyPoints)
//                {
//                    polyFile<<"    "<<node.transpose()<<std::endl;
//                    //                        poly3DFile<<node.transpose()<<std::endl;
//                }
//
//                std::ofstream pointsFile("points.txt");
//                for(const auto& node : this->nodes())
//                {
//                    if(!node.second.expired())
//                    {
//                        pointsFile<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
//                    }
//                }
//
//                std::ofstream edgesFile("edges.txt");
//                for(const auto& patch : patches())
//                {
//                    for(const auto& edge : patch.second->edges())
//                    {
//                        edgesFile<<edge->tag()<<std::endl;
//                    }
//                }
//                assert(this->isCompact() && "Plane not compact");
//            }
//        }
//    }
//
//
//    template<int dim>
//    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints,const VectorDim& shift)
//    {
//        std::vector<VectorLowerDim> lowerPolyPoints;
//        for(const auto& point : polyPoints)
//        {
//            assert(this->referencePlane->contains(point-shift) && "reference plane does not cointain point");
//            lowerPolyPoints.push_back(this->getLocalPosition(point,shift));
//        }
//        addPatchesContainingPolygon(lowerPolyPoints);
//    }
//
//    template<int dim>
//    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints)
//    {
//        std::vector<VectorLowerDim> lowerPolyPoints;
//        for(const auto& point : polyPoints)
//        {
//            assert(this->referencePlane->contains(point) && "reference plane does not cointain point");
//            lowerPolyPoints.push_back(this->getLocalPosition(point,VectorDim::Zero()));
//        }
//        addPatchesContainingPolygon(lowerPolyPoints);
//    }
    
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
        
//        std::cout<<"PeriodiPlane patches:"<<std::endl;
//        for(const auto& patch : patches())
//        {
//            std::cout<<"Patch "<<patch.shift.transpose()<<std::endl;
//            for(const auto& edge : patch.second.edges())
//            {
//                std::cout<<edge->source->sID<<"->"<<edge->sink->sID<<std::endl;
//            }
//        }
        
        
    }
    


    
//    template<int dim>
//    std::shared_ptr<PeriodicPlanePatch<dim>> PeriodicGlidePlane<dim>::getPatch(const VectorDim& shift)
//    {
//        return patches().get(shift);
//    }
    
    template<int dim>
    typename PeriodicGlidePlane<dim>::VectorDim PeriodicGlidePlane<dim>::getGlidePlaneShiftfromReferencePlane(const GlidePlane<dim> *gp) const
    {
        const long int t(gp->key.planeIndex()-this->referencePlane->key.planeIndex());
        const VectorDimI alphas(gp->key.reciprocalDirectionComponents().transpose()*periodicGlidePlaneFactory.N);
        VectorDimI sol(VectorDimI::Zero());
        solveDiophantine3vars(alphas,t,sol);
        //            VectorDim shift(periodicGlidePlaneFactory.B*sol.template cast<double>());
        //            return shift;
        return periodicGlidePlaneFactory.B*sol.template cast<double>();
    }
    
//    template<int dim>
//    GlidePlaneKey<dim> PeriodicGlidePlane<dim>::getGlidePlaneKey(const VectorDim& shift)
//    {
//        return GlidePlaneKey<dim>(referencePlane->P+shift,referencePlane->n);
//    }
//
//    template<int dim>
//    std::shared_ptr<GlidePlane<dim>> PeriodicGlidePlane<dim>::getGlidePlane(const VectorDim& shift)
//    {
//        return glidePlaneFactory.getFromKey(getGlidePlaneKey(shift));
//    }
    
    /**********************************************************************/
        /**********************************************************************/
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
    
    template <int dim,class T>
    T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
    {
        os  << ds.glidePlaneID<<"\t"
        /**/<< ds.patchBoundaryID<<"\t"
        /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
        return os;
    }
    
}
#endif
