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




namespace model
{

    template<int dim>
    PeriodicPlaneEdge<dim>::PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                                              const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in) :
    /* init */ patch(patch_in)
    /* init */,source(source_in)
    /* init */,sink(sink_in)
    /* init */,meshIntersection(meshIntersection_in)
    /* init */,next(nullptr)
    /* init */,prev(nullptr)
    /* init */,twin(nullptr)
    {
        //        std::cout<<"Creating PeriodicPlaneEdge "<<this->tag()<<std::endl;
        source->addLink(this);
        sink->addLink(this);
        if(twin)
        {
            patch->periodicPlane->removeUntwinnedEdge(twin);
        }
        else
        {
            patch->periodicPlane->addUntwinnedEdge(this);
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
            patch->periodicPlane->addUntwinnedEdge(twin);
        }
        else
        {
            patch->periodicPlane->removeUntwinnedEdge(this);
        }
    }
    
    template<int dim>
    std::string PeriodicPlaneEdge<dim>::tag() const
    {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
      */
        return std::to_string(source->sID) + " " + std::to_string(sink->sID)+" "+std::to_string(patch->sID);
    }
    
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlaneNode<dim>::PeriodicPlaneNode(const VectorLowerDim& pos) :
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
    PeriodicPlanePatch<dim>::PeriodicPlanePatch(PeriodicGlidePlane<dim>& periodicPlane_in,
                       const VectorDim& shift_in) :
    /* init */ periodicPlane(&periodicPlane_in)
    /* init */,shift(shift_in)
    /* init */,glidePlane(periodicPlane->getGlidePlane(shift))
    {
        //            std::cout<<cyanColor<<"Creating patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
        
        if(isRightHandedBoundary(glidePlane->meshIntersections,*periodicPlane->referencePlane))
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
        
        periodicPlane->updateBoundaries();
        
    }
    
    template<int dim>
    PeriodicPlanePatch<dim>::~PeriodicPlanePatch()
    {
        //            std::cout<<cyanColor<<"Destroying patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
        this->clear(); // call before updateBoundaries
        periodicPlane->updateBoundaries();
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
        for(size_t k=0;k<bms.size();++k)
        {// compute 2d points on PeriodicPlane, and coneect them
            std::shared_ptr<PeriodicPlaneNode<dim>> source(periodicPlane->getSharedNode(bms[k]->P0,shift));
            std::shared_ptr<PeriodicPlaneNode<dim>>   sink(periodicPlane->getSharedNode(bms[k]->P1,shift));
            this->emplace_back(new PeriodicPlaneEdge<dim>(this,source,sink,bms[k]));
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
    
    /**********************************************************************/
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::MatrixDim PeriodicGlidePlaneBase<dim>::getL2G(VectorDim z)
    {
        //            const double xNorm(x.norm());
        const double zNorm(z.norm());
        assert(zNorm>FLT_EPSILON);
        z/=zNorm;
        
        VectorDim x(VectorDim::UnitX().cross(z));
        double xNorm(x.norm());
        if(xNorm>FLT_EPSILON)
        {
            x=x/xNorm;
        }
        else
        {
            x=VectorDim::UnitY().cross(z);
            xNorm=x.norm();
            if(xNorm>FLT_EPSILON)
            {
                x=x/xNorm;
            }
            else
            {
                x=VectorDim::UnitZ().cross(z);
                xNorm=x.norm();
                if(xNorm>FLT_EPSILON)
                {
                    x=x/xNorm;
                }
                else
                {
                    assert(false && "CANNOT FIND VECTOR ORTHOGONAL TO z");
                }
            }
        }
        
        assert(std::fabs(x.norm()-1.0)<FLT_EPSILON);
        assert(std::fabs(z.norm()-1.0)<FLT_EPSILON);
        assert(fabs(x.dot(z)<FLT_EPSILON));
        MatrixDim temp(Eigen::Matrix3d::Identity());
        temp.col(2)=z;
        temp.col(0)=x;
        temp.col(1)=temp.col(2).cross(temp.col(0));
        return temp;
    }
    
    template<int dim>
    PeriodicGlidePlaneBase<dim>::PeriodicGlidePlaneBase(GlidePlaneFactory<dim>& glidePlaneFactory_in,
                           const GlidePlaneKey<dim>& referencePlaneKey) :
    /* init */ glidePlaneFactory(glidePlaneFactory_in)
    /* init */,referencePlane(glidePlaneFactory.get(referencePlaneKey))
    /* init */,L2G(getL2G(referencePlane->unitNormal))
    {
        
        assert(referencePlane->meshIntersections.size()>=3);
    }
    
    //        ~PeriodicGlidePlaneBase()
    //        {
    //            std::cout<<"DESTROYING PeriodicGlidePlaneBase"<<std::endl;
    //        }
    
    
    
    
    
    
    template<int dim>
    GlidePlaneKey<dim> PeriodicGlidePlaneBase<dim>::getGlidePlaneKey(const VectorDim& shift)
    {
        return GlidePlaneKey<dim>(referencePlane->P+shift,referencePlane->n);
    }
    
    template<int dim>
    std::shared_ptr<GlidePlane<dim>> PeriodicGlidePlaneBase<dim>::getGlidePlane(const VectorDim& shift)
    {
        return glidePlaneFactory.get(getGlidePlaneKey(shift));
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::BoundariesContainerType& PeriodicGlidePlaneBase<dim>::outerBoundaries()
    {
        return _outerBoundaries;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneBase<dim>::BoundariesContainerType& PeriodicGlidePlaneBase<dim>::outerBoundaries() const
    {
        return _outerBoundaries;
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::BoundariesContainerType& PeriodicGlidePlaneBase<dim>::innerBoundaries()
    {
        return _innerBoundaries;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneBase<dim>::BoundariesContainerType& PeriodicGlidePlaneBase<dim>::innerBoundaries() const
    {
        return _innerBoundaries;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneBase<dim>::UntwinnedEdgeContainerType& PeriodicGlidePlaneBase<dim>::untwinnedEdges() const
    {
        return *this;
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::UntwinnedEdgeContainerType& PeriodicGlidePlaneBase<dim>::untwinnedEdges()
    {
        return *this;
    }
    
    template<int dim>
    void PeriodicGlidePlaneBase<dim>::addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
    {
        untwinnedEdges().insert(link);
    }
    
    template<int dim>
    void PeriodicGlidePlaneBase<dim>::removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link)
    {
        const size_t erased(untwinnedEdges().erase(link));
        assert(erased==1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::NodeCointainerType& PeriodicGlidePlaneBase<dim>::nodes()
    {
        return *this;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneBase<dim>::NodeCointainerType& PeriodicGlidePlaneBase<dim>::nodes() const
    {
        return *this;
    }
    
    template<int dim>
    int PeriodicGlidePlaneBase<dim>::isInsideOuterBoundary(const VectorLowerDim& test)
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
    bool PeriodicGlidePlaneBase<dim>::isCompact() const
    {
        return outerBoundaries().size()==1 && innerBoundaries().size()==0;
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::VectorDim PeriodicGlidePlaneBase<dim>::rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
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
    double PeriodicGlidePlaneBase<dim>::rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
    {
        return rightHandedNormal(bnd).dot(referencePlane->unitNormal);
    }
    
    template<int dim>
    bool PeriodicGlidePlaneBase<dim>::isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd)
    {
        return rightHandedArea(bnd)>=0.0;
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::VectorLowerDim PeriodicGlidePlaneBase<dim>::getLocalPosition(const VectorDim& point,const VectorDim& shift) const
    {
        const VectorDim pointLocal(L2G.transpose()*(point-shift-referencePlane->P));
        if(fabs(pointLocal(2))>FLT_EPSILON)
        {
            std::cout<<"point"<<point.transpose()<<std::endl;
            std::cout<<"shift"<<shift.transpose()<<std::endl;
            std::cout<<"referencePlane->P"<<referencePlane->P.transpose()<<std::endl;
            std::cout<<"L2G=\n"<<L2G<<std::endl;
            std::cout<<"pointLocal"<<pointLocal.transpose()<<std::endl;
            assert(false && "local point has non-zero z-coordinate");
        }
        return pointLocal.template segment<2>(0);
    }
    
    template<int dim>
    typename PeriodicGlidePlaneBase<dim>::VectorDim PeriodicGlidePlaneBase<dim>::getGlobalPosition(const VectorLowerDim& point) const
    {// terurns the position on the plane in global goordinates
        return L2G.template block<dim,dim-1>(0,0)*point+referencePlane->P;
    }
    
    template<int dim>
    std::shared_ptr<PeriodicPlaneNode<dim> > PeriodicGlidePlaneBase<dim>::getSharedNode(const VectorDim& pointDim,const VectorDim& shift)
    {
        const VectorLowerDim point(getLocalPosition(pointDim,shift));
        const auto iter(nodes().find(point));
        if(iter==nodes().end())
        {// point does not exist
            //                std::shared_ptr<PeriodicPlaneNode<dim>> newNode(new PeriodicPlaneNode<dim>(point));
            //                nodes().emplace(point,newNode.get());
            return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
            //                return newNode;
        }
        else
        {// point exists
            if(iter->second.expired())
            {// node deleted elsewhere
                nodes().erase(iter);
                return nodes().emplace(point,std::shared_ptr<PeriodicPlaneNode<dim> >(new PeriodicPlaneNode<dim>(point))).first->second.lock();
            }
            else
            {
                return iter->second.lock();
            }
        }
    }


    /**********************************************************************/
    template<int dim>
    PeriodicGlidePlane<dim>::PeriodicGlidePlane(PeriodicGlidePlaneFactoryType& pgpf,
                       const GlidePlaneKeyType& key_in) :
    /* init */ PeriodicGlidePlaneBase<dim>(pgpf.glidePlaneFactory,key_in)
    /* init */,periodicGlidePlaneFactory(pgpf)
    {
        
    }
    
    template<int dim>
    const typename PeriodicGlidePlane<dim>::PatchContainerType& PeriodicGlidePlane<dim>::patches() const
    {
        return *this;
    }
    
    template<int dim>
    typename PeriodicGlidePlane<dim>::PatchContainerType& PeriodicGlidePlane<dim>::patches()
    {
        return *this;
    }
    
    template<int dim>
    void PeriodicGlidePlane<dim>::createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy)
    {
        //            std::cout<<"createNewBoundary"<<std::endl;
        BoundaryContainerType temp;
        temp.reserve(untwinnedCopy.size());
        while(true)
        {
            //                std::cout<<"currentEdge "<<currentEdge->tag()<<std::endl;
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
                
                
                std::ofstream edgesFile("edges.txt");
                for(const auto& patch : patches())
                {
                    for(const auto& edge : patch.second->edges())
                    {
                        edgesFile<<edge->tag()<<std::endl;
                    }
                }
                
                
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
    
    template<int dim>
    void PeriodicGlidePlane<dim>::updateBoundaries()
    {
        //            std::cout<<"Updating OuterBoundary"<<std::endl;
        this->outerBoundaries().clear();
        this->innerBoundaries().clear();
        UntwinnedEdgeContainerType untwinnedCopy(this->untwinnedEdges());
        while(untwinnedCopy.size())
        {
            for(const auto& edge : untwinnedCopy)
            {
                if(edge->source->outEdges().size()==1 && edge->source->inEdges().size()==1)
                {// must start from node with only one edge in and one edge out
                    createNewBoundary(edge,untwinnedCopy);
                    break;
                }
            }
        }
    }
    
    template<int dim>
    typename PeriodicGlidePlane<dim>::VectorDim PeriodicGlidePlane<dim>::getShift(const PeriodicPlaneEdge<dim>& edge) const
    {
        
        VectorDim shift(VectorDim::Zero());
        for(const PlanarMeshFace<dim>* const face : edge.meshIntersection->faces)
        {
            const auto parallelFaceID(edge.patch->glidePlane->grain.region.parallelFaces().at(face->sID));
            const auto parallelFace(edge.patch->glidePlane->grain.region.faces().at(parallelFaceID));
            shift+=(parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal();
            
        }
        return shift;
    }
    
    template<int dim>
    void PeriodicGlidePlane<dim>::fillHoles()
    {
        while(this->innerBoundaries().size())
        {
            const PeriodicPlaneEdge<dim>* const holeEdge(*this->innerBoundaries().front().begin());
            getPatch(getShift(*holeEdge)+holeEdge->patch->shift);
        }
    }
    
    template<int dim>
    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints,const VectorDim& shift)
    {
        std::vector<VectorLowerDim> lowerPolyPoints;
        for(const auto& point : polyPoints)
        {
            assert(this->referencePlane->contains(point-shift) && "reference plane does not cointain point");
            lowerPolyPoints.push_back(this->getLocalPosition(point,shift));
        }
        addPatchesContainingPolygon(lowerPolyPoints);
    }
    
    template<int dim>
    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorDim>& polyPoints)
    {
        std::vector<VectorLowerDim> lowerPolyPoints;
        for(const auto& point : polyPoints)
        {
            assert(this->referencePlane->contains(point) && "reference plane does not cointain point");
            lowerPolyPoints.push_back(this->getLocalPosition(point,VectorDim::Zero()));
        }
        addPatchesContainingPolygon(lowerPolyPoints);
    }
    
    template<int dim>
    void PeriodicGlidePlane<dim>::addPatchesContainingPolygon(const std::vector<VectorLowerDim>& polyPoints)
    {
        
        if(polyPoints.size()>=3)
        {
            
            //                std::set<long int> pointsHeights;
            //                for(const auto& point : polyPoints)
            //                {
            //                    const auto heightPair(LatticePlane::computeHeight(this->referencePlane->n,point));
            //                    assert(heightPair.first && "Point not on a lattice plane");
            //                    pointsHeights.insert(heightPair.second);
            //                }
            //                assert(pointsHeights.size()==1 && "polyPoints on different planes");
            ////                const GlidePlaneKey<dim> pointsPlaneKey(this->referencePlane->grain.grainID,polyPoints[0],this->referencePlane->n);
            //                const GlidePlaneKey<dim> pointsPlaneKey(polyPoints[0],this->referencePlane->n);
            
            //                const auto pointsPlane(this->glidePlaneFactory.get(pointsPlaneKey));
            //                const VectorDim pointsShift(pointsPlane->P-this->referencePlane->P);
            getPatch(VectorDim::Zero());
            
            VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
            for(const auto& seg : this->referencePlane->meshIntersections) // problem is here for referencePlane not cutting mesh
            {
                insideReferencePoint+=this->getLocalPosition(seg->P0,VectorDim::Zero());
            }
            insideReferencePoint/=this->referencePlane->meshIntersections.size();
            assert(this->isInsideOuterBoundary(insideReferencePoint));
            
            
            //                std::deque<std::shared_ptr<PeriodicPlanePatch<dim>>> tempPatches;
            
            const VectorLowerDim P0(polyPoints[0]);
            
            while(!this->isInsideOuterBoundary(P0))
            {
                std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
                for(const auto& bndEdge : this->untwinnedEdges())
                {// loop over outer boundaries and holes
                    SegmentSegmentDistance<dim-1> ssd(P0,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
                    if(ssd.dMin<FLT_EPSILON)
                    {// intersection with current boundary found
                        crossdEdges.insert(bndEdge);
                    }
                }
                
                switch (crossdEdges.size())
                {
                    case 1:
                    {
                        getPatch(getShift(**crossdEdges.begin())+(*crossdEdges.begin())->patch->shift);
                        break;
                    }
                        
                    case 2:
                    {
                        assert((*crossdEdges.begin())->patch==(*crossdEdges.rbegin())->patch);
                        getPatch(getShift(**crossdEdges.begin())+getShift(**crossdEdges.rbegin())+(*crossdEdges.begin())->patch->shift);
                        break;
                    }
                        
                    default:
                    {
                        std::cout<<"crossdEdges.size()="<<crossdEdges.size()<<std::endl;
                        assert(false && "1 or 2 edges must be crossed");
                        break;
                    }
                }
                
                assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
            }
            
            // erase patches needed to find P0
            std::list<VectorDim> eraseKeys;
            for(const auto& patch : patches())
            {
                if(!patch.second->contains(P0))
                {
                    eraseKeys.push_back(patch.first);
                }
            }
            for(const auto& key : eraseKeys)
            {
                patches().erase(key);
            }
            
            for(size_t k=0;k<polyPoints.size();++k)
            {
                const VectorLowerDim startPoint(polyPoints[k]);
                const VectorLowerDim endPoint(k==polyPoints.size()-1? polyPoints[0] : polyPoints[k+1]);
                while(true)
                {
                    //                        std::set<const PeriodicPlaneEdge<dim>*> crossdEdges;
                    std::map<const PeriodicPlanePatch<dim>*,std::set<const PeriodicPlaneEdge<dim>*>> crossdEdges;
                    for(const auto& bndEdge : this->untwinnedEdges())
                    {// loop over outer boundaries and holes
                        SegmentSegmentDistance<dim-1> ssd(startPoint,endPoint,*bndEdge->source,*bndEdge->sink);
                        if(ssd.dMin<FLT_EPSILON)
                        {// intersection with current boundary found
                            crossdEdges[bndEdge->patch].insert(bndEdge);
                        }
                    }
                    
                    if(crossdEdges.size()==0)
                    {// No untwinned edges have been crossed, break and move to next pair of polyPoints
                        break;
                    }
                    else
                    {
                        for(const auto& pair : crossdEdges)
                        {
                            VectorDim shift(pair.first->shift);
                            for(const auto& edge : pair.second)
                            {
                                shift+=getShift(*edge);
                            }
                            getPatch(shift); // this will change untwinnedEdges. Remain in while loop
                            
                        }
                    }
                    assert(this->outerBoundaries().size()==1 && "THERE MUST BE ONLY ONE OUTER BOUNDARY");
                }
            }
            
            fillHoles();
            
            if(true)
            {
                
                
                std::ofstream polyFile("poly.txt");
                //                    std::ofstream poly3DFile("poly3D.txt");
                
                polyFile<<insideReferencePoint.transpose()<<std::endl;
                for(const auto& node : polyPoints)
                {
                    polyFile<<"    "<<node.transpose()<<std::endl;
                    //                        poly3DFile<<node.transpose()<<std::endl;
                }
                
                std::ofstream pointsFile("points.txt");
                for(const auto& node : this->nodes())
                {
                    if(!node.second.expired())
                    {
                        pointsFile<<"    "<<node.second.lock()->sID<<" "<<node.first.transpose()<<std::endl;
                    }
                }
                
                std::ofstream edgesFile("edges.txt");
                for(const auto& patch : patches())
                {
                    for(const auto& edge : patch.second->edges())
                    {
                        edgesFile<<edge->tag()<<std::endl;
                    }
                }
                assert(this->isCompact() && "Plane not compact");
            }
        }
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
        
        std::cout<<"PeriodiPlane patches:"<<std::endl;
        for(const auto& patch : patches())
        {
            std::cout<<"Patch "<<patch.shift.transpose()<<std::endl;
            for(const auto& edge : patch.second.edges())
            {
                std::cout<<edge->source->sID<<"->"<<edge->sink->sID<<std::endl;
            }
        }
        
        
    }
    
    template<int dim>
    std::shared_ptr<PeriodicPlanePatch<dim>> PeriodicGlidePlane<dim>::getPatch(const VectorDim& shift)
    {
        return patches().get(shift);
    }
    
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
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO() :
    /* init */ glidePlaneID(0)
    /* init */,referencePlaneID(0)
    /* init */,shift(VectorDim::Zero())
    {
        
    }
    
    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch) :
    /* init */ glidePlaneID(patch.glidePlane->sID)
    /* init */,referencePlaneID(patch.periodicPlane->referencePlane->sID)
    /* init */,shift(patch.shift)
    {
        
    }
    
    template<int dim>
    PeriodicPlanePatchIO<dim>::PeriodicPlanePatchIO(std::stringstream& ss) :
    /* init */ glidePlaneID(0)
    /* init */,referencePlaneID(0)
    /* init */,shift(VectorDim::Zero())
    {
        ss>>glidePlaneID;
        ss>>referencePlaneID;
        for(int d=0;d<dim;++d)
        {
            ss>>shift(d);
        }
    }
    
    template <int dim,class T>
    T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
    {
        os  << ds.glidePlaneID<<"\t"
        /**/<< ds.referencePlaneID<<"\t"
        /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
        return os;
    }
    
}
#endif
