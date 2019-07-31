/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_H_
#define model_PeriodicGlidePlane_H_


#include<Eigen/StdVector>
#include <memory>
#include <string>


#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
#include <TerminalColors.h>




namespace model
{
    
    template<int dim>
    class PeriodicGlidePlane;
    
    template<int dim>
    class PeriodicPlanePatch;
    
    template<int dim>
    class PeriodicPlaneNode;
    
    template<int dim>
    class PeriodicPlaneEdge;
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlaneEdge
    {
        
        const PeriodicPlanePatch<dim>* const patch;
        const std::shared_ptr<PeriodicPlaneNode<dim>> source;
        const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
        const MeshBoundarySegment<dim>* const meshIntersection;
        PeriodicPlaneEdge<dim>* next;
        PeriodicPlaneEdge<dim>* prev;
        PeriodicPlaneEdge<dim>* twin;
        
        /**********************************************************************/
        PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                          const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                          const MeshBoundarySegment<dim>* const meshIntersection_in)
        ;
        ~PeriodicPlaneEdge();
        
        std::string tag() const;
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct NodalConnectivity
    {
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        
        PeriodicPlaneEdgeType* inEdge;
        PeriodicPlaneEdgeType* outEdge;
        
        NodalConnectivity() :
        /* init */ inEdge(nullptr)
        /* init */,outEdge(nullptr)
        {
            
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicPlaneNode : public StaticID<PeriodicPlaneNode<dim>>
    /*                   */,public Eigen::Matrix<double,dim-1,1>
    {
        
    public:
        
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef PeriodicPlaneNode<dim> PeriodicPlaneNodeType;
        typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
        typedef NodalConnectivity<dim> NodalConnectivityType;
        typedef std::map<size_t,NodalConnectivity<dim>> NodalConnectivityContainerType;
        typedef std::set<PeriodicPlaneEdgeType*> InOutEdgeContainerType;
        
    private:
        
        NodalConnectivityContainerType _patchConnectivities;
        NodalConnectivityContainerType _neighborConnectivities;
        
        
    public:
        
        /**********************************************************************/
        PeriodicPlaneNode(const VectorLowerDim& pos) :
        /* init*/ VectorLowerDim(pos)
        {
                        std::cout<<"Creating PeriodicPlaneNode "<<this->sID<<std::endl;
        }
        
        /**********************************************************************/
        ~PeriodicPlaneNode()
        {
                        std::cout<<"Destroying PeriodicPlaneNode "<<this->sID<<std::endl;
        }
        
        /**********************************************************************/
        void addLink(PeriodicPlaneEdgeType* const link)
        {
                        std::cout<<"PeriodicPlaneNode "<<this->sID<<" adding PeriodicPlaneEdge "<<link->tag()<<std::endl;
            
            
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
        
        /**********************************************************************/
        void removeLink(PeriodicPlaneEdgeType* const link)
        {
                        std::cout<<"PeriodicPlaneNode "<<this->sID<<" removing PeriodicPlaneEdge "<<link->tag()<<std::endl;
            
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
        
        /**********************************************************************/
        const NodalConnectivityContainerType& patchConnectivities() const
        {
            return _patchConnectivities;
        }
        
        /**********************************************************************/
        NodalConnectivityContainerType& patchConnectivities()
        {
            return _patchConnectivities;
        }
        
        /**********************************************************************/
        const NodalConnectivityContainerType& neighborConnectivities() const
        {
            return _neighborConnectivities;
        }
        
        /**********************************************************************/
        NodalConnectivityContainerType& neighborConnectivities()
        {
            return _neighborConnectivities;
        }
        
        
        
        /**********************************************************************/
        PeriodicPlaneEdgeType* inEdge() const
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
            return temp.size()==1? *temp.begin() : nullptr;
        }
        
        /**********************************************************************/
        PeriodicPlaneEdgeType* outEdge() const
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
            return temp.size()==1? *temp.begin() : nullptr;
        }
        
        
        
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
    /*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>,Eigen::aligned_allocator<std::shared_ptr<PeriodicPlaneEdge<dim>>>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        PeriodicGlidePlane<dim>* const periodicPlane;
        const VectorDim shift;
        const std::shared_ptr<GlidePlane<dim>> glidePlane;
        typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>,Eigen::aligned_allocator<std::shared_ptr<PeriodicPlaneEdge<dim>>>> PeriodicPlaneEdgeContainerType;
        
        /**********************************************************************/
        PeriodicPlanePatch(PeriodicGlidePlane<dim>& periodicPlane_in,
                           const VectorDim& shift_in) :
        /* init */ periodicPlane(&periodicPlane_in)
        /* init */,shift(shift_in)
        /* init */,glidePlane(periodicPlane->getGlidePlane(shift))
        {
            std::cout<<cyanColor<<"Creating patch "<<this->sID<<" with shift="<<std::setprecision(15)<<std::scientific<<shift.transpose()<<defaultColor<<std::endl;
            //            std::cout<<"meshIntersections.size()="<<glidePlane->meshIntersections.size()<<std::endl;
            
            if(BoundingMeshSegments<dim>::isRightHandedBoundary(glidePlane->meshIntersections,*periodicPlane->referencePlane))
            {
                addMeshIntersections(glidePlane->meshIntersections);
            }
            else
            {
                BoundingMeshSegments<dim> flippedTemp;
                for (auto it = glidePlane->meshIntersections.rbegin(); it != glidePlane->meshIntersections.rend(); ++it)
                {
                    flippedTemp.emplace_back(it->P1,it->P0,it->faces);
                }
                addMeshIntersections(flippedTemp);
            }
            
            periodicPlane->updateOuterBoundary();
            
        }
        
        ~PeriodicPlanePatch()
        {
                        std::cout<<"Destroying patch "<<this->sID<<std::endl;
            periodicPlane->updateOuterBoundary();
        }
        
        void addMeshIntersections(const BoundingMeshSegments<dim>& bms)
        {
            for(size_t k=0;k<bms.size();++k)
            {// compute 2d points on PeriodicPlane, and coneect them
                //size_t k1(k0!=bms.size()-1? k0+1 : 0);
                std::shared_ptr<PeriodicPlaneNode<dim>> source(periodicPlane->getSharedNode(bms[k].P0-shift));
                std::shared_ptr<PeriodicPlaneNode<dim>>   sink(periodicPlane->getSharedNode(bms[k].P1-shift));
                this->emplace_back(new PeriodicPlaneEdge<dim>(this,source,sink,&bms[k]));
            }
        }
        
        /**********************************************************************/
        const PeriodicPlaneEdgeContainerType& edges() const
        {
            return *this;
        }
        
    };
    
    template<int dim>
    struct TypeTraits<PeriodicGlidePlane<dim>>
    {
        typedef Eigen::Matrix<double,dim,1> KeyType;
        typedef PeriodicPlanePatch<dim> ValueType;
        typedef CompareVectorsByComponent<double,dim,float> CompareType;
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class PeriodicGlidePlane : private std::map<Eigen::Matrix<double,dim-1,1>,PeriodicPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>>
    /*                      */,private std::set<const PeriodicPlaneEdge<dim>*>
    /*                      */,private std::vector<const PeriodicPlaneEdge<dim>*>
    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>>
//    /*                      */,private std::map<Eigen::Matrix<double,dim,1>,PeriodicPlanePatch<dim>,CompareVectorsByComponent<double,dim,float>>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
//        typedef std::map<Eigen::Matrix<double,dim,1>,PeriodicPlanePatch<dim>,CompareVectorsByComponent<double,dim,float>> PatchContainerType;
        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>> PatchContainerType;
        typedef std::map<Eigen::Matrix<double,dim-1,1>,PeriodicPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        
        /**********************************************************************/
        static MatrixDim getL2G(const VectorDim& x,
                                const VectorDim& z)
        {
            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(xNorm>FLT_EPSILON);
            assert(zNorm>FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
            MatrixDim temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z/zNorm;
            temp.col(0)=x/xNorm;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp;
        }
        
    public:
        
        GlidePlaneFactory<dim>& glidePlaneFactory;
        const std::shared_ptr<GlidePlane<dim>> referencePlane;
        const MatrixDim L2G;
        
        PeriodicGlidePlane(GlidePlaneFactory<dim>& glidePlaneFactory_in,
                           const GlidePlaneKey<dim>& referencePlaneKey,
                           const VectorDim& globalX) :
        /* init */ glidePlaneFactory(glidePlaneFactory_in)
        /* init */,referencePlane(glidePlaneFactory.get(referencePlaneKey))
        /* init */,L2G(getL2G(globalX,referencePlane->unitNormal))
        {
            getPatch(VectorDim::Zero());
        }
        
        ~PeriodicGlidePlane()
        {
                        std::cout<<"DESTROYING PeriodicGlidePlane"<<std::endl;
        }
        
        bool isCompact() const
        {
            
            //            if(boundaryLinks().size()>=3)
            //            {
            //                const auto firstLink(*boundaryLinks().begin());
            //                auto currentLink(firstLink);
            //                size_t counter(1);
            //                while(currentLink->sink->sID!=firstLink->source->sID)
            //                {
            //                    currentLink=currentLink->sink->outEdge();
            //                    if(!currentLink)
            //                    {
            //                        return false;
            //                    }
            //                    counter++;
            //                }
            //                return boundaryLinks().size()==counter;
            //            }
            //            else
            //            {
            //                return false;
            //            }
            return outerBoundary().size()>=3;
        }
        
        void updateOuterBoundary()
        {
            std::cout<<"Updating OuterBoundary"<<std::endl;
            std::cout<<"old outerBoundarySize="<<outerBoundary().size()<<std::endl;
            outerBoundary().clear();
            outerBoundary().reserve(boundaryLinks().size());
            if(boundaryLinks().size()>=3)
            {
                outerBoundary().push_back(*boundaryLinks().begin());
                while(outerBoundary().back()->sink->sID!=outerBoundary().front()->source->sID)
                {
                    outerBoundary().push_back(outerBoundary().back()->sink->outEdge());
                    if(!outerBoundary().back())
                    {
                        outerBoundary().clear();
                        break;
                    }
                }
                if(boundaryLinks().size()!=outerBoundary().size())
                {
                    outerBoundary().clear();
                }
            }
                        std::cout<<"new outerBoundarySize="<<outerBoundary().size()<<std::endl;
        }
        
        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift)
        {
            return GlidePlaneKey<dim>(referencePlane->grain.grainID,referencePlane->P+shift,referencePlane->n);
        }
        
        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift)
        {
            return glidePlaneFactory.get(getGlidePlaneKey(shift));
        }
        
        /**********************************************************************/
        const PatchContainerType& patches() const
        {
            return *this;
        }
        
        PatchContainerType& patches()
        {
            return *this;
        }
        
        /**********************************************************************/
        NodeCointainerType& nodes()
        {
            return *this;
        }
        
        const NodeCointainerType& nodes() const
        {
            return *this;
        }
        
        std::vector<const PeriodicPlaneEdge<dim>*>& outerBoundary()
        {
            return *this;
        }
        
        const std::vector<const PeriodicPlaneEdge<dim>*>& outerBoundary() const
        {
            return *this;
        }
        
        const std::set<const PeriodicPlaneEdge<dim>*>& boundaryLinks() const
        {
            return *this;
        }
        
        std::set<const PeriodicPlaneEdge<dim>*>& boundaryLinks()
        {
            return *this;
        }
        
        void addBoundaryLink(const PeriodicPlaneEdge<dim>* link)
        {
            boundaryLinks().insert(link);
        }
        
        void removeBoundaryLink(const PeriodicPlaneEdge<dim>* link)
        {
            const size_t erased(boundaryLinks().erase(link));
            assert(erased==1 && "COULD NOT ERASE LINK FROM BOUNDARYLINKS");
        }
        
        
        
        
        //        int pnpoly(const Vector2dContainer& polygon, const Vector2d& test)
        //        {
        //            int i, j, c = 0;
        //            for (i = 0, j = polygon.size()-1; i < polygon.size(); j = i++)
        //            {
        //                if ( ((polygon[i](1)>test(1)) != (polygon[j](1)>test(1))) &&
        //                    (test(0) < (polygon[j](0)-polygon[i](0)) * (test(1)-polygon[i](1)) / (polygon[j](1)-polygon[i](1)) + polygon[i](0)) )
        //                {
        //                    c = !c;
        //                }
        //            }
        //            return c;
        //        }
        
        int isInsideOuterBoundary(const VectorLowerDim& test)
        {
            if(isCompact())
            {
                size_t i, j, c = 0;
                for (i = 0, j = outerBoundary().size()-1; i < outerBoundary().size(); j = i++)
                {
                    const VectorLowerDim& Pi(*outerBoundary()[i]->source);
                    const VectorLowerDim& Pj(*outerBoundary()[j]->source);

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
        
        /**********************************************************************/
        VectorDim getShift(const PeriodicPlaneEdge<dim>& edge) const
        {
            std::cout<<"getting shift for edge "<<edge.tag()<<std::endl;
            std::cout<<"edge meshIntersection->faces.size()="<<edge.meshIntersection->faces.size()<<std::endl;

            VectorDim shift(VectorDim::Zero());
            for(const PlanarMeshFace<dim>* const face : edge.meshIntersection->faces)
            {
                std::cout<<"faceID="<<face->sID<<std::endl;
                const auto parallelFaceID(edge.patch->glidePlane->grain.region.parallelFaces().at(face->sID));
                std::cout<<"paralleFaceID="<<parallelFaceID<<std::endl;
                std::cout<<"edge.patch->glidePlane->sID="<<edge.patch->glidePlane->sID<<std::endl;
                std::cout<<"edge.patch->glidePlane->grain.grainID="<<edge.patch->glidePlane->grain.grainID<<std::endl;
                std::cout<<"edge.patch->glidePlane->grain.region.faces().size()="<<edge.patch->glidePlane->grain.region.faces().size()<<std::endl;
//                std::cout<<"edge.patch->glidePlane->grain.region.faces().size()="<<edge.patch->glidePlane->grain.region.faces().size()<<std:endl;

                const auto parallelFace(edge.patch->glidePlane->grain.region.faces().at(parallelFaceID));
                std::cout<<"parallelFaceID="<<parallelFace->sID<<std::endl;
                std::cout<<"parallelFace->center()="<<parallelFace->center().transpose()<<std::endl;
                std::cout<<"parallelFace->outNormal()="<<parallelFace->outNormal().transpose()<<std::endl;
                std::cout<<"face->outNormal()="<<face->outNormal().transpose()<<std::endl;
                std::cout<<"face->center()="<<face->center().transpose()<<std::endl;

                shift+=(parallelFace->center()-face->center()).dot(face->outNormal())*face->outNormal();
                std::cout<<"partial_shift="<<shift.transpose()<<std::endl;

            }
            std::cout<<"shift="<<shift.transpose()<<std::endl;
            return shift;
        }
        

        /**********************************************************************/
        template<typename NodeType>
//        void addPatchesContainingPolygon(const std::vector<std::shared_ptr<NodeType>>& polyPoints)
        void addPatchesContainingPolygon(const std::vector<NodeType>& polyPoints)
        {
            
            //const double dt(1.0);
            
            // Compute a reference point inside the outerBoundary
            VectorLowerDim insideReferencePoint(VectorLowerDim::Zero());
            for(const auto& seg : referencePlane->meshIntersections)
            {
                insideReferencePoint+=getLocalPosition(seg.P0);
            }
            insideReferencePoint/=referencePlane->meshIntersections.size();
            assert(isInsideOuterBoundary(insideReferencePoint));
            
            for(const auto& polyPoint : polyPoints)
            {
                std::cout<<"Finding node "<<polyPoint.sID<<std::endl;
                const VectorLowerDim P(getLocalPosition(polyPoint.P));
//                const VectorLowerDim P(getLocalPosition(polyPoint.P+polyPoint.V*dt));
                while(!isInsideOuterBoundary(P))
                {
                    std::cout<<"point "<<polyPoint.sID<<" outside"<<std::endl;

                    const std::vector<const PeriodicPlaneEdge<dim>*> currentOuterBnd(outerBoundary());
                    for(const auto& bndEdge : currentOuterBnd)
                    {
                        std::cout<<"intersecting edge "<<bndEdge->tag()<<std::endl;

                        SegmentSegmentDistance<dim-1> ssd(P,insideReferencePoint,*bndEdge->source,*bndEdge->sink);
                        if(ssd.dMin<FLT_EPSILON)
                        {// intersection with current boundary found
                            std::cout<<"intersection found "<<std::endl;

                            getPatch(getShift(*bndEdge));
                            break;
                        }
                    }
                }
                std::cout<<"point "<<polyPoint.sID<<" inside"<<std::endl;

                
            }


            
        }
        
        
        /**********************************************************************/
        void print()
        {
            std::cout<<"PeriodiPlane nodes:"<<std::endl;
            for(const auto& node : nodes())
            {
                std::cout<<node.second->sID<<" "<<node.first.transpose()<<std::endl;
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
        
        /**********************************************************************/
        std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift)
        {
            std::cout<<"Getting patch "<<shift.transpose()<<std::endl;
            std::cout<<"patch key exists "<<(patches().find(shift)!=patches().end())<<std::endl;
            return patches().get(shift);
//            const auto
//            const auto newPatch(patches().emplace(std::piecewise_construct,
//                                                  std::forward_as_tuple(shift),
//                                                  std::forward_as_tuple(this,
//                                                                        shift)).first->second);
        }
        
        /**********************************************************************/
        VectorLowerDim getLocalPosition(const VectorDim& point) const
        {
            const VectorDim pointLocal(L2G.transpose()*(point-referencePlane->P));
            assert(fabs(pointLocal(2))<FLT_EPSILON);
            return pointLocal.template segment<2>(0);
        }
        
        
        /**********************************************************************/
        std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim)
        {
            const VectorLowerDim point(getLocalPosition(pointDim));
            const auto iter(nodes().find(point));
            if(iter==nodes().end())
            {// point does not exist
                std::shared_ptr<PeriodicPlaneNode<dim>> newNode(new PeriodicPlaneNode<dim>(point));
                nodes().emplace(point,newNode.get());
                return newNode;
            }
            else
            {
                assert(iter->second->patchConnectivities().size()>0 && "EXISTING NODE IS ISOLATED");
                const auto patchConnectivity(iter->second->patchConnectivities().begin()->second);
                if(patchConnectivity.outEdge)
                {
                    return patchConnectivity.outEdge->source;
                }
                else
                {
                    if(patchConnectivity.inEdge)
                    {
                        return patchConnectivity.inEdge->sink;
                    }
                    else
                    {
                        assert(false && "EXISTING POSITION IS NOT CONNECTED");
                        return std::shared_ptr<PeriodicPlaneNode<dim> >(nullptr);
                    }
                }
            }
        }
        
    };
    
    
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlaneEdge<dim>::PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                                              const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                                              const MeshBoundarySegment<dim>* const meshIntersection_in) :
    /* init */ patch(patch_in)
    /* init */,source(source_in)
    /* init */,sink(sink_in)
    /* init */,meshIntersection(meshIntersection_in)
    /* init */,next(nullptr)
    /* init */,prev(nullptr)
    /* init */,twin(nullptr)
    {
                std::cout<<"Creating PeriodicPlaneEdge "<<this->tag()<<std::endl;
        
        
        source->addLink(this);
        sink->addLink(this);
        if(twin)
        {
            patch->periodicPlane->removeBoundaryLink(twin);
        }
        else
        {
            patch->periodicPlane->addBoundaryLink(this);
        }
    }
    
    /**********************************************************************/
    template<int dim>
    PeriodicPlaneEdge<dim>::~PeriodicPlaneEdge()
    {
                std::cout<<"Destroying PeriodicPlaneEdge "<<this->tag()<<std::endl;
        
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
            patch->periodicPlane->addBoundaryLink(twin);
        }
        else
        {
            patch->periodicPlane->removeBoundaryLink(this);
        }
    }
    
    /**********************************************************************/
    template<int dim>
    std::string PeriodicPlaneEdge<dim>::tag() const
    {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
      */
        return std::to_string(source->sID) + "->" + std::to_string(sink->sID);
    }
}
#endif
