/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNODE_H_
#define model_DISLOCATIONNODE_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <tuple>

#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/Geometry/Splines/SplineNode.h>
#include <model/Math/GramSchmidt.h>
#include <model/Mesh/Simplex.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/PlanePlaneIntersection.h>
#include <model/Geometry/PlaneLineIntersection.h>
#include <model/Geometry/LineSegment.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/DefectiveCrystalParameters.h>


#ifndef NDEBUG
#define VerboseDislocationNode(N,x) if(verboseDislocationNode>=N){model::cout<<x;}
#else
#define VerboseDislocationNode(N,x)
#endif

namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationNode :
    /*          */ public SplineNode<DislocationNode<_dim,corder,InterpolationType>,_dim,corder,InterpolationType>,
    /*          */ private std::set<const GrainBoundary<_dim>*>,
    /*          */ private std::set<const Grain<_dim>*>,
    /*          */ private std::set<const MeshPlane<_dim>*>,
    /*          */ private BoundingLineSegments<_dim>
    {
        
    public:
        
        constexpr static int dim=_dim; // make dim available outside class
        typedef DislocationNode       <dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment    <dim,corder,InterpolationType> LinkType;
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        typedef typename NodeBaseType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef typename TypeTraits<NodeType>::LoopNetworkType LoopNetworkType;
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef typename TypeTraits<NodeType>::MeshLocation MeshLocation;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef std::set<const MeshPlaneType*> MeshPlaneContainerType;
        typedef std::set<const Grain<dim>*> GrainContainerType;
        typedef std::set<const GrainBoundary<dim>*> GrainBoundaryContainerType;
        typedef typename NodeBaseType::NeighborContainerType NeighborContainerType;
        typedef typename NodeBaseType::LoopLinkContainerType LoopLinkContainerType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        static const double bndTol;
        static int verboseDislocationNode;
        
    private:
        
        /**********************************************************************/
        void updateMeshPlaneIntersections(const MeshPlaneType& lastGlidePlane)
        {
            BoundingLineSegments<dim> temp;
            
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections"<<std::endl;);
            VerboseDislocationNode(2,"  lastGlidePlane.P="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.P.transpose()<<std::endl;);
            VerboseDislocationNode(2,"  lastGlidePlane.unitNormal="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.unitNormal.transpose()<<std::endl;);
            
            switch (meshPlanes().size())
            {
                case 0:
                {// there must be at least one glide plane
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                    VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
                    _glidePlaneIntersections.clear();
                    break;
                }
                    
                case 2:
                {// a second plane is being added, so we must have no _glidePlaneIntersections
                    VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                    
                    // Grab the infinite line of intersection between the two planes
                    const PlanePlaneIntersection<dim>& ppi(this->network().glidePlaneIntersection(&meshPlane(0),&meshPlane(1)));
                    
                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                    {/* Two distinct glide planes can be coincident only if they belong to different grains
                      * In that case, the intersection of their bounding boxes should be one line segment
                      */
                        VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2a"<<std::endl;);
                        if(boundingBoxSegments().size()!=1)
                        {
                            model::cout<<"DislocationNode "<<this->sID<<std::endl;
                            model::cout<<"glidePlane(0) is "<<meshPlane(0).P.transpose()<<","<<meshPlane(0).unitNormal.transpose()<<std::endl;
                            model::cout<<"glidePlane(1) is "<<meshPlane(1).P.transpose()<<","<<meshPlane(1).unitNormal.transpose()<<std::endl;
                            assert(false && "There should be only one line in boundingBoxSegments()");
                        }
                        //assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
                        _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
                    }
                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                    {/* If the two planes are incident then the intersection of
                      * their bounding boxes is either a pair of singluar segments (2 points)
                      * or a line segment on the boundary
                      */
                        switch (boundingBoxSegments().size())
                        {
                            case 1:
                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
                                VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2b"<<std::endl;);
                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
                                VerboseDislocationNode(4,"  boundingBoxSegments()[0].first="<<boundingBoxSegments()[0].first.transpose()<<std::endl;);
                                VerboseDislocationNode(4,"  boundingBoxSegments()[0].second="<<boundingBoxSegments()[0].second.transpose()<<std::endl;);
                                
                                break;
                            }
                                
                            case 2:
                            {// The two intersections must be degenerate (2 boundary points)
                                //                                std::cout<<boundingBoxSegments()[0].first.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments()[0].second.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments()[1].first.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments()[1].second.transpose()<<std::endl;
                                VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2c"<<std::endl;);
                                assert((boundingBoxSegments()[0].first-boundingBoxSegments()[0].second).squaredNorm()<FLT_EPSILON);
                                assert((boundingBoxSegments()[1].first-boundingBoxSegments()[1].second).squaredNorm()<FLT_EPSILON);
                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[1].first);
                                break;
                            }
                                
                            default:
                            {
                                model::cout<<"DislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
                                for(const auto& pair : boundingBoxSegments())
                                {
                                    model::cout<<"("<<pair.first.transpose()<<","<<pair.second.transpose()<<")"<<std::endl;
                                }
                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
                            }
                        }
                    }
                    else
                    {
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                    // Now we must have exactly one _glidePlaneIntersections
                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must have size 1");
                    
                    break;
                }
                    
                default:
                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must exist");
                    
                    // intersect the _glidePlaneIntersections with the new plane
                    PlaneLineIntersection<dim> pli(lastGlidePlane.P,
                                                   lastGlidePlane.unitNormal,
                                                   _glidePlaneIntersections[0].first, // origin of line
                                                   _glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first // line direction
                                                   );
                    
                    if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
                    {// nothing to do, _glidePlaneIntersections remains unchanged
                        
                    }
                    else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
                    {// _glidePlaneIntersections becomes a singular point
                        _glidePlaneIntersections[0].first =pli.P;
                        _glidePlaneIntersections[0].second=pli.P;
                    }
                    else
                    {
                        model::cout<<"DislocationNode "<<this->sID<<std::endl;
                        model::cout<<"MeshPlanes are:"<<std::endl;
                        for(const auto& plane : meshPlanes())
                        {
                            model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
                        }
                        model::cout<<"MeshPlane intersection is:"<<std::endl;
                        model::cout<<std::setprecision(15)<<std::scientific<<"  P1="<<_glidePlaneIntersections[0].first.transpose()<<", P2="<<_glidePlaneIntersections[0].second.transpose()<<std::endl;
                        
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                }
                    
            }
            
            
            VerboseDislocationNode(2,"  _glidePlaneIntersections are"<<std::endl;);
            for(const auto& pair : _glidePlaneIntersections)
            {
                VerboseDislocationNode(2,"P1="<<std::setprecision(15)<<std::scientific<<pair.first.transpose()<<", P2="<<pair.second.transpose()<<std::endl;);
                
            }
            
            assert(_glidePlaneIntersections.size()<=1 && "_glidePlaneIntersections can have at the most size 1");
        }
        
        /**********************************************************************/
        bool addMeshPlane(const MeshPlaneType& gp)
        {
            const bool success=meshPlanes().insert(&gp).second;
            if(success)
            {
                VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" addGlidePlane. meshPlanes().size()="<<meshPlanes().size()<<std::endl;);
                VerboseDislocationNode(4,"Current bounding box"<<std::endl;);
                VerboseDislocationNode(4,boundingBoxSegments()<<std::endl;);
                VerboseDislocationNode(4,"Plane bounding box"<<std::endl;);
                VerboseDislocationNode(4,BoundingLineSegments<dim>(gp)<<std::endl;);
                VerboseDislocationNode(4,"Plane meshIntersections: "<<gp.meshIntersections.size()<<std::endl;);
                //VerboseDislocationNode(4,BoundingLineSegments<dim>(gp)<<std::endl;);
                
                
                
                assert(gp.contains(this->get_P()) && "Glide Plane does not contain DislocationNode");
                boundingBoxSegments().updateWithMeshPlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                assert((boundingBoxSegments().size() || !_isOnBoundingBox) && "EMPTY boundingBoxSegments");
                updateMeshPlaneIntersections(gp);
                
                VerboseDislocationNode(4,"new bounding box"<<std::endl;);
                VerboseDislocationNode(4,boundingBoxSegments()<<std::endl;);
                
                
                //                grains().insert(&this->network().poly.grain(gp.regionIDs.first));    // Insert new grain in grainSet
                //                grains().insert(&this->network().poly.grain(gp.regionIDs.second));   // Insert new grain in grainSet
            }
            return success;
        }
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes()
        {
            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" adding GrainBoundaryPlanes"<<std::endl;);
            
            size_t addedGp=0;
            // Check if node is on a GB
            for(const auto& grain : grains())
            {
                for(const auto& gb : grain->grainBoundaries())
                {
                    VerboseDislocationNode(4,"GB "<<gb.second->tag()<<", d="<<gb.second->distanceTo(this->get_P())<<", contained="<<gb.second->contains(this->get_P())<<std::endl;);
                    if(gb.second->contains(this->get_P()))
                    {
                        grainBoundaries().insert(gb.second);
                        addedGp+=addMeshPlane(*gb.second);
                    }
                }
            }
            //            VerboseDislocationNode(3,"added "<<addedGp<<" planes"<<std::endl;);
            
            
            //            for(const auto& gb : this->network().poly.grainBoundaries())
            //            {
            //                if(gb.second.contains(this->get_P()))
            //                {
            //                    grainBoundaries().insert(&gb.second);
            //                    addedGp+=addMeshPlane(gb.second);
            //                }
            //            }
            
            if(isGrainBoundaryNode())
            {
                VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" added "<<addedGp<<" GrainBoundaryPlanes"<<std::endl;);
                for(const auto& pair : this->neighbors())
                {
                    std::get<1>(pair.second)->addGrainBoundaryPlanes();
                }
            }
            
            return addedGp;
        }
        
        
        
        /**********************************************************************/
        VectorDim snapToBoundingBox(const VectorDim& P) const
        {/*!\param[in] P position to be snapped to the bounding box
          * \returns a point on the bounding box close to P. The returned point
          * is the closest to the bounding box, unless the closest point causes
          * boundarySegments to become interior. In that case the closest boundary
          * vertex is returned.
          */
            
            VerboseDislocationNode(4,"snapping P="<<P.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            std::multimap<double,VectorDim> snapMap;
            
            // Collect possible snap points, sorted by distance to P
            for(size_t k=0;k<boundingBoxSegments().size();++k)
            {
                const std::pair<VectorDim,VectorDim>& vertexPair(boundingBoxSegments()[k]);
                const VectorDim segm(vertexPair.second-vertexPair.first);
                const double segmNorm2(segm.squaredNorm());
                if(segmNorm2>FLT_EPSILON)
                {
                    snapMap.emplace((P-vertexPair.first).squaredNorm(),vertexPair.first);
                    snapMap.emplace((P-vertexPair.second).squaredNorm(),vertexPair.second);
                    
                    double u((P-vertexPair.first).dot(segm)/segmNorm2);
                    if(u>0.0 && u<1.0)
                    {
                        const VectorDim x(vertexPair.first+u*segm);
                        snapMap.emplace((P-x).squaredNorm(),x);
                    }
                }
                else
                {
                    const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
                    snapMap.emplace((P-x).squaredNorm(),x);
                }
            }
            
            VerboseDislocationNode(4,"there are "<<snapMap.size()<<" possible snap points."<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            // Return the first point to which we can snap
            for(const auto& pair : snapMap)
            {
                if(isMovableTo(pair.second))
                {
                    return pair.second;
                }
            }
            
            assert(false && "snapToBoundingBox FAILED.");
            
            
            //
            //
            //            const VectorDim pL=std::get<0>(boundingBoxSegments().snap(P));
            //            const VectorDim pV=boundingBoxSegments().snapToVertex(P).second;
            //
            //
            //            bool pLcontained=true;
            //            bool pVcontained=true;
            //
            //            std::deque<VectorDim> bndChords;
            //
            //            for(const auto& pair : this->neighbors())
            //            {
            //                VerboseDislocationNode(3,"checking "<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<std::endl;);
            //
            //                if(std::get<1>(pair.second)->isBoundarySegment())
            //                {// boundary segments must not become internal
            //                    VerboseDislocationNode(3,std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<" is boundary"<<std::endl;);
            //
            //                    pLcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pL+std::get<0>(pair.second)->get_P())).first;
            //                    pVcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pV+std::get<0>(pair.second)->get_P())).first;
            //                    VerboseDislocationNode(3,"pLcontained="<<pLcontained<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            //                    VerboseDislocationNode(3,"pVcontained="<<pVcontained<<std::endl;);//", lineID="<<pVcontained.second<<std::endl;
            //
            //                    bndChords.push_back(std::get<1>(pair.second)->chord().normalized());
            //                }
            //
            //                if(std::get<1>(pair.second)->isGrainBoundarySegment())
            //                {// grainBoundary segments must not become internal
            //                    for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
            //                    {
            //                        pLcontained*=gb->contains(pL);
            //                        pVcontained*=gb->contains(pV);
            //
            //                        bndChords.push_back(std::get<1>(pair.second)->chord().normalized());
            //                    }
            //                }
            //            }
            //
            //            bool parallelBndChords=true;
            //            for(size_t i=0;i<bndChords.size();++i)
            //            {
            //                for(size_t j=i+1;j<bndChords.size();++j)
            //                {
            //                    parallelBndChords*=(bndChords[i].cross(bndChords[j]).norm()<FLT_EPSILON);
            //                }
            //            }
            //
            //            if(pLcontained && parallelBndChords)
            //            {
            //                VerboseDislocationNode(4,"snapping to pL="<<pL.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            //                return pL;
            //            }
            //            else
            //            {
            //                if(pVcontained)
            //                {
            //                    VerboseDislocationNode(4,"snapping to pV="<<pV.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            //                    return pV;
            //                }
            //                else
            //                {
            //                    model::cout<<"DislocationNode "<<this->sID<<" snapToBoundingBox FAILED."<<std::endl;
            //                    assert(false && "snapToBoundingBox FAILED.");
            //                    return VectorDim::Zero();
            //                }
            //            }
            
        }
        
        //        /**********************************************************************/
        //        VectorDim snapToBoundingBox(const VectorDim& P) const
        //        {/*!\param[in] P position to be snapped to the bounding box
        //          * \returns a point on the bounding box close to P. The returned point
        //          * is the closest to the bounding box, unless the closest point causes
        //          * boundarySegments to become interior. In that case the closest boundary
        //          * vertex is returned.
        //          */
        //
        //            VerboseDislocationNode(4,"snapping P="<<P.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
        //
        //
        //            const VectorDim pL=std::get<0>(boundingBoxSegments().snap(P));
        //            const VectorDim pV=boundingBoxSegments().snapToVertex(P).second;
        //
        //
        //            bool pLcontained=true;
        //            bool pVcontained=true;
        //
        //            std::deque<VectorDim> bndChords;
        //
        //            for(const auto& pair : this->neighbors())
        //            {
        //                VerboseDislocationNode(3,"checking "<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<std::endl;);
        //
        //                if(std::get<1>(pair.second)->isBoundarySegment())
        //                {// boundary segments must not become internal
        //                    VerboseDislocationNode(3,std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<" is boundary"<<std::endl;);
        //
        //                    pLcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pL+std::get<0>(pair.second)->get_P())).first;
        //                    pVcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pV+std::get<0>(pair.second)->get_P())).first;
        //                    VerboseDislocationNode(3,"pLcontained="<<pLcontained<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
        //                    VerboseDislocationNode(3,"pVcontained="<<pVcontained<<std::endl;);//", lineID="<<pVcontained.second<<std::endl;
        //
        //                    bndChords.push_back(std::get<1>(pair.second)->chord().normalized());
        //                }
        //
        //                if(std::get<1>(pair.second)->isGrainBoundarySegment())
        //                {// grainBoundary segments must not become internal
        //                    for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
        //                    {
        //                        pLcontained*=gb->contains(pL);
        //                        pVcontained*=gb->contains(pV);
        //
        //                        bndChords.push_back(std::get<1>(pair.second)->chord().normalized());
        //                    }
        //                }
        //            }
        //
        //            bool parallelBndChords=true;
        //            for(size_t i=0;i<bndChords.size();++i)
        //            {
        //                for(size_t j=i+1;j<bndChords.size();++j)
        //                {
        //                    parallelBndChords*=(bndChords[i].cross(bndChords[j]).norm()<FLT_EPSILON);
        //                }
        //            }
        //
        //            if(pLcontained && parallelBndChords)
        //            {
        //                VerboseDislocationNode(4,"snapping to pL="<<pL.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
        //                return pL;
        //            }
        //            else
        //            {
        //                if(pVcontained)
        //                {
        //                    VerboseDislocationNode(4,"snapping to pV="<<pV.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
        //                    return pV;
        //                }
        //                else
        //                {
        //                    model::cout<<"DislocationNode "<<this->sID<<" snapToBoundingBox FAILED."<<std::endl;
        //                    assert(false && "snapToBoundingBox FAILED.");
        //                    return VectorDim::Zero();
        //                }
        //            }
        //
        //        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
//            if (this->network().use_boundary)
//            {
                if (guess==NULL)
                {
                    temp=this->network().mesh.search(this->get_P());
                }
                else
                {
                    if(grains().size()==1)
                    {// node only in one region
                        if((*grains().begin())->grainID!=guess->region->regionID)
                        {
                            temp=this->network().mesh.searchRegion((*grains().begin())->grainID,this->get_P());
                        }
                        else
                        {
                            temp=this->network().mesh.searchRegionWithGuess(this->get_P(),guess);
                        }
                    }
                    else
                    {
                        temp=this->network().mesh.searchWithGuess(this->get_P(),guess);
                    }
                }
                
                if(!temp.first) // DislocationNode not found inside mesh
                {
                    // Detect if the DislocationNode is sligtly outside the boundary
                    int faceID;
                    const double baryMin(temp.second->pos2bary(this->get_P()).minCoeff(&faceID));
                    const bool isApproxOnBoundary(std::fabs(baryMin)<1.0e3*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
                    if(!isApproxOnBoundary)
                    {
                        model::cout<<"DislocationNode "<<this->sID<<" @ "<<this->get_P().transpose()<<std::endl;
                        model::cout<<"Simplex "<<temp.second->xID<<std::endl;
                        model::cout<<"bary "<<temp.second->pos2bary(this->get_P())<<std::endl;
                        model::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                        model::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                        model::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                        assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
                    }
                }
//            }
            
            return temp.second;
        }
        
        
        //        /**********************************************************************/
        //        void make_projectionMatrix()
        //        {
        //
        //            Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
        //            VectorOfNormalsType  CN;
        //            for(const auto& plane : meshPlanes())
        //            {
        //                CN.push_back(plane->unitNormal);
        //            }
        //
        //            if(isBoundaryNode())
        //            {
        //                CN.push_back(boundaryNormal);
        //
        //            }
        //
        //            // Add normal to region boundary
        //            //            CN.push_back(regionBndNormal);
        //
        //            // Find independent vectors
        //            GramSchmidt::orthoNormalize(CN);
        //
        //            // Assemble projection matrix (prjM)
        //            this->prjM.setIdentity();
        //            for (size_t k=0;k<CN.size();++k)
        //            {
        //                this->prjM*=( I-CN[k]*CN[k].transpose() );
        //            }
        //
        //        }
        
        /**********************************************************************/
        BoundingLineSegments<dim> _glidePlaneIntersections; //
        
        
        bool _isGlissile;
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        //! The previous velocity vector of *this DislocationNode
        VectorDofType vOld;
        double velocityReductionCoeff;
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        bool _isOnBoundingBox;
        //        bool _isGrainBoundaryNode;
        VectorDim boundaryNormal;
        //        VectorDim C;
        
        std::shared_ptr<NodeType> virtualNode;
        
        
    public:
        
        const bool isVirtualBoundaryNode;

        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            use_velocityFilter=TextFileParser(fileName).readScalar<double>("use_velocityFilter",true);
            velocityReductionFactor=TextFileParser(fileName).readScalar<double>("velocityReductionFactor",true);
            assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
            //            EDR.readScalarInFile(fullName.str(),"verboseDislocationNode",NodeType::verboseDislocationNode);
            verboseDislocationNode=TextFileParser(fileName).readScalar<int>("verboseDislocationNode",true);
        }
        
        
        /**********************************************************************/
        DislocationNode(LoopNetworkType* const ln,
                        const VectorDim& Pin,
                        const VectorDofType& Vin,
                        const double& vrc) :
        /* base constructor */ NodeBaseType(ln,Pin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex((const Simplex<dim,dim>*) NULL)),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(vrc),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol))
        /* init list        */,isVirtualBoundaryNode(false)
        {/*! Constructor from DOF
          */
            VerboseDislocationNode(1,"Creating DislocationNode "<<this->sID<<" from position"<<std::endl;);
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const VectorDim& Pin) :
        /* base constructor */ NodeBaseType(pL.loopNetwork,Pin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        //        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ velocityReductionCoeff(std::min(pL.source->velocityReduction(),pL.sink->velocityReduction())),
        /* init list        */ _isOnBoundingBox(pL.isBoundarySegment()),
        /* init list        */ boundaryNormal(SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) )
        /* init list        */,isVirtualBoundaryNode(false)
        {/*! Constructor from ExpandingEdge and DOF
          */
            VerboseDislocationNode(1,"Creating DislocationNode "<<this->sID<<" from expanding "<<pL.source->sID<<"->"<<pL.sink->sID<<std::endl;);
        }
        
        /**********************************************************************/
        DislocationNode(LoopNetworkType* const ln,
                        const VectorDim& Pin) :
        /* base constructor */ NodeBaseType(ln,Pin),
        /* init list        */ _isGlissile(false),
        /* init list        */ p_Simplex(this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC? get_includingSimplex((const Simplex<dim,dim>*) NULL) : NULL),
        /* init list        */ velocity(VectorDim::Zero()),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(VectorDim::Zero())
        /* init list        */,isVirtualBoundaryNode(true)
        {/*! Constructor from DOF
          */
            VerboseDislocationNode(1,"Creating VirtualDislocationNode "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~DislocationNode()
        {
            VerboseDislocationNode(1,"Destroying DislocationNode "<<this->sID<<" ("<<this<<")"<<std::endl;);
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<", virtual node count="<<virtualNode.use_count()<<std::endl;);
            
            if(virtualNode)
            {
//                VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" removing virtual "<<virtualNode->sID<<std::endl;);
                this->network().remove(virtualNode->sID);
            }
//
        }
        
        /**********************************************************************/
        const std::shared_ptr<NodeType>& virtualBoundaryNode() const
        {
            return virtualNode;
            
//            std::shared_ptr<NodeType> temp(nullptr);
//            
//            if(this->network().useVirtualExternalLoops /* && this->network().use_bvp */ && isBoundaryNode()) //ENABLE THIS USE_BVP
//            {
//                std::vector<LinkType*> virtualNeighbors;
//                for (const auto& neighborIter : this->neighbors())
//                {
//                    if(   std::get<1>(neighborIter.second)->isVirtualBoundarySegment()
//                       && !std::get<1>(neighborIter.second)->hasZeroBurgers())
//                    {
//                        virtualNeighbors.push_back(std::get<1>(neighborIter.second));
//                    }
//                }
//                
//                switch (virtualNeighbors.size())
//                {// a unique virtual node does not exist, create a new one
//                    case 0:
//                    {
//                        temp.reset(new NodeType(&this->network(),this->get_P()+100.0*boundaryNormal));
//                        break;
//                    }
//                        
//                    case 1:
//                    {// a unique virtual node already exists, return that
//                        temp=virtualNeighbors[0]->source->sID==this->sID? virtualNeighbors[0]->sink : virtualNeighbors[0]->source;
//                        assert(temp->isVirtualBoundaryNode);
//                        break;
//                    }
//                        
//                    default:
//                    {
//                        assert(false && "THERE CAN BE AT THE MOST ONE VIRTUAL NODE FOR EACH BOUNDARY NODE");
//                        break;
//                    }
//                }
//            
//            }
//            
//            return temp;
        }
        
        /**********************************************************************/
        const MeshPlaneContainerType& meshPlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        MeshPlaneContainerType& meshPlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const MeshPlaneType& meshPlane(const size_t& n) const
        {
            assert(n<meshPlanes().size());
            auto iter=meshPlanes().begin();
            std::advance(iter,n);
            return **iter;
        }
        
        /**********************************************************************/
        const GrainContainerType& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainContainerType& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        GrainBoundaryContainerType& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundaryContainerType& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        /**********************************************************************/
        BoundingLineSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& glidePlaneIntersections() const
        {
            
            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        VectorDim snapToMeshPlaneIntersection(const VectorDim& P)
        {
            
            switch (_glidePlaneIntersections.size())
            {
                case 0:
                {
                    //                    assert(glidePlanes().size()>0);
                    VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 0"<<std::endl;);
                    assert(meshPlanes().size()==1);
                    return meshPlane(0).snapToPlane(P);
                    break;
                }
                    
                case 1:
                {
                    const VectorDim D=_glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first;
                    const double normD2(D.squaredNorm());
                    if(normD2>FLT_EPSILON)
                    {
                        const double u=(P-_glidePlaneIntersections[0].first).dot(D)/normD2;
                        VerboseDislocationNode(3,"u="<<u<<std::endl;);
                        
                        if(u<0.0)
                        {
                            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1a"<<std::endl;);
                            return _glidePlaneIntersections[0].first;
                        }
                        else if(u>1.0)
                        {
                            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1b"<<std::endl;);
                            return _glidePlaneIntersections[0].second;
                        }
                        else
                        {
                            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1c"<<std::endl;);
                            return _glidePlaneIntersections[0].first+u*D;
                        }
                    }
                    else
                    {
                        VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1d"<<std::endl;);
                        return _glidePlaneIntersections[0].first;
                    }
                    
                    break;
                }
                    
                default:
                {
                    VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 2"<<std::endl;);
                    assert(0 && "THERE CAN BE AT MOST ONE LINE OF INTERSECTION");
                    return VectorDim::Zero();
                    break;
                }
            }
            
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" addLoopLink "<<pL->tag()<<std::endl;);
            
            NodeBaseType::addLoopLink(pL); // forward to base class
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            if(!pL->loop()->isVirtualBoundaryLoop())
            {
                const bool success = addMeshPlane(*pL->loop()->glidePlane.get());
                if(success)
                {
                    grains().insert(&this->network().poly.grain(pL->loop()->grain.grainID));    // Insert new grain in grainSet
                    _isGlissile*=pL->loop()->isGlissile;
                }
                
                addGrainBoundaryPlanes();
                pL->pLink->addGrainBoundaryPlanes();
                
                if(boundingBoxSegments().contains(this->get_P()).first)
                {
                    _isOnBoundingBox=true;
                    setToBoundary(this->get_P());
                }
                
            }
            
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" finished addLoopLink "<<pL->tag()<<std::endl;);

            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
            
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            if(!pL->loop()->isVirtualBoundaryLoop())
            {// Re-construct nodeConfinement
                
                _isGlissile=true;
                meshPlanes().clear();
                boundingBoxSegments().clear();
                _glidePlaneIntersections.clear();
                grains().clear();
                
                for(const auto& loopLink : this->loopLinks())
                {
                    if(loopLink->loop()->glidePlane)
                    {
                        const bool success = addMeshPlane(*loopLink->loop()->glidePlane.get());
                        if(success)
                        {
                            grains().insert(&this->network().poly.grain(loopLink->loop()->grain.grainID));    // Insert new grain in grainSet
                            _isGlissile*=loopLink->loop()->isGlissile;
                        }
                    }
                }
                
                addGrainBoundaryPlanes();
                pL->pLink->addGrainBoundaryPlanes();
                
                if(boundingBoxSegments().contains(this->get_P()).first)
                {
                    _isOnBoundingBox=true;
                    setToBoundary(this->get_P());
                }
            }
            
            VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" finished removeLoopLink "<<pL->tag()<<std::endl;);

        }
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const __attribute__ ((deprecated)) // REMOVE THIS FUNCTION AFTER CHANGING REMESH AND DISLOCATION NETWORK COMPONENT
        {
            VectorOfNormalsType temp;
            
            if(_isGlissile)
            {
                for(const auto& plane : meshPlanes())
                {
                    temp.push_back(plane->unitNormal);
                }
                temp.push_back(boundaryNormal);
                GramSchmidt::orthoNormalize(temp);
                assert(temp.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
            }
            else
            {
                temp.push_back((VectorDim()<<1.0,0.0,0.0).finished());
                temp.push_back((VectorDim()<<0.0,1.0,0.0).finished());
                temp.push_back((VectorDim()<<0.0,0.0,1.0).finished());
            }
            
            return temp;
        }
        
        
        //        /**********************************************************************/
        //        void set_V(const VectorDofType& vNew)
        //        {
        //            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
        //        }
        
        /**********************************************************************/
        void projectVelocity()
        {
            //            if(velocity.squaredNorm()>FLT_EPSILON)
            //            {
            VectorOfNormalsType temp;
            
            for(const auto& loop : this->loops())
            {
                if(loop->isGlissile)
                {
                    if(loop->glidePlane)
                    {
                        temp.push_back(loop->glidePlane->unitNormal);
                    }
                }
                else
                {
                    velocity.setZero();
                    break;
                }
            }
            temp.push_back(boundaryNormal);
            
            if(velocity.squaredNorm()>FLT_EPSILON)
            {
                GramSchmidt::orthoNormalize(temp);
                
                for(const auto& vec : temp)
                {
                    velocity-=velocity.dot(vec)*vec;
                }
                
            }
            //            }
        }
        
        
        /**********************************************************************/
        void set_V(const VectorDofType& vNew)
        {
            vOld=velocity; // store current value of velocity before updating
            
            //            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
            velocity=vNew;
            projectVelocity();
            
            if(use_velocityFilter)
            {
                const double filterThreshold=0.05*velocity.norm()*vOld.norm();
                
                if(velocity.dot(vOld)<-filterThreshold)
                {
                    velocityReductionCoeff*=velocityReductionFactor;
                }
                else if(velocity.dot(vOld)>filterThreshold)
                {
                    velocityReductionCoeff/=velocityReductionFactor;
                }
                else
                {
                    // don't change velocityReductionCoeff
                }
                if(velocityReductionCoeff>1.0)
                {
                    velocityReductionCoeff=1.0;
                }
                if(velocityReductionCoeff<0.005)
                {
                    velocityReductionCoeff=0.005;
                }
                velocity*=velocityReductionCoeff;
            }
        }
        
        /**********************************************************************/
        const VectorDofType& get_V() const
        {/*! The nodal velocity vector
          */
            return velocity;
        }
        
        //        /**********************************************************************/
        //        bool isOscillating() const
        //        {
        //            return velocityReductionCoeff<std::pow(velocityReductionFactor,3);
        //        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* includingSimplex() const
        {/*!\returns A pointer to the const Simplex imcluding *this DislocationNode
          */
            return p_Simplex;
        }
        
        /**********************************************************************/
        const VectorDim& bndNormal() const
        {
            return boundaryNormal;
        }
        
        /**********************************************************************/
        MeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            MeshLocation temp = MeshLocation::outsideMesh;
            
            if(!isVirtualBoundaryNode)
            {
                if(isBoundaryNode())
                {
                    temp=MeshLocation::onMeshBoundary;
                }
                else
                {
                    if(isGrainBoundaryNode())
                    {
                        temp=MeshLocation::onRegionBoundary;
                    }
                    else
                    {
                        temp=MeshLocation::insideMesh;
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        const bool& isOnBoundingBox() const
        {
            return _isOnBoundingBox;
        }
        
        /**********************************************************************/
        bool isBoundaryNode() const
        {
            bool isBndNode(_isOnBoundingBox);
            if(isBndNode)
            {
                isBndNode=(boundaryNormal.squaredNorm()>FLT_EPSILON); // otherwise isGrainBoundaryNode() must be true
                if(!isBndNode && !isGrainBoundaryNode())
                {
                    model::cout<<"DislocationNode "<<this->sID<<", P="<<this->get_P().transpose()<<std::endl;
                    model::cout<<"_isOnBoundingBox="<<_isOnBoundingBox<<std::endl;
                    model::cout<<"isGrainBoundaryNode()="<<isGrainBoundaryNode()<<std::endl;
                    model::cout<<"NODE ON BoundingBox MUST BE EITHER A BOUNDARY NODE OR A GB NODE"<<std::endl;
                    //                    assert(0 && "NODE ON BoundingBox MUST BE EITHER A BOUNDARY NODE OR A GB NODE");
                }
            }
            return isBndNode ;
        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return grainBoundaries().size();
        }
        
        //        /**********************************************************************/
        //        bool isPureBoundaryNode() const
        //        {
        //            return isBoundaryNode() && isConnectedToBoundaryNodes();
        //        }
        
        /**********************************************************************/
        bool isPureBoundaryNode() const
        {
            bool temp=isBoundaryNode();
            if(temp)
            {
                for (const auto& neighborIter : this->neighbors())
                {
                    temp*=std::get<1>(neighborIter.second)->isBoundarySegment();
                    if(!temp)
                    {
                        break;
                    }
                }
            }
            return  temp;
        }
        
        /**********************************************************************/
        bool isConnectedToBoundaryNodes() const
        {
            bool temp(true);
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=(std::get<0>(neighborIter.second)->isBoundaryNode() || std::get<1>(neighborIter.second)->hasZeroBurgers());
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleBoundaryNode() const
        {
            bool temp=false;
            if(isOnBoundingBox())
            {
                temp=true; // true if all non-virtual neighbors are boundary
                std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> chordDeq;
                
                for (const auto& neighborIter : this->neighbors())
                {
                    if (!std::get<1>(neighborIter.second)->isVirtualBoundarySegment()
                        //                            !std::get<1>(neighborIter.second)->hasZeroBurgers()
                        )
                    {
                        temp*=std::get<1>(neighborIter.second)->isBoundarySegment();
                        chordDeq.push_back(std::get<1>(neighborIter.second)->chord());
                    }
                }
                
                if(temp && chordDeq.size())
                {
                    for(const auto& chord : chordDeq)
                    {
                        temp*=(chord.cross(chordDeq[0]).squaredNorm()<FLT_EPSILON);
                    }
                }
                else
                {
                    temp=false;
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleGrainBoundaryNode() const
        {
            
            bool temp=false;
            if(isGrainBoundaryNode())
            {
                if(this->isSimple())
                {
                    temp=true;
                    std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> chordDeq;
                    
                    for (const auto& neighborIter : this->neighbors())
                    {
                        if (!std::get<1>(neighborIter.second)->hasZeroBurgers())
                        {
                            temp*=std::get<1>(neighborIter.second)->isGrainBoundarySegment();
                            chordDeq.push_back(std::get<1>(neighborIter.second)->chord());
                        }
                    }
                    
                    if(temp && chordDeq.size())
                    {
                        for(const auto& chord : chordDeq)
                        {
                            temp*=(chord.cross(chordDeq[0]).squaredNorm()<(FLT_EPSILON*chord.squaredNorm()*chordDeq[0].squaredNorm()));
                        }
                    }
                    else
                    {
                        temp=false;
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSessileNode() const
        {
            
            bool temp=true;
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=std::get<1>(neighborIter.second)->isSessile();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isSimpleSessileNode() const
        {
            bool temp(this->isSimple() && isSessileNode());
            if(temp)
            {// make sure attached sessile segments are aligned
                const LinkType* firstLink(std::get<1>(this->neighbors().begin()->second));
                const LinkType* secondLink(std::get<1>(this->neighbors().rbegin()->second));
                temp*=(firstLink->chord().normalized().cross(secondLink->chord().normalized()).norm()<FLT_EPSILON);
            }
            return temp;
        }
        
        /**********************************************************************/
        NeighborContainerType nonZeroNeighbors() const
        {
            NeighborContainerType temp;
            for (const auto& neighborIter : this->neighbors())
            {
                if(!std::get<1>(neighborIter.second)->hasZeroBurgers())
                {
                    temp.emplace(neighborIter.first,neighborIter.second);
                }
            }
            return temp;
        }
        
        
        
        /**********************************************************************/
        bool isRemovable(const double& Lmin,const double& cosRemove) const
        {
        
            bool temp(   !isVirtualBoundaryNode
                      && (isSimpleBoundaryNode() || isSimpleGrainBoundaryNode() || isSimpleSessileNode() || isGeometricallyRemovable(Lmin,cosRemove))
                      );
            
//            if(temp && usePeriodic && virtualNode)
//            {
//                    temp*=(virtualNode->isSimpleBoundaryNode() ) ;
//            }
            
            return temp;
        }
        
        
        /**********************************************************************/
        bool isGeometricallyRemovable(const double& Lmin,const double& cosRemove) const
        {
            bool temp=false;
            const auto linksMap=this->linksByLoopID();
            if(linksMap.size()==1)
            {
                const LoopLinkContainerType& linkSet(linksMap.begin()->second);
                assert(linkSet.size()==2);
                const LoopLinkType& link0(**linkSet. begin());
                const LoopLinkType& link1(**linkSet.rbegin());
                
                if(link0.loop()->isGlissile)
                {
                    const VectorDim chord0(link0.sink()->get_P()-link0.source()->get_P());
                    const VectorDim chord1(link1.sink()->get_P()-link1.source()->get_P());
                    const double chord0Norm(chord0.norm());
                    const double chord1Norm(chord1.norm());
                    
                    if(chord0Norm<Lmin || chord1Norm<Lmin)
                    {
                        if(chord0.dot(chord1)>cosRemove*chord0Norm*chord1Norm)
                        {
                            const VectorDim dv0(link0.sink()->get_V()-link0.source()->get_V());
                            const VectorDim dv1(link1.sink()->get_V()-link1.source()->get_V());
                            if(chord0.dot(dv0)<0.0 || chord1.dot(dv1)<0.0) // at least one of the two segments is getting shorter
                            {
                                temp=true;
                            }
                        }
                    }
                    
                    if(chord0Norm+chord1Norm<Lmin)
                    {
                        temp=true;
                    }
                }
            }
            
            return temp;
        }
        
        
        /**********************************************************************/
        const double& velocityReduction() const
        {
            return velocityReductionCoeff;
        }
        
        /**********************************************************************/
        const bool& isGlissile() const
        {
            return _isGlissile;
        }
        
        /**********************************************************************/
        void setToBoundary(const VectorDim& X)
        {
            _isOnBoundingBox=true;
            boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
            
            if(isBoundaryNode() && !isVirtualBoundaryNode) //ENABLE THIS USE_BVP
            {
                
                switch (this->network().simulationParameters.simulationType)
                {
                    case DefectiveCrystalParameters::FINITE_FEM:
                    {
                        if(this->network().useVirtualExternalLoops)
                        {
                            if(virtualNode)
                            {
                                static_cast<NodeBaseType*>(virtualNode.get())->set_P(X+100.0*boundaryNormal);
                            }
                            else
                            {
                                VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
                                virtualNode.reset(new NodeType(&this->network(),X+100.0*boundaryNormal));
                            }
                        }
                        break;
                    }
                        
                    case DefectiveCrystalParameters::PERIODIC:
                    {
                            if(virtualNode)
                            {
                                static_cast<NodeBaseType*>(virtualNode.get())->set_P(X - ((this->network().mesh.xMax()-this->network().mesh.xMin()).cwiseProduct(boundaryNormal)));
                            }
                            else
                            {
                                VerboseDislocationNode(2,"DislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
                                virtualNode.reset(new NodeType(&this->network(),X-((this->network().mesh.xMax()-this->network().mesh.xMin()).cwiseProduct(boundaryNormal))));
                            }
                        break;
                    }
                        
                    default:
                        break;
                }
                

            }
            NodeBaseType::set_P(X); // in turn this calls PlanarDislocationSegment::updateGeometry, so the boundaryNormal must be computed before this line 
            assert(boundingBoxSegments().contains(this->get_P()).first);
        }
        
        /**********************************************************************/
        bool set_P(const VectorDim& newP)
        {
            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" current P="<< this->get_P().transpose()<<"set_P to "<<newP.transpose()<<std::endl;);
            // make sure that node is on glide planes
            bool glidePlanesContained=true;
            for(const auto& gp : meshPlanes())
            {
                glidePlanesContained*=gp->contains(newP);
            }
            
            
            if(glidePlanesContained)
            {
                //                NodeBaseType::set_P(newP);
//                if(this->network().use_boundary) // using confining mesh
//                {
                
                    if(_isOnBoundingBox)
                    {// node was on bounding box, it must remain on bounding box
                        const VectorDim X(snapToBoundingBox(newP));
                        setToBoundary(X);
                        //                        boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                        //                        NodeBaseType::set_P(X);
                        
                        if(addGrainBoundaryPlanes())
                        {// GB-planes were added, the bounding box has changed, so snap again
                            VerboseDislocationNode(3,"case 4"<<std::endl;);
                            const VectorDim X1(snapToBoundingBox(this->get_P()));
                            setToBoundary(X1);
                            //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(X1,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                            //                            NodeBaseType::set_P(X1); // put node back on bounding box
                        }
                    }
                    else
                    {// internal node
                        
                        std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                        
                        if(temp.first)
                        {// internal node, and newP is inside current grain
                            if(   isConnectedToBoundaryNodes()
                               && boundingBoxSegments().size()==2
                               && glidePlaneIntersections().size()==1)
                            {// force special case to boundary to get rid of small debris
                                if((newP-glidePlaneIntersections()[0].first).norm()<this->network().surfaceAttractionDistance)
                                {
                                    setToBoundary(glidePlaneIntersections()[0].first);
                                    //                                    _isOnBoundingBox=true;
                                    //                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(glidePlaneIntersections()[0].first,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                                    //                                    NodeBaseType::set_P(glidePlaneIntersections()[0].first);
                                }
                                else if((newP-glidePlaneIntersections()[0].second).norm()<this->network().surfaceAttractionDistance)
                                {
                                    setToBoundary(glidePlaneIntersections()[0].second);
                                    //                                    _isOnBoundingBox=true;
                                    //                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(glidePlaneIntersections()[0].second,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                                    //                                    NodeBaseType::set_P(glidePlaneIntersections()[0].second);
                                }
                                else
                                {
                                    NodeBaseType::set_P(newP);
                                }
                            }
                            else
                            {
                                NodeBaseType::set_P(newP);
                            }
                        }
                        else
                        {// internal node, and newP is outside current grain
                            
                            const VectorDim X(snapToBoundingBox(newP));
                            //                            _isOnBoundingBox=true;
                            //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                            //                            NodeBaseType::set_P(X); // put node back on bounding box
                            setToBoundary(X);
                            
                            
                            if(addGrainBoundaryPlanes())
                            {// GB-planes were added, the bounding box has changed, and it may now be a set of degenerate lines
                                if(boundingBoxSegments().contains(this->get_P()).first)
                                {// new bounding box contains node
                                    VerboseDislocationNode(3,"case 5"<<std::endl;);
                                    //                                    _isOnBoundingBox=true;
                                    //                                    const VectorDim X(snapToBoundingBox(newP));
                                    //                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                                    //                                    NodeBaseType::set_P(X); // put node back on bounding box
                                    //                                    NodeBaseType::set_P(snapToBoundingBox(newP)); // kill numerical errors
                                }
                                else
                                {// new bounding box does not contain node
                                    VerboseDislocationNode(3,"case 6"<<std::endl;);
                                    const VectorDim X1(snapToMeshPlaneIntersection(this->get_P()));
                                    _isOnBoundingBox=false;
                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(X1,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                                    NodeBaseType::set_P(X1); // kill numerical errors
                                    //                                    setToBoundary(X1);
                                }
                            }
                            //                            else
                            //                            {// node is now on bounding box, and no GB planes were added
                            //                                VerboseDislocationNode(3,"case 7"<<std::endl;);
                            //                            }
                        }
                        
                        
                    }
                    
                    ///////////////
                    
                    //                    std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                    //                    if(temp.first)
                    //                    {// newP is inside current grain
                    //                        if(_isOnBoundingBox)
                    //                        {// if the node is already on the the bounding box, keep it there
                    //                            VerboseDislocationNode(3,"case 1"<<std::endl;);
                    //                            const VectorDim X(snapToBoundingBox(newP));
                    //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                            NodeBaseType::set_P(X);
                    //                        }
                    //                        else
                    //                        {// node was internal to the grain and remains internal
                    //                            VerboseDislocationNode(3,"case 2"<<std::endl;);
                    //                            if(   isConnectedToBoundaryNodes()
                    //                               && boundingBoxSegments().size()==2
                    //                               && glidePlaneIntersections().size()==1)
                    //                            {
                    //                                if((newP-glidePlaneIntersections()[0].first).norm()<this->network().surfaceAttractionDistance)
                    //                                {
                    //                                    NodeBaseType::set_P(glidePlaneIntersections()[0].first);
                    //                                }
                    //                                else if((newP-glidePlaneIntersections()[0].second).norm()<this->network().surfaceAttractionDistance)
                    //                                {
                    //                                    NodeBaseType::set_P(glidePlaneIntersections()[0].second);
                    //                                }
                    //                                else
                    //                                {
                    //                                    NodeBaseType::set_P(newP);
                    //                                }
                    //                            }
                    //                            else
                    //                            {
                    //                                NodeBaseType::set_P(newP);
                    //                            }
                    //
                    //                            if(boundingBoxSegments().contains(this->get_P()).first)
                    //                            {// there is a chance that newP is exactly on the bounding box
                    //                                _isOnBoundingBox=true;
                    //                                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                                NodeBaseType::set_P(this->get_P());
                    //                            }
                    //                        }
                    //                    }
                    //                    else
                    //                    {// newP is outside current grain (or on grain boundary)
                    //                        VerboseDislocationNode(3,"case 3"<<std::endl;);
                    //                        const VectorDim newPbox(snapToBoundingBox(newP));
                    //                        boundaryNormal=SimplexBndNormal::get_boundaryNormal(newPbox,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                        NodeBaseType::set_P(newPbox); // put node on the bouding box
                    //                        if(_isOnBoundingBox)
                    //                        {// node was on bounding box, and exited the grain
                    //                            if(addGrainBoundaryPlanes())
                    //                            {// GB-planes were added, the bounding box has changed, so snap again
                    //                                VerboseDislocationNode(3,"case 4"<<std::endl;);
                    //                                const VectorDim X(snapToBoundingBox(newP));
                    //                                boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                                NodeBaseType::set_P(X); // put node back on bounding box
                    //                            }
                    //                        }
                    //                        else
                    //                        {// node was internal and exited the grain
                    //                            if(addGrainBoundaryPlanes())
                    //                            {// GB-planes were added, the bounding box has changed, and it may now be a set of degenerate lines
                    //                                if(boundingBoxSegments().contains(this->get_P()).first)
                    //                                {// new bounding box contains node
                    //                                    VerboseDislocationNode(3,"case 5"<<std::endl;);
                    //                                    _isOnBoundingBox=true;
                    //                                    const VectorDim X(snapToBoundingBox(newP));
                    //                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                                    NodeBaseType::set_P(X); // put node back on bounding box
                    ////                                    NodeBaseType::set_P(snapToBoundingBox(newP)); // kill numerical errors
                    //                                }
                    //                                else
                    //                                {// new bounding box does not contain node
                    //                                    VerboseDislocationNode(3,"case 6"<<std::endl;);
                    //                                    _isOnBoundingBox=false;
                    //                                    const VectorDim X(snapToMeshPlaneIntersection(newP));
                    //                                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                                    NodeBaseType::set_P(X); // kill numerical errors
                    //                                }
                    //                            }
                    //                            else
                    //                            {// node is now on bounding box, and no GB planes were added
                    //                                VerboseDislocationNode(3,"case 7"<<std::endl;);
                    //                                _isOnBoundingBox=true;
                    //                                const VectorDim X(snapToBoundingBox(newP));
                    //                                boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                    //                                NodeBaseType::set_P(X); // put node back on bounding box
                    //                            }
                    //                        }
                    //
                    //                    }
                    
                    p_Simplex=get_includingSimplex(p_Simplex); // update including simplex
                    //                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // check if node is now on a boundary
                    //                    if(boundaryNormal.squaredNorm()<FLT_EPSILON)
                    //                    {
                    //                        model::cout<<"DislocationNode "<<this->sID<<", @"<<this->get_P().transpose()<<std::endl;
                    //                        assert(false && "BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL");
                    //                    }
                    
                    if(_isOnBoundingBox)
                    {
                        assert(boundingBoxSegments().contains(this->get_P()).first);
                        boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // check if node is now on a boundary
                        if(boundaryNormal.squaredNorm()<FLT_EPSILON && !isGrainBoundaryNode())
                        {
                            model::cout<<"DislocationNode "<<this->sID<<", @"<<std::setprecision(15)<<std::scientific<<this->get_P().transpose()<<std::endl;
                            model::cout<<"BoundingBox Lines:"<<std::endl;
                            model::cout<<boundingBoxSegments();
                            //                            for (const auto& pair : boundingBoxSegments())
                            //                            {
                            //                                model::cout<<std::setprecision(15)<<std::scientific<<"("<<pair.first.transpose()<<"), ("<<pair.second.transpose()<<std::endl;
                            //                            }
                            std::cout<<"BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL"<<std::endl;
                            //assert(false && "BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL");
                        }
                    }
                    else
                    {
                        boundaryNormal.setZero();
                    }
                    
//                }
//                else
//                {
//                    VerboseDislocationNode(3,"case 7"<<std::endl;);
//                    NodeBaseType::set_P(newP);
//                }
                
                
                //                make_projectionMatrix();
            }
            else
            {
                model::cout<<"DislocationNode "<<this->sID<<std::endl;
                assert(0 && "new position outside glide planes.");
            }
            
            
            
            VerboseDislocationNode(3,"DislocationNode "<<this->sID<<" _isOnBoundingBox="<<_isOnBoundingBox<<std::endl;);
            
            return (this->get_P()-newP).norm()<FLT_EPSILON;
        }
        
        /**********************************************************************/
        bool isMovableTo(const VectorDim& X) const
        {
            bool isMovable=true;
            
            VerboseDislocationNode(4,"checking if DislocationNode "<<this->sID<< " isMovable:"<<std::endl;);
            
            for(const auto& gp : meshPlanes())
            {// X must be contained by all glidePlanes
                isMovable*=gp->contains(X);
            }
            VerboseDislocationNode(4,"  meshPlanes contains X? "<<isMovable<<std::endl;);
            
            if(isMovable)
            {
                
                //                HERE CHECK THAT X IS IN REGION
                
                for(const auto& pair : this->neighbors())
                {
                    if(std::get<1>(pair.second)->isBoundarySegment())
                    {// boundary segments other than must remain boundary if this node is moved
                        const bool bndNeighborMovable=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
                        VerboseDislocationNode(4,"  boundaryNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<bndNeighborMovable<<std::endl;);
                        isMovable*=bndNeighborMovable;
                        if(!isMovable)
                        {
                            break;
                        }
                    }
                    
                    if(std::get<1>(pair.second)->isGrainBoundarySegment())
                    {// grain-boundary segments must remain grain-boundary if this node is moved
                        for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
                        {
                            const bool gbNeighborMovable=gb->contains(X);
                            VerboseDislocationNode(4,"  gbNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<gbNeighborMovable<<std::endl;);
                            isMovable*=gbNeighborMovable;
                            if(!isMovable)
                            {
                                break;
                            }
                        }
                        if(!isMovable)
                        {
                            break;
                        }
                    }
                    
                    if(std::get<1>(pair.second)->isSessile())
                    {// sessile segments cannot change direction if this node is moved
                        const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*(std::get<0>(pair.second)->get_P()-X).norm()*(std::get<0>(pair.second)->get_P()-this->get_P()).norm());
                        VerboseDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
                        isMovable*=sessileNeighborMovable;
                        if(!isMovable)
                        {
                            break;
                        }
                        //                        isMovable*=LineSegment<dim>(std::get<0>(pair.second)->get_P(),X).contains(this->get_P());
                        //                        isMovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
                    }
                }
            }
            
            return isMovable;
        }
        
        /**********************************************************************/
        void move(const double & dt,const double& dxMax)
        {
            
            //velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
            
            
            
            
            
            const VectorDim P_old(this->get_P());
            
            //VectorDim dX=2.0*A2-A1-this->get_P();
            
            VectorDim dX=velocity.template segment<dim>(0)*dt;
            //VectorDim dX=velocity.template segment<dim>(0)*dt*velocityReductionCoeff;
            //VectorDim dX=(1.5*velocity.template segment<dim>(0)-0.5*vOld)*dt; // Adams–Bashforth
            //VectorDim dX=(0.5*velocity.template segment<dim>(0)+0.5*vOld)*dt; // Adams–Bashforth
            
            
            //            //Limit dX for boundaryNodes bec
            //            const double dXnorm(dX.norm());
            //            if((isBoundaryNode() || isConnectedToBoundaryNodes()) && dXnorm>dxMax)
            //            {
            //                dX*=dxMax/dXnorm;
            //            }
            
            if (dX.squaredNorm()>0.0 && _isGlissile) // move a node only if |v|!=0
            {
                
                // Make sure that new position is at intersection of glidePlanes
                const VectorDim newP=snapToMeshPlaneIntersection(this->get_P()+dX);
                set_P(newP);
            }
            else
            {
                velocity.setZero();
            }
            
            // Store actual velocity
            if(dt>0.0)
            {
                velocity=(this->get_P()-P_old)/dt;
            }
            
            //            vOld=velocity; // store current value of velocity before updating
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const NodeType& ds)
        {
            os<< DislocationNodeIO<dim>(ds);
            return os;
        }
        
    };
    
    
    // static data
    template <int _dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<_dim,corder,InterpolationType>::use_velocityFilter=true;
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    double DislocationNode<_dim,corder,InterpolationType>::velocityReductionFactor=0.75;
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    const double DislocationNode<_dim,corder,InterpolationType>::bndTol=1.0e-4;
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    int DislocationNode<_dim,corder,InterpolationType>::verboseDislocationNode=0;
    
}
#endif


//        /**********************************************************************/
//        void applyVelocityFilter(double vMaxGood)
//        {
//            if(use_velocityFilter)
//            {
//                const double Rc=3.0*10.0;
//                const double filterThreshold=0.05*vMaxGood*vMaxGood;
//
//                if((this->get_P()-C).norm()<Rc)
//                {
//                    if(velocity.dot(vOld)<-filterThreshold)
//                    {
//                        velocityReductionCoeff*=velocityReductionFactor;
//                    }
//                }
//                else
//                {
//                    if(velocity.dot(vOld)<-filterThreshold)
//                    {
//                        velocityReductionCoeff*=velocityReductionFactor;
//                        C=this->get_P();
//                    }
//                    else if(velocity.dot(vOld)>0.0)
//                    {
//                        velocityReductionCoeff/=velocityReductionFactor;
//                    }
//                    else
//                    {
//                        // don't change velocityReductionCoeff
//                    }
//                }
//
//                //                const double filterThreshold=0.05*velocity.norm()*vOld.norm();
//
//                //                const double VdotVold=
//                //                if(velocity.dot(vOld)<-filterThreshold)
//                //                {
//                //                    velocityReductionCoeff*=velocityReductionFactor;
//                //                }
//                //                else if(velocity.dot(vOld)>filterThreshold)
//                //                {
//                //                    velocityReductionCoeff/=velocityReductionFactor;
//                //                }
//                //                else
//                //                {
//                //                    // don't change velocityReductionCoeff
//                //                }
//                if(velocityReductionCoeff>1.0)
//                {
//                    velocityReductionCoeff=1.0;
//                }
//                if(velocityReductionCoeff<0.0001)
//                {
//                    velocityReductionCoeff=0.0001;
//                }
//
//
//
//                if( isOscillating()
//                   //&& velocity.norm()>vMaxGood /*velocity.norm()>0.0*/
//                   )
//                { // oscillating node
//                    std::cout<<"node "<<this->sID<<" BAD, ";
//                    //                    velocity*=(velocityReductionCoeff*vMaxGood/velocity.norm());
//                    if(velocity.norm()>vMaxGood)
//                    {
//                        velocity=velocityReductionCoeff*vMaxGood*velocity.normalized();
//                    }
//                    //                    std::cout<<velocity.norm()<<std::endl;
//                }
//                else
//                { // good node
//                    std::cout<<"node "<<this->sID<<" GOOD, ";
//                    velocity*=velocityReductionCoeff;
//                }
//
//                //                velocity*=velocityReductionCoeff;
//                std::cout<<"velocityReductionCoeff="<<velocityReductionCoeff<<", vMaxGood="<<vMaxGood<<", velocity.norm()="<<velocity.norm()<<", newVelocity="<<velocity.norm()<<std::endl;
//
//            }
//        }

//        /**********************************************************************/
//        bool is_simple() const
//        {
//            size_t nonZeroLink=0;
//            for (const auto& neighborIter : this->neighbors())
//            {
////                if (!std::get<2>(neighborIter.second)==0)
////                {
//                    if (!std::get<1>(neighborIter.second)->hasZeroBurgers())
//                    {  // neighbor not searched
//                        nonZeroLink++;
//                    }
////                }
//            }
//            return (nonZeroLink==2);
//        }

//        /**********************************************************************/
//        const BoundingLineSegments<dim>& glidePlaneIntersections() const
//        {
//            return glidePlaneIntersections();
//        }

//        /**********************************************************************/
//        const BoundingLineSegments<dim>& boundingBoxSegments() const
//        {
//            return boundingBoxSegments();
//        }

//
//        /**********************************************************************/
//        bool isPureGBNode() const
//        {
//            bool temp(!this->is_isolated());
//            for (const auto& neighborIter : this->neighbors())
//            {
////                if (std::get<2>(neighborIter.second)) // not self
////                {
//                    temp*=std::get<0>(neighborIter.second)->isGrainBoundaryNode();
////                }
//            }
//
//            return (isGrainBoundaryNode()&&temp);
//        }

//        /******************************************************************************/
//        void neighborsAt(const LatticeVectorType& L0, std::set<size_t>& temp) const
//        {/*!\param[in] P0 position to be serached
//          * \param[out]temp set of IDs of neighbors of this which are located at P0 (possibly including *this)
//          * \param[in] tol tolerance used to detect position overlap
//          */
//            for (const auto& nIiter : this->neighbors())
//            { // loop over neighborhood
//                if((std::get<0>(nIiter.second)->get_L()-L0).squaredNorm()==0)
//                { // a neighbor of I exists at P0
//                    temp.insert(std::get<0>(nIiter.second)->sID);
//                }
//            }
//        }
