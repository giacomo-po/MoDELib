/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfinedDislocationObject_H_
#define model_ConfinedDislocationObject_H_

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <GlidePlane.h>
#include <Grain.h>
#include <PlanarMeshFace.h>
#include <LineSegment.h>
//#include <BoundingLineSegments.h>

namespace model
{
    
    template <int dim>
    struct BoundingMeshSegments : public MeshPlane<dim>::MeshBoundarySegmentContainerType
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;

        
        BoundingMeshSegments(){}
        
        BoundingMeshSegments(const BoundingMeshSegments<dim>& bb1,const BoundingMeshSegments<dim>& bb2)
        {
            assert(0 && "FINISH HERE");
        }
        
        /**********************************************************************/
        std::set<const MeshBoundarySegment<dim>*> containingSegments(const VectorDim& P) const
        {
            std::set<const MeshBoundarySegment<dim>*> temp;
            
            for(const auto& seg : *this)
            {
                //                std::cout<<pair.first<<", d="<<pair.second.distanceTo(P)<<std::endl;
                if(seg.contains(P))
                {
                    temp.insert(&seg);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool contains(const VectorDim& P) const
        {
            return containingSegments(P).size();
        }
        
        /**********************************************************************/
        const BoundingMeshSegments<dim>& boundingBoxSegments() const
        {
            return *this;
        }
        
        BoundingMeshSegments<dim>& boundingBoxSegments()
        {
            return *this;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const BoundingMeshSegments<dim>& bls)
        {
            for(const auto& seg : bls)
            {
                os<<seg<<std::endl;
            }
            return os;
        }
        
    };
    
    
    template <int dim>
    struct ConfinedDislocationObject : private std::set<const GlidePlane<dim>*>
    /*                              */,private std::set<const PlanarMeshFace<dim>*>
        /*                              */,public BoundingMeshSegments<dim>
    //    /*                              */,private MeshPlane<dim>::MeshBoundarySegmentContainerType
//    /*                              */,private BoundingLineSegments<dim>
//    /*                              */,private std::map<std::set<const PlanarMeshFace<dim>*>,MeshBoundarySegment<dim>>
//    /*                              */,private std::set<MeshBoundarySegment<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType*> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType*> PlanarMeshFaceContainerType;
        
        typedef MeshBoundarySegment<dim> MeshBoundarySegmentType;
//        typedef std::map<std::set<const PlanarMeshFaceType*>,MeshBoundarySegmentType> MeshBoundarySegmentContainerType;
//        typedef BoundingMeshSegments<dim> MeshBoundarySegmentContainerType;
        
        typedef std::vector<VectorDim, Eigen::aligned_allocator<VectorDim>> PositionCointainerType;
        
        
    private:
        
        /**********************************************************************/
        void updateBoundingBoxWithMeshFace(const PlanarMeshFace<dim>& face)
        {
            assert(this->boundingBoxSegments().size() && "boundingBoxSegments() cannot be empty");

            BoundingMeshSegments<dim> temp;

            Plane<dim> plane(face.asPlane());
            for(const auto& seg : this->boundingBoxSegments())
            {
                if(seg.find(&face)!=seg.end())
                {// intersection segment is defined on face
                    const bool P0contained(plane.contains(seg.P0));
                    const bool P1contained(plane.contains(seg.P1));
//                    auto faces(seg.faces());
//                    faces.insert(face);
                    if(P0contained && P1contained)
                    {// the whole segment is contained
                        const bool success(temp.inser(seg).second);
                        assert(success && "COULD NOT INSERT IN TEMP DURING FACE UPDATE");
                    }
                    else if(P0contained)
                    {
                        assert(false && "FINISH HERE1");
                    }
                    else if(P1contained)
                    {
                        assert(false && "FINISH HERE2");
                    }
                    else
                    {
                                                assert(false && "FINISH HERE3");
                    }
                }
                
                
                
//                const MeshBoundarySegment<dim>& seg(pair.second);
//                const bool P0contained(plane.contains(seg.P0));

            }

            this->swap(temp);
            assert(this->boundingBoxSegments().size() && "boundingBoxSegments cannot be empty");
        }

        /**********************************************************************/
        void updateBoundingBoxWithGlidePlane(const GlidePlane<dim>& lastGlidePlane)
        {
            if(glidePlanes().size())
            {
                assert(0 && "FINISH HERE");
            }
            else
            {
                for(const auto& seg : lastGlidePlane.meshIntersections)
                {
                    assert(!seg.hasZeroLength() && "Plane-Face intersection has zero length");
                    const bool success(this->boundingBoxSegments().insert(seg).second);
                    assert(success && "COULD NOT INSERT MeshBoundarySegment from first GlidePlane");
                }
            }
        }
        
        
        
//        {
//
//
////            boundingBoxSegments().clear();
////
////            int planeID(0);
////            for(const auto& glidePlane : glidePlanes())
////            {
////                if(planeID==0)
////                {//first plane
////                    for(const auto: seg : lastGlidePlane.meshIntersections)
////                    {
////                        const bool success(boundingBoxSegments().insert(seg).second);
////                        assert(success && "COULD NOT INSERT MeshBoundarySegment from first GlidePlane");
////                    }
////                }
////                else
////                {
////
////                }
////                planeID++;
////            }
//
//
//
//            if(glidePlanes().size())
//            {
//
//                MeshBoundarySegmentContainerType newBox;
//
////                std::vector<MeshBoundarySegmentType,Eigen::aligned_allocator<MeshBoundarySegmentType>> temp;
//
////                MeshBoundarySegmentContainerType temp;
//
////                MeshBoundarySegmentContainerType lastBox;
////                for(const auto: seg : lastGlidePlane.meshIntersections)
////                {
////                    assert(!seg.hasZeroLength() && "Plane-Face intersection has zero length");
////                    const bool success(lastBox.insert(seg).second);
////                    assert(success && "COULD NOT INSERT MeshBoundarySegment from lastGlidePlane");
////                }
//
//                for(const auto& s1 : lastGlidePlane.meshIntersections())
//                {// loop over all boundary segments of lastGlidePlane
//
//                    bool foundIntersection=false;
////                    std::unique_ptr<MeshBoundarySegment<dim>> newSeg;
//                    for(const auto& s2 : boundingBoxSegments())
//                    {// for each s1 we need to keep the "best" intersection with s2
//
//                        std::set<const PlanarMeshFaceType*> commonFaces;
//                        std::set_intersection(s1.faces().begin(),s1.faces().end(),s2.faces().begin(),s2.faces.endend(),std::inserter(commonFaces,commonFaces.begin()));
//
//                        if(commonFaces.size())
//                        {// the two segments have at least a common face
//                            SegmentSegmentDistance<dim> ssd(segCopy.P0,segCopy.P1,
//                                                            s2.P0,s2.P1);
//
//                            std::set<const PlanarMeshFaceType*> unionFaces;
//                            std::set_union(s1.faces().begin(),s1.faces().end(),s2.faces().begin(),s2.faces.endend(),std::inserter(unionFaces,unionFaces.begin()));
//
//
//                            const auto iSeg(ssd.intersectionSegment());
//                            std::unique_ptr<MeshBoundarySegment<dim>> newSeg;
//
//                            switch (iSeg.size())
//                            {
//
//                                default:
//                                {// no intersection, do nothing
//                                    break;
//                                }
//
//                                case 1:
//                                {// single point intersection
////                                    auto faces(s1.faces());
////                                    faces.insert(s2.faces().begin(), s2.faces().end());
//                                    newSeg.reset(new MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]),unionFaces));
//                                }
//
//                                case 2:
//                                {// single point intersection
////                                    auto faces(s1.faces());
////                                    faces.insert(s2.faces().begin(), s2.faces().end());
//                                    newSeg.reset(new MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]),unionFaces));
//                                }
//
//
//                            }
//
//                            if(newSeg)
//                            {// an intersection has been found
//
//                            }
//
//                        }
//
//
//
//
//                    }
//
////                    if(s1.hasZeroLength())
////                    {
////
////                    }
////                    else
////                    {// s1 is an extended segment defined on one or more faces. Note that lastBox cannot have degenerate lines (points)
////
////                    }
//
////                    const MeshBoundarySegment<dim>& s1(pair.second);
//
//
////                    if(iter!=lastBox.end())
////                    {// a segment in  lastBox on the same faces has been found
////
////                        const MeshBoundarySegment<dim>& s2(*iter);
////
////
//                        SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,
//                                                        s2.P0,s2.P1);
////
//                        const auto iSeg(ssd.intersectionSegment());
////
////                        switch (iSeg.size())
////                        {
////
////                            default:
////                            {// no intersection, do nothing
////                                break;
////                            }
////
////                            case 1:
////                            {// single-point intersection
////
////                                PROBLEM HERE IS THAT THIS POINT COULD ALREADY BELONG TO ANOTHER FACE
////
////                                for(const auto existingSeg : temp)
////                                {
////                                    if(existingSeg.hasZeroLength())
////                                    {
////                                        if(existingSeg.contains(iSeg[0]))
////                                        {
////                                            auto faces(existingSeg.faces());
////                                            const auto success(temp.erase(faces).second);
////                                            assert(success && "COULD NOT ERASE FROM TEMP");
////
////                                            for(const auto& face : s1)
////                                            {
////                                                faces.insert(face);
////                                            }
////
////                                            const auto succes2(temp.emplace(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]),faces).second);
////                                            assert(success2 && "COULD NOT INSERT MeshBoundarySegment in temp");
////
////                                        }
////                                    }
////                                }
////
////                                const auto success(temp.emplace(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]),s1.faces()).second);
////                                assert(success && "COULD NOT INSERT MeshBoundarySegment in temp");
////                                break;
////                            }
////
////                            case 2:
////                            {// line intersection
////
////                                const auto success(temp.emplace(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]),s1.faces()).second);
////                                assert(success && "COULD NOT INSERT MeshBoundarySegment in temp");
////                                break;
////                            }
////                        }
////                    }
//                }
//
//                this->swap(temp);
//            }
//            else
//            {
//                for(const auto: seg : lastGlidePlane.meshIntersections)
//                {
//                    assert(!seg.hasZeroLength() && "Plane-Face intersection has zero length");
//                    const bool success(boundingBoxSegments().insert(seg).second);
//                    assert(success && "COULD NOT INSERT MeshBoundarySegment from first GlidePlane");
//                }
//            }
//
//            assert(boundingBoxSegments().size() && "boundingBoxSegments cannot be empty");
//
//        }
        
//        /**********************************************************************/
//        void updateBoundingBoxWithGlidePlane(const GlidePlane<dim>& lastGlidePlane)
//        {
//            if(glidePlanes().size())
//            {
//                MeshBoundarySegmentContainerType temp;
//
//                MeshBoundarySegmentContainerType lastBox;
//                for(const auto: seg : lastGlidePlane.meshIntersections)
//                {
//                    const bool success(lastBox.emplace(seg.face->sID,seg).second);
//                    assert(success && "COULD NOT INSERT MeshBoundarySegment in lastBox");
//                }
//
//                for(const auto& pair : boundingBoxSegments())
//                {
//
//                    const MeshBoundarySegment<dim>& s1(pair.second);
//
//                    const auto iter(lastBox.find(pair.first));
//                    if(iter!=lastBox.end())
//                    {// a segment on temp on the same face has been found
//
//                        const MeshBoundarySegment<dim>& s2(iter->second);
//
//
//                        SegmentSegmentDistance<dim> ssd(s1.P0,s1.P1,
//                                                        s2.P0,s2.P1);
//
//                        const auto iSeg(ssd.intersectionSegment());
//
//                        switch (iSeg.size())
//                        {
//
//                            default:
//                            {// no intersection, do nothing
//                                break;
//                            }
//
//                            case 1:
//                            {// single-point intersection
//                                const auto success(temp.emplace(s1.face->sID,MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[0]),s1.face )));
//                                assert(success && "COULD NOT INSERT MeshBoundarySegment in temp");
//                                break;
//                            }
//
//                            case 2:
//                            {// single-point intersection
//                                const auto success(temp.emplace(s1.face->sID,MeshBoundarySegment<dim>(std::get<0>(iSeg[0]),std::get<0>(iSeg[1]),s1.face )));
//                                assert(success && "COULD NOT INSERT MeshBoundarySegment in temp");
//                                break;
//                            }
//                        }
//                    }
//                }
//
//                this->swap(temp);
//            }
//            else
//            {
//                for(const auto: seg : lastGlidePlane.meshIntersections)
//                {
//                    const bool success(boundingBoxSegments().emplace(seg.face->sID,seg).second);
//                    assert(success && "COULD NOT INSERT MeshBoundarySegment in boundingBoxSegments");
//                }
//            }
//
//            assert(boundingBoxSegments().size() && "boundingBoxSegments cannot be empty");
//
//        }
        
        /**********************************************************************/
        void updateConfinement()
        {
            for(const auto& glidePlane : glidePlanes())
            {
                
                for(const auto& pos : posCointainer)
                {
                    assert(glidePlane->contains(pos) && "glidePlane MUS CONTAIN POSITION");
                }
                
                for(const auto& face : glidePlane->grain.region.faces())
                {
                    if(meshFaces().find(face.second.get())!=meshFaces().end())
                    {// face is already a current confining face
                        for(const auto& pos : posCointainer)
                        {
                            assert(face.second->asPlane().contains(pos) && "FACE MUS CONTAIN POSITION");
                        }
                    }
                    else
                    {// face not a current confining face
                        bool cointained(true);
                        for(const auto& pos : posCointainer)
                        {
                            cointained*=face.second->asPlane().contains(pos);
                        }
                        if(cointained)
                        {// faces contains all positions
                            meshFaces().insert(face.second.get());
                            updateBoundingBoxWithMeshFace(*face.second);
                        }
                    }
                }
            }
            
            
            // Update _isOnExternalBoundary, _isOnInternalBoundary, and _outNormal
            _isOnExternalBoundary=false;
            _isOnInternalBoundary=false;
            _outNormal.setZero();
            for(const auto& face : meshFaces())
            {
                if(face->regionIDs.first==face->regionIDs.second)
                {
                    _isOnExternalBoundary=true;
                }
                else
                {
                    _isOnInternalBoundary=true;
                }
                _outNormal+=face->outNormal();
            }
            const double _outNormalNorm(_outNormal.norm());
            if(_outNormalNorm>FLT_EPSILON)
            {
                _outNormal/=_outNormalNorm;
            }
            else
            {
                _outNormal.setZero();
            }
        }
        
        /**********************************************************************/
        void updateMeshPlaneIntersections(const GlidePlaneType& lastGlidePlane)
        {
            //            BoundingLineSegments<dim> temp;
            
            //            //VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections"<<std::endl;);
            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.P="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.P.transpose()<<std::endl;);
            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.unitNormal="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.unitNormal.transpose()<<std::endl;);
            
            switch (glidePlanes().size())
            {
                case 0:
                {// there must be at least one glide plane
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
                    _glidePlaneIntersections.reset(nullptr);
                    break;
                }
                    
                case 2:
                {// a second plane is being added, so we must have no _glidePlaneIntersections
                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
                    //                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                    assert(!_glidePlaneIntersections && "_glidePlaneIntersections must be empty");
                    
                    // Grab the infinite line of intersection between the two planes
                    const PlanePlaneIntersection<dim>& ppi(this->network().glidePlaneIntersection(*glidePlanes().begin(),*glidePlanes().rbegin()));
                    
                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                    {/* Two distinct glide planes can be coincident only if they belong to different grains
                      * In that case, the intersection of their bounding boxes should be one line segment
                      */
                        //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2a"<<std::endl;);
                        if(this->boundingBoxSegments().size()!=1)
                        {
                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                            model::cout<<"glidePlane(0) is "<<*glidePlanes().begin()->P.transpose()<<","<<*glidePlanes().begin()->unitNormal.transpose()<<std::endl;
                            model::cout<<"glidePlane(1) is "<<*glidePlanes().rbegin()->P.transpose()<<","<<*glidePlanes().rbegin()->unitNormal.transpose()<<std::endl;
                            assert(false && "There should be only one line in boundingBoxSegments()");
                        }
                        //assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
                        //                        _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
                        //                        _glidePlaneIntersections.push_back(boundingBoxSegments()[0]);
                        _glidePlaneIntersections.reset(new LineSegment<dim>(this->boundingBoxSegments().begin()->second));
                    }
                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                    {/* If the two planes are incident then the intersection of
                      * their bounding boxes is either a pair of singluar segments (2 points)
                      * or a line segment on the boundary
                      */
                        switch (this->boundingBoxSegments().size())
                        {
                            case 1:
                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2b"<<std::endl;);
                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
                                _glidePlaneIntersections.reset(new LineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().begin()->second.P1));
                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P0="<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;);
                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P1="<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;);
                                
                                break;
                            }
                                
                            case 2:
                            {// The two intersections must be degenerate (2 boundary points)
                                //                                std::cout<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P0.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P1.transpose()<<std::endl;
                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2c"<<std::endl;);
                                assert((this->boundingBoxSegments(). begin()->second.P0-this->boundingBoxSegments(). begin()->second.P1).squaredNorm()<FLT_EPSILON);
                                assert((this->boundingBoxSegments().rbegin()->second.P0-this->boundingBoxSegments().rbegin()->second.P1).squaredNorm()<FLT_EPSILON);
                                _glidePlaneIntersections.reset(new LineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().rbegin()->second.P0));
                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().rbegin()->second.P0);
                                break;
                            }
                                
                            default:
                            {
                                model::cout<<"PlanarDislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
                                model::cout<<this->boundingBoxSegments();
                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
                            }
                        }
                    }
                    else
                    {
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
                    
                    break;
                }
                    
                default:
                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
                    
                    
                    PlaneSegmentIntersection<dim> pli(lastGlidePlane.P,
                                                      lastGlidePlane.unitNormal,
                                                      _glidePlaneIntersections->P0, // origin of line
                                                      _glidePlaneIntersections->P1 // line direction
                                                      );
                    
                    
                    switch (pli.type)
                    {
                        case PlaneSegmentIntersection<dim>::COINCIDENT:
                        {// nothing to do, _glidePlaneIntersections remains unchanged
                            break;
                        }
                            
                        case PlaneSegmentIntersection<dim>::INCIDENT:
                        {// _glidePlaneIntersections becomes a point (degenerate line)
                            const VectorDim x(0.5*(pli.x0+pli.x1));
                            _glidePlaneIntersections.reset(new LineSegment<dim>(x,x));
                            break;
                        }
                            
                        default:
                        {
                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                            model::cout<<"MeshPlanes are:"<<std::endl;
                            for(const auto& plane : glidePlanes())
                            {
                                model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
                            }
                            model::cout<<"MeshPlane intersection is:"<<std::endl;
                            model::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P2="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
                            
                            assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                            break;
                        }
                    }
                }
                    
            }
            
            if(_glidePlaneIntersections)
            {
                //VerbosePlanarDislocationNode(2,"  _glidePlaneIntersections are: "<<_glidePlaneIntersections->P0.transpose()<<", "<<_glidePlaneIntersections->P1.transpose()<<std::endl;);
            }
        }
        
        GlidePlaneObserver<dim>& gpObserver;
        PositionCointainerType posCointainer;
        std::unique_ptr<LineSegment<dim>> _glidePlaneIntersections;
        bool _isOnExternalBoundary;
        bool _isOnInternalBoundary;
        VectorDim _outNormal;
        
        
    public:
        
        /**********************************************************************/
        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        /* init */ gpObserver(gpo)
        /* init */,_isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            
        }
        
        /**********************************************************************/
        void clear()
        {
            glidePlanes().clear();
            meshFaces().clear();
            this->boundingBoxSegments().clear();
        }
        
        /**********************************************************************/
        const GlidePlaneContainerType& glidePlanes() const
        {
            return *this;
        }
        
        GlidePlaneContainerType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const PlanarMeshFaceContainerType& meshFaces() const
        {
            return *this;
        }
        
        PlanarMeshFaceContainerType& meshFaces()
        {
            return *this;
        }
        
//        /**********************************************************************/
//        const BoundingMeshSegments<dim>& boundingBoxSegments() const
//        {
//            return *this;
//        }
//
//        BoundingMeshSegments<dim>& boundingBoxSegments()
//        {
//            return *this;
//        }
        
        /**********************************************************************/
        const std::unique_ptr<LineSegment<dim>>& glidePlaneIntersections() const
        {
            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        const bool& isOnExternalBoundary() const
        {/*!\returns _isOnExternalBoundarySegment.
          */
            return _isOnExternalBoundary;
        }
        
        /**********************************************************************/
        const bool& isOnInternalBoundary() const
        {
            return _isOnInternalBoundary;
        }
        
        /**********************************************************************/
        bool isOnBoundary() const
        {
            return _isOnExternalBoundary || _isOnInternalBoundary;
        }
        
        /**********************************************************************/
        const VectorDim& bndNormal() const
        {
            return _outNormal;
        }
        
        /**********************************************************************/
        void updateGeometry(PositionCointainerType&& temp)
        {
            posCointainer=temp;
            updateConfinement();
        }
        
        /**********************************************************************/
        void addGlidePlane(const GlidePlaneType* const glidePlane)
        {
            if(glidePlane)
            {// a glidePlane exists
                const bool success(glidePlanes().insert(glidePlane).second);
                if(success)
                {// A new glide plane was added
                    updateBoundingBoxWithGlidePlane(*glidePlane);
                    updateMeshPlaneIntersections(*glidePlane);
                    updateConfinement();
                }
            }
        }
        
        /**********************************************************************/
        GrainContainerType grains() const
        {
            GrainContainerType temp;
            for(const auto& glidePlane : glidePlanes())
            {
                temp.insert(&glidePlane->grain);
            }
            return temp;
        }
        

        
        
    };
}
#endif


//        /**********************************************************************/
//        const GrainContainerType& grains() const
//        {
//            return *this;
//        }
//
//        GrainContainerType& grains()
//        {
//            return *this;
//        }
