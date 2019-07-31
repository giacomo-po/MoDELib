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
#include <GlidePlane.h>
#include <Grain.h>
#include <PlanarMeshFace.h>
#include <FiniteLineSegment.h>
#include <LineLineIntersection.h>

//#include <BoundingFiniteLineSegments.h>

namespace model
{
    
    
    
    template <int dim>
    struct ConfinedDislocationObject : private std::set<const GlidePlane<dim>*>
    /*                              */,private std::set<const PlanarMeshFace<dim>*>
    /*                              */,public BoundingMeshSegments<dim>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType*> GrainContainerType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType*> PlanarMeshFaceContainerType;
        typedef MeshBoundarySegment<dim> MeshBoundarySegmentType;
        typedef std::vector<VectorDim> PositionCointainerType;
        
        
    private:
        
        /**********************************************************************/
        void updateBoundingBoxWithMeshFace(const PlanarMeshFace<dim>& face)
        {
            //            std::cout<<"pdateBoundingBoxWithMeshFace "<<face.sID<<std::endl;
            //            std::cout<<this->boundingBoxSegments()<<std::endl;
            if(this->boundingBoxSegments().size())
            {// A bounding box already exists
//                assert(this->boundingBoxSegments().size() && "boundingBoxSegments() cannot be empty");
                
                BoundingMeshSegments<dim> temp;
                
                Plane<dim> plane(face.asPlane());
                for(const auto& seg : this->boundingBoxSegments())
                {// PRBLEM HERE IS THAT THIS ALGORITHM DOES NOT CHECK THAT TEMP ALREADY ONCLUDE "OVERLAPPING" BOUNDING BOX POINTS
                    auto newFaces(seg->faces);
                    newFaces.insert(&face);
                    
                    //                if(seg.faces.find(&face)!=seg.faces.end())
                    //                {// intersection segment is defined on face
                    const bool P0contained(plane.contains(seg->P0));
                    const bool P1contained(plane.contains(seg->P1));
                    //                    auto faces(seg.faces());
                    //                    faces.insert(face);
                    if(P0contained && P1contained)
                    {// the whole segment is contained
                        temp.emplace_back(new MeshBoundarySegment<dim>(seg->P0,seg->P1,newFaces));
                        //                        const bool success(.second);
                        //                        assert(success && "COULD NOT INSERT IN TEMP DURING FACE UPDATE");
                    }
                    else if(P0contained)
                    {
                        temp.emplace_back(new MeshBoundarySegment<dim>(seg->P0,seg->P0,newFaces));
                        //                        assert(false && "FINISH HERE1");
                    }
                    else if(P1contained)
                    {
                        temp.emplace_back(new MeshBoundarySegment<dim>(seg->P1,seg->P1,newFaces));
                        //                        assert(false && "FINISH HERE2");
                    }
                    else
                    {// do nothing, this excludes seg from new bounding box
                        //                                               assert(false && "FINISH HERE3");
                    }
                    //                }
                }
                
                this->boundingBoxSegments().swap(temp);
                //            std::cout<<"Now bounding box is"<<std::endl;
                //            std::cout<<this->boundingBoxSegments()<<std::endl;
                assert(this->boundingBoxSegments().size() && "boundingBoxSegments cannot be empty");
            }

        }
        

        
        
        
        //        GlidePlaneObserver<dim>& gpObserver;
        PositionCointainerType posCointainer;
        std::unique_ptr<FiniteLineSegment<dim>> _glidePlaneIntersections;
        bool _isOnExternalBoundary;
        bool _isOnInternalBoundary;
        VectorDim _outNormal;
        
        
    public:
        
        /**********************************************************************/
        ConfinedDislocationObject(const PositionCointainerType& temp) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateGeometry(temp);
        }
        
        /**********************************************************************/
        ConfinedDislocationObject(const VectorDim& P0) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateGeometry(P0);
        }
        
        /**********************************************************************/
        ConfinedDislocationObject(const VectorDim& P0,const VectorDim& P1) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateGeometry(P0,P1);
        }
        
        
        /**********************************************************************/
        ConfinedDislocationObject(const ConfinedDislocationObject<dim>& A,
                                  const ConfinedDislocationObject<dim>& B) :
        //        /* init */ gpObserver(A.gpObserver)
        /* init */ _isOnExternalBoundary(A.isOnExternalBoundary() || B.isOnExternalBoundary())
        /* init */,_isOnInternalBoundary(A.isOnInternalBoundary() || B.isOnInternalBoundary())
        /* init */,_outNormal(VectorDim::Zero())
        {
            for(const auto& glidePlane : A.glidePlanes())
            {
                addGlidePlane(glidePlane);
            }
            for(const auto& glidePlane : B.glidePlanes())
            {
                addGlidePlane(glidePlane);
            }
            for(const auto& face : A.meshFaces())
            {
                meshFaces().insert(face);
            }
            for(const auto& face : B.meshFaces())
            {
                meshFaces().insert(face);
            }
            
        }
        
        /**********************************************************************/
        void clear()
        {
            glidePlanes().clear();
//            meshFaces().clear(); // NEVER CLEAR CONFINING FACES !!!
            _glidePlaneIntersections.reset(nullptr);
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
//        std::set<size_t> imageMeshFaceIDs(const std::set<size_t>& mirroringFaceIDs) const
//        {
//            const auto faces(imageMeshFaces(mirroringFaceIDs));
//            std::set<size_t> temp;
//            for(const auto& face : faces)
//            {
//                temp.insert(face->sID);
//            }
//            return temp;
//        }
        
        /**********************************************************************/
        std::set<size_t> meshFaceIDs() const
        {
            std::set<size_t> temp;
            for(const auto& face : this->meshFaces())
            {
                temp.insert(face->sID);
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isOnMeshFaces(const std::set<size_t>& faceIDs) const
        {
            std::set<size_t> theseFaceIDs(meshFaceIDs());
            bool temp(true);
            for(const auto& val : faceIDs)
            {
                temp*=(theseFaceIDs.find(val)!=theseFaceIDs.end());
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        VectorDim snapToGlidePlanes(const VectorDim& P)
        {/*!@param[in] P input positions
          *\returns the position P snapped to the (intersection of) glide planes
          */
            if(_glidePlaneIntersections)
            {
                return _glidePlaneIntersections->snap(P);
            }
            else
            {
                //                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToGlidePlanes, case 0"<<std::endl;);
                assert(this->glidePlanes().size()==1);
                return (*this->glidePlanes().begin())->snapToPlane(P);
            }
            
        }
        
        /**********************************************************************/
        const std::unique_ptr<FiniteLineSegment<dim>>& glidePlaneIntersections() const
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
        void updateGeometry(const VectorDim& P0)
        {
            posCointainer.clear();
            posCointainer.push_back(P0);
            updateConfinement();
        }
        
        /**********************************************************************/
        void updateGeometry(const VectorDim& P0,const VectorDim& P1)
        {
            posCointainer.clear();
            posCointainer.push_back(P0);
            posCointainer.push_back(P1);
            updateConfinement();
        }
        
        /**********************************************************************/
        void addGlidePlane(const GlidePlaneType* const lastGlidePlane)
        {
            if(lastGlidePlane)
            {// a glidePlane exists
                const bool success(glidePlanes().insert(lastGlidePlane).second);
                if(success)
                {// A new glide plane was added
                    //                    updateBoundingBoxWithGlidePlane(*glidePlane);
                    //                    updateMeshPlaneIntersections(*glidePlane);
                    
                    
                    switch (glidePlanes().size())
                    {
                        case 0:
                        {// there must be at least one glide plane
                            _glidePlaneIntersections.reset(nullptr);
                            this->boundingBoxSegments().clear();
                            assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                            break;
                        }
                            
                        case 1:
                        {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                            //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
                            _glidePlaneIntersections.reset(nullptr);
                            this->boundingBoxSegments().clear();
                            for(const auto& seg : lastGlidePlane->meshIntersections)
                            {// copy boundary segments from plane
                                this->boundingBoxSegments().push_back(seg);
                            }
                            break;
                        }
                            
                        case 2:
                        {// a second plane is being added, so we must have no _glidePlaneIntersections
                            //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
                            //                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                            _glidePlaneIntersections.reset(nullptr);
                            this->boundingBoxSegments().clear();
                            //
                            //                            assert(!_glidePlaneIntersections && "_glidePlaneIntersections must be empty");
                            
                            // Grab the infinite line of intersection between the two planes
                            //                            const PlanePlaneIntersection<dim>& ppi(gpObserver.glidePlaneIntersection(*glidePlanes().begin(),*glidePlanes().rbegin()));
                            const PlanePlaneIntersection<dim> ppi(**glidePlanes().begin(),**glidePlanes().rbegin());
                            const GlidePlane<dim>& glidePlane0(**glidePlanes().begin());
                            const GlidePlane<dim>& glidePlane1(**glidePlanes().rbegin());
                            
                            if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                            {/* Two distinct glide planes can be coincident only if they belong to different grains
                              */
                                
                                const Grain<dim>& grain0(glidePlane0.grain);
                                const Grain<dim>& grain1(glidePlane1.grain);
                                
                                
                                assert(grain0.grainID!=grain1.grainID);
                                const GrainBoundary<dim>& gb(*grain0.grainBoundaries().at(std::make_pair(grain0.grainID,grain1.grainID)));
                                
                                //                                std::vector<VectorDim> roots;
                                std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > roots;
                                for(const auto& meshInt : gb.meshIntersections)
                                {
                                    PlaneSegmentIntersection<dim> psi(glidePlane0,*meshInt);
                                    if(psi.type==PlaneSegmentIntersection<dim>::INCIDENT)
                                    {
                                        //                                        roots.push_back(psi.x0);
                                        roots.insert(psi.x0);
                                    }
                                }
                                assert(roots.size()==2 && "THERE MUST BE 2 INTERSECTION POINTS BETWEEN GLIDEPLANE(s) and GRAIN-BOUNDARY PERIMETER");
                                _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(*roots.begin(),*roots.rbegin()));
                                this->boundingBoxSegments().emplace_back(new MeshBoundarySegment<dim>(*roots.begin(),*roots.rbegin(),gb.face.get()));
                            }
                            else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                            {/* If the two planes are incident then the intersection of
                              * their bounding boxes is either a pair of singluar segments (2 points)
                              * or a line segment on the boundary
                              */
                                
                                std::set<Eigen::Matrix<double,dim,1>,CompareVectorsByComponent<double,dim,float> > roots;
                                
                                //                                std::vector<VectorDim> roots;
                                for(const auto& meshInt : glidePlane0.meshIntersections)
                                {
                                    const double segLength((meshInt->P1-meshInt->P0).norm());
                                    const VectorDim D0((meshInt->P1-meshInt->P0)/segLength);
                                    LineLineIntersection<dim> lli(meshInt->P0,D0,ppi.P,ppi.d);
                                    if(lli.type==LineLineIntersection<dim>::INCIDENT)
                                    {
                                        const double u0((lli.x0-meshInt->P0).dot(D0));
                                        //                                        std::cout<<u0<<std::endl;
                                        if(u0>=0.0 && u0<=segLength)
                                        {
                                            //                                            roots.push_back(lli.x0);
                                            roots.insert(lli.x0);
                                        }
                                    }
                                    else if(lli.type==LineLineIntersection<dim>::COINCIDENT)
                                    {// a coincident line was found, which means that the glide planes intersec on a boundary face
                                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(meshInt->P0,meshInt->P1));
                                        this->boundingBoxSegments().push_back(meshInt);
                                        break;
                                    }
                                }
                                
                                if(!_glidePlaneIntersections)
                                {// no coincident intersection was found
                                    if(roots.size()!=2)
                                    {
                                        model::cout<<"BOUNDARY POINTS OF INCIDENT PLANES ARE:"<<std::endl;
                                        for(const auto& root : roots)
                                        {
                                            model::cout<<root.transpose()<<std::endl;
                                        }
                                        assert(false && "THERE MUST BE 2 INTERSECTION POINTS BETWEEN GLIDEPLANE(s) and GRAIN-BOUNDARY PERIMETER");
                                    }
                                    //                                    _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(roots[0],roots[1]));
                                    _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(*roots.begin(),*roots.rbegin()));
                                    
                                    for(const auto& root : roots)
                                    {
                                        std::set<const PlanarMeshFace<dim>*> faces;
                                        const auto containingSegments0(glidePlane0.meshIntersections.containingSegments(root));
                                        const auto containingSegments1(glidePlane1.meshIntersections.containingSegments(root));
                                        assert(containingSegments0.size() && "glidePlane0 must contain root");
                                        assert(containingSegments1.size() && "glidePlane1 must contain root");
                                        
                                        for(const auto& meshInt : containingSegments0)
                                        {
                                            for(const auto& curFace : meshInt->faces)
                                            {
                                                faces.insert(curFace);
                                            }
                                        }
                                        for(const auto& meshInt : containingSegments1)
                                        {
                                            for(const auto& curFace : meshInt->faces)
                                            {
                                                faces.insert(curFace);
                                            }
                                        }
                                        
                                        this->boundingBoxSegments().emplace_back(new MeshBoundarySegment<dim>(root,root,faces));
                                    }
                                }
                            }
                            else
                            {
                                assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                            }
                            
                            //                            assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
                            
                            break;
                        }
                            
                        default:
                        {// Case of more that 2 planes.
                            if(_glidePlaneIntersections)
                            {// A _glidePlaneIntersections must exist, otherwise interseciton of glide planes was already found empty
                                
                                
                                PlaneSegmentIntersection<dim> pli(lastGlidePlane->P,
                                                                  lastGlidePlane->unitNormal,
                                                                  _glidePlaneIntersections->P0, // segment start
                                                                  _glidePlaneIntersections->P1 // segment end
                                                                  );
                                
                                
                                switch (pli.type)
                                {
                                    case PlaneSegmentIntersection<dim>::COINCIDENT:
                                    {// nothing to do, _glidePlaneIntersections and bounding box remains unchanged
                                        break;
                                    }
                                        
                                    case PlaneSegmentIntersection<dim>::INCIDENT:
                                    {// _glidePlaneIntersections becomes a point (degenerate line)
                                        const VectorDim x(0.5*(pli.x0+pli.x1));
                                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(x,x));
                                        
                                        // Update this->boundingBoxSegments()
                                        std::set<const PlanarMeshFace<dim>*> faces;
                                        const auto containingSegments(this->boundingBoxSegments().containingSegments(x));
                                        for(const auto& meshInt : containingSegments)
                                        {
                                            for(const auto& curFace : meshInt->faces)
                                            {
                                                faces.insert(curFace);
                                            }
                                        }
                                        this->boundingBoxSegments().clear();
                                        if(faces.size())
                                        {
                                            this->boundingBoxSegments().emplace_back(new MeshBoundarySegment<dim>(x,x,faces));
                                        }
                                        break;
                                    }
                                        
                                    default:
                                    {
                                        _glidePlaneIntersections.reset(nullptr);
                                        this->boundingBoxSegments().clear();
                                        
                                        if(posCointainer.size())
                                        {// an actual dislocation object is being confined
                                            //                                    model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                                            model::cout<<"MeshPlanes are:"<<std::endl;
                                            for(const auto& plane : glidePlanes())
                                            {
                                                model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
                                                model::cout<<"Plane bounding box:"<<std::endl;
                                                std::cout<<plane->meshIntersections<<std::endl;
                                                
                                            }
                                            
                                            model::cout<<"lastGlidePlane is:"<<std::endl;
                                            model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<lastGlidePlane->P.transpose()<<", n="<<lastGlidePlane->unitNormal.transpose()<<std::endl;
                                            
                                            model::cout<<"MeshPlane intersection is:"<<std::endl;
                                            model::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P1="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
                                            
                                            assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                                        }
                                        break;
                                    }
                                }
                                break;
                            }
                            else
                            {
                                this->boundingBoxSegments().clear();
                            }
                        }
                            
                    }
                    
                    updateConfinement();
                }
            }
        }
        
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
                    //                    if(meshFaces().find(face.second.get())!=meshFaces().end())
                    //                    {// face is already a current confining face
                    //                        for(const auto& pos : posCointainer)
                    //                        {
                    //                            assert(face.second->asPlane().contains(pos) && "FACE MUS CONTAIN POSITION");
                    //                        }
                    //                    }
                    if(meshFaces().find(face.second.get())==meshFaces().end())
                    {// face not a current confining face
                        bool cointained(posCointainer.size()); // if posCointainer is empty set cointained to false
                        for(const auto& pos : posCointainer)
                        {
                            cointained*=face.second->asPlane().contains(pos);
                        }
                        if(cointained)
                        {// faces contains all positions
                            meshFaces().insert(face.second.get());
                            //                            updateBoundingBoxWithMeshFace(*face.second);
                        }
                    }
                }
            }
            
            
            //            updateBoundingBoxWithMeshFaces();
            
            
            
            // Update _isOnExternalBoundary, _isOnInternalBoundary, and _outNormal
            _isOnExternalBoundary=false;
            _isOnInternalBoundary=false;
            _outNormal.setZero();
            for(const auto& face : meshFaces())
            {
                
                for(const auto& pos : posCointainer)
                {// A face must include all positions
                    assert(face->asPlane().contains(pos) && "FACE MUS CONTAIN POSITION");
                }
                
                updateBoundingBoxWithMeshFace(*face);
                
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


//        /**********************************************************************/
//        void updateMeshPlaneIntersections(const GlidePlaneType& lastGlidePlane)
//        {
//            //            BoundingFiniteLineSegments<dim> temp;
//
//            //            //VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections"<<std::endl;);
//            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.P="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.P.transpose()<<std::endl;);
//            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.unitNormal="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.unitNormal.transpose()<<std::endl;);
//
//            switch (glidePlanes().size())
//            {
//                case 0:
//                {// there must be at least one glide plane
//                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
//                    break;
//                }
//
//                case 1:
//                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
//                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
//                    _glidePlaneIntersections.reset(nullptr);
//                    break;
//                }
//
//                case 2:
//                {// a second plane is being added, so we must have no _glidePlaneIntersections
//                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
//                    //                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
//                    assert(!_glidePlaneIntersections && "_glidePlaneIntersections must be empty");
//
//                    // Grab the infinite line of intersection between the two planes
//                    const PlanePlaneIntersection<dim>& ppi(this->network().glidePlaneIntersection(*glidePlanes().begin(),*glidePlanes().rbegin()));
//
//                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
//                    {/* Two distinct glide planes can be coincident only if they belong to different grains
//                      * In that case, the intersection of their bounding boxes should be one line segment
//                      */
//                        //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2a"<<std::endl;);
//                        if(this->boundingBoxSegments().size()!=1)
//                        {
//                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
//                            model::cout<<"glidePlane(0) is "<<*glidePlanes().begin()->P.transpose()<<","<<*glidePlanes().begin()->unitNormal.transpose()<<std::endl;
//                            model::cout<<"glidePlane(1) is "<<*glidePlanes().rbegin()->P.transpose()<<","<<*glidePlanes().rbegin()->unitNormal.transpose()<<std::endl;
//                            assert(false && "There should be only one line in boundingBoxSegments()");
//                        }
//                        //assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
//                        //                        _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
//                        //                        _glidePlaneIntersections.push_back(boundingBoxSegments()[0]);
//                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second));
//                    }
//                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
//                    {/* If the two planes are incident then the intersection of
//                      * their bounding boxes is either a pair of singluar segments (2 points)
//                      * or a line segment on the boundary
//                      */
//                        switch (this->boundingBoxSegments().size())
//                        {
//                            case 1:
//                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
//                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2b"<<std::endl;);
//                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
//                                _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().begin()->second.P1));
//                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P0="<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;);
//                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P1="<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;);
//
//                                break;
//                            }
//
//                            case 2:
//                            {// The two intersections must be degenerate (2 boundary points)
//                                //                                std::cout<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P0.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P1.transpose()<<std::endl;
//                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2c"<<std::endl;);
//                                assert((this->boundingBoxSegments(). begin()->second.P0-this->boundingBoxSegments(). begin()->second.P1).squaredNorm()<FLT_EPSILON);
//                                assert((this->boundingBoxSegments().rbegin()->second.P0-this->boundingBoxSegments().rbegin()->second.P1).squaredNorm()<FLT_EPSILON);
//                                _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().rbegin()->second.P0));
//                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().rbegin()->second.P0);
//                                break;
//                            }
//
//                            default:
//                            {
//                                model::cout<<"PlanarDislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
//                                model::cout<<this->boundingBoxSegments();
//                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
//                            }
//                        }
//                    }
//                    else
//                    {
//                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                    }
//
//                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
//
//                    break;
//                }
//
//                default:
//                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
//                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
//
//
//                    PlaneSegmentIntersection<dim> pli(lastGlidePlane.P,
//                                                      lastGlidePlane.unitNormal,
//                                                      _glidePlaneIntersections->P0, // origin of line
//                                                      _glidePlaneIntersections->P1 // line direction
//                                                      );
//
//
//                    switch (pli.type)
//                    {
//                        case PlaneSegmentIntersection<dim>::COINCIDENT:
//                        {// nothing to do, _glidePlaneIntersections remains unchanged
//                            break;
//                        }
//
//                        case PlaneSegmentIntersection<dim>::INCIDENT:
//                        {// _glidePlaneIntersections becomes a point (degenerate line)
//                            const VectorDim x(0.5*(pli.x0+pli.x1));
//                            _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(x,x));
//                            break;
//                        }
//
//                        default:
//                        {
//                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
//                            model::cout<<"MeshPlanes are:"<<std::endl;
//                            for(const auto& plane : glidePlanes())
//                            {
//                                model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
//                            }
//                            model::cout<<"MeshPlane intersection is:"<<std::endl;
//                            model::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P2="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
//
//                            assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                            break;
//                        }
//                    }
//                }
//
//            }
//
//            if(_glidePlaneIntersections)
//            {
//                //VerbosePlanarDislocationNode(2,"  _glidePlaneIntersections are: "<<_glidePlaneIntersections->P0.transpose()<<", "<<_glidePlaneIntersections->P1.transpose()<<std::endl;);
//            }
//        }

//        /**********************************************************************/
//        void updateMeshPlaneIntersections(const GlidePlaneType& lastGlidePlane)
//        {
//            //            BoundingFiniteLineSegments<dim> temp;
//
//            //            //VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections"<<std::endl;);
//            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.P="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.P.transpose()<<std::endl;);
//            //VerbosePlanarDislocationNode(2,"  lastGlidePlane.unitNormal="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.unitNormal.transpose()<<std::endl;);
//
//            switch (glidePlanes().size())
//            {
//                case 0:
//                {// there must be at least one glide plane
//                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
//                    break;
//                }
//
//                case 1:
//                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
//                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
//                    _glidePlaneIntersections.reset(nullptr);
//                    break;
//                }
//
//                case 2:
//                {// a second plane is being added, so we must have no _glidePlaneIntersections
//                    //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
//                    //                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
//                    assert(!_glidePlaneIntersections && "_glidePlaneIntersections must be empty");
//
//                    // Grab the infinite line of intersection between the two planes
//                    const PlanePlaneIntersection<dim>& ppi(this->network().glidePlaneIntersection(*glidePlanes().begin(),*glidePlanes().rbegin()));
//
//                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
//                    {/* Two distinct glide planes can be coincident only if they belong to different grains
//                      * In that case, the intersection of their bounding boxes should be one line segment
//                      */
//                        //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2a"<<std::endl;);
//                        if(this->boundingBoxSegments().size()!=1)
//                        {
//                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
//                            model::cout<<"glidePlane(0) is "<<*glidePlanes().begin()->P.transpose()<<","<<*glidePlanes().begin()->unitNormal.transpose()<<std::endl;
//                            model::cout<<"glidePlane(1) is "<<*glidePlanes().rbegin()->P.transpose()<<","<<*glidePlanes().rbegin()->unitNormal.transpose()<<std::endl;
//                            assert(false && "There should be only one line in boundingBoxSegments()");
//                        }
//                        //assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
//                        //                        _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
//                        //                        _glidePlaneIntersections.push_back(boundingBoxSegments()[0]);
//                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second));
//                    }
//                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
//                    {/* If the two planes are incident then the intersection of
//                      * their bounding boxes is either a pair of singluar segments (2 points)
//                      * or a line segment on the boundary
//                      */
//                        switch (this->boundingBoxSegments().size())
//                        {
//                            case 1:
//                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
//                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2b"<<std::endl;);
//                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
//                                _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().begin()->second.P1));
//                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P0="<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;);
//                                //VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P1="<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;);
//
//                                break;
//                            }
//
//                            case 2:
//                            {// The two intersections must be degenerate (2 boundary points)
//                                //                                std::cout<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P0.transpose()<<std::endl;
//                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P1.transpose()<<std::endl;
//                                //VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2c"<<std::endl;);
//                                assert((this->boundingBoxSegments(). begin()->second.P0-this->boundingBoxSegments(). begin()->second.P1).squaredNorm()<FLT_EPSILON);
//                                assert((this->boundingBoxSegments().rbegin()->second.P0-this->boundingBoxSegments().rbegin()->second.P1).squaredNorm()<FLT_EPSILON);
//                                _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(this->boundingBoxSegments().begin()->second.P0,this->boundingBoxSegments().rbegin()->second.P0));
//                                //                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().rbegin()->second.P0);
//                                break;
//                            }
//
//                            default:
//                            {
//                                model::cout<<"PlanarDislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
//                                model::cout<<this->boundingBoxSegments();
//                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
//                            }
//                        }
//                    }
//                    else
//                    {
//                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                    }
//
//                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
//
//                    break;
//                }
//
//                default:
//                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
//                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");
//
//
//                    PlaneSegmentIntersection<dim> pli(lastGlidePlane.P,
//                                                      lastGlidePlane.unitNormal,
//                                                      _glidePlaneIntersections->P0, // origin of line
//                                                      _glidePlaneIntersections->P1 // line direction
//                                                      );
//
//
//                    switch (pli.type)
//                    {
//                        case PlaneSegmentIntersection<dim>::COINCIDENT:
//                        {// nothing to do, _glidePlaneIntersections remains unchanged
//                            break;
//                        }
//
//                        case PlaneSegmentIntersection<dim>::INCIDENT:
//                        {// _glidePlaneIntersections becomes a point (degenerate line)
//                            const VectorDim x(0.5*(pli.x0+pli.x1));
//                            _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(x,x));
//                            break;
//                        }
//
//                        default:
//                        {
//                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
//                            model::cout<<"MeshPlanes are:"<<std::endl;
//                            for(const auto& plane : glidePlanes())
//                            {
//                                model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
//                            }
//                            model::cout<<"MeshPlane intersection is:"<<std::endl;
//                            model::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P2="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
//
//                            assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                            break;
//                        }
//                    }
//                }
//
//            }
//
//            if(_glidePlaneIntersections)
//            {
//                //VerbosePlanarDislocationNode(2,"  _glidePlaneIntersections are: "<<_glidePlaneIntersections->P0.transpose()<<", "<<_glidePlaneIntersections->P1.transpose()<<std::endl;);
//            }
//        }
