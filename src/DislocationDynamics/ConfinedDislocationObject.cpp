/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConfinedDislocationObject_cpp_
#define model_ConfinedDislocationObject_cpp_

#include <ConfinedDislocationObject.h>


//#include <BoundingFiniteLineSegments.h>

namespace model
{

        /**********************************************************************/
        template <int dim>
        void ConfinedDislocationObject<dim>::updateBoundingBoxWithMeshFace(const PlanarMeshFace<dim>& face)
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
        
        /**********************************************************************/
        template <int dim>
        ConfinedDislocationObject<dim>::ConfinedDislocationObject(const PositionCointainerType& temp) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateConfinement(temp);
        }
        
        /**********************************************************************/
        template <int dim>
        ConfinedDislocationObject<dim>::ConfinedDislocationObject(const VectorDim& P0) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateConfinement(P0);
        }
        
        /**********************************************************************/
        template <int dim>
        ConfinedDislocationObject<dim>::ConfinedDislocationObject(const VectorDim& P0,const VectorDim& P1) :
        //        ConfinedDislocationObject(GlidePlaneObserver<dim>& gpo) :
        //        /* init */ gpObserver(gpo)
        /* init */ _isOnExternalBoundary(false)
        /* init */,_isOnInternalBoundary(false)
        /* init */,_outNormal(VectorDim::Zero())
        {
            updateConfinement(P0,P1);
        }
        
        
        /**********************************************************************/
        template <int dim>
        ConfinedDislocationObject<dim>::ConfinedDislocationObject(const ConfinedDislocationObject<dim>& A,
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
        
        template <int dim>
        ConfinedDislocationObject<dim>& ConfinedDislocationObject<dim>::confinedObject()
        {
            return *this;
        }

        template <int dim>
        const ConfinedDislocationObject<dim>& ConfinedDislocationObject<dim>::confinedObject() const
        {
            return *this;
        }

        
        /**********************************************************************/
        template <int dim>
        void ConfinedDislocationObject<dim>::clear()
        {
            glidePlanes().clear();
//            meshFaces().clear(); // NEVER CLEAR CONFINING FACES !!!
            _glidePlaneIntersections.reset(nullptr);
            this->boundingBoxSegments().clear();
        }
        
        /**********************************************************************/
        template <int dim>
        const typename ConfinedDislocationObject<dim>::GlidePlaneContainerType& ConfinedDislocationObject<dim>::glidePlanes() const
        {
            return *this;
        }
        
        template <int dim>
        typename ConfinedDislocationObject<dim>::GlidePlaneContainerType& ConfinedDislocationObject<dim>::glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        template <int dim>
        const typename ConfinedDislocationObject<dim>::PlanarMeshFaceContainerType& ConfinedDislocationObject<dim>::meshFaces() const
        {
            return *this;
        }
        
        template <int dim>
        typename ConfinedDislocationObject<dim>::PlanarMeshFaceContainerType& ConfinedDislocationObject<dim>::meshFaces()
        {
            return *this;
        }
        
//        const GlidePlane<dim>* glidePlane() const
//        {
//            return glidePlanes().size()==1? *glidePlanes().begin() : nullptr;
//        }
        
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
        template <int dim>
        std::set<size_t> ConfinedDislocationObject<dim>::meshFaceIDs() const
        {
            std::set<size_t> temp;
            for(const auto& face : this->meshFaces())
            {
                temp.insert(face->sID);
            }
            return temp;
        }
        
        /**********************************************************************/
        template <int dim>
        bool ConfinedDislocationObject<dim>::isOnMeshFaces(const std::set<size_t>& faceIDs) const
        {
            std::set<size_t> theseFaceIDs(meshFaceIDs());
            bool temp(true);
            for(const auto& val : faceIDs)
            {
                temp= (temp && (theseFaceIDs.find(val)!=theseFaceIDs.end()));
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        template <int dim>
        typename ConfinedDislocationObject<dim>::VectorDim ConfinedDislocationObject<dim>::snapToGlidePlanes(const VectorDim& P)
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
        template <int dim>
        typename ConfinedDislocationObject<dim>::VectorDim ConfinedDislocationObject<dim>::snapToGlidePlanesinPeriodic(const VectorDim& P)
        {/*!@param[in] P input positions
          *\returns the position P snapped to the (intersection of) glide planes
          */
            if(_glidePlaneIntersections)
            {
                //Add here forsnapping to point
                if ((_glidePlaneIntersections->P0-_glidePlaneIntersections->P1).norm()>FLT_EPSILON)
                {
                    return _glidePlaneIntersections->snapToInfiniteLine(P);
                }
                else
                {
                    //Snap to the point
                    return 0.5*(_glidePlaneIntersections->P0+_glidePlaneIntersections->P1);
                }
            }
            else
            {
                //                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToGlidePlanes, case 0"<<std::endl;);
                // std::cout<<" glidePlanes size "<<this->glidePlanes().size()<<std::endl;
                assert(this->glidePlanes().size()==1);
                return (*this->glidePlanes().begin())->snapToPlane(P);
            }
            
        }
        
        /**********************************************************************/
        template <int dim>
        const std::unique_ptr<FiniteLineSegment<dim>>& ConfinedDislocationObject<dim>::glidePlaneIntersections() const
        {
            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        template <int dim>
        const bool& ConfinedDislocationObject<dim>::isOnExternalBoundary() const
        {/*!\returns _isOnExternalBoundarySegment.
          */
            return _isOnExternalBoundary;
        }
        
        /**********************************************************************/
        template <int dim>
        const bool& ConfinedDislocationObject<dim>::isOnInternalBoundary() const
        {
            return _isOnInternalBoundary;
        }
        
        /**********************************************************************/
        template <int dim>
        bool ConfinedDislocationObject<dim>::isOnBoundary() const
        {
            return _isOnExternalBoundary || _isOnInternalBoundary;
        }
        
        /**********************************************************************/
        template <int dim>
        const typename ConfinedDislocationObject<dim>::VectorDim& ConfinedDislocationObject<dim>::bndNormal() const
        {
            return _outNormal;
        }
        
        /**********************************************************************/
        template <int dim>
        void ConfinedDislocationObject<dim>::updateConfinement(const VectorDim& P0)
        {
            posCointainer.clear();
            posCointainer.push_back(P0);
            updateConfinement();
        }
        
        /**********************************************************************/
        template <int dim>
        void ConfinedDislocationObject<dim>::updateConfinement(const VectorDim& P0,const VectorDim& P1)
        {
            posCointainer.clear();
            posCointainer.push_back(P0);
            posCointainer.push_back(P1);
            updateConfinement();
        }
        
        template <int dim>
        void ConfinedDislocationObject<dim>::updateConfinement(const PositionCointainerType& temp)
        {
            posCointainer=temp;
            updateConfinement();
        }
        
        /**********************************************************************/
        template <int dim>
        void ConfinedDislocationObject<dim>::addGlidePlane(const GlidePlaneType* const lastGlidePlane)
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
//                            const PlanePlaneIntersection<dim> ppi(**glidePlanes().begin(),**glidePlanes().rbegin());
                            const GlidePlane<dim>& glidePlane0(**glidePlanes().begin());
                            const GlidePlane<dim>& glidePlane1(**glidePlanes().rbegin());
                            const PlanePlaneIntersection<dim> ppi(glidePlane0,glidePlane1);

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
//                                        std::cout<<"incident u0="<<u0<<std::endl;
//                                        std::cout<<"segLength="<<segLength<<std::endl;
                                        if(u0>=0.0 && u0<=segLength)
                                        {
                                            roots.insert(lli.x0);
                                        }
                                    }
                                    else if(lli.type==LineLineIntersection<dim>::COINCIDENT)
                                    {// a coincident line was found, which means that the glide planes intersec on a boundary face
//                                        std::cout<<"coincident"<<std::endl;
                                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(meshInt->P0,meshInt->P1));
                                        this->boundingBoxSegments().push_back(meshInt);
                                        break;
                                    }
                                }
                                
                                if(!_glidePlaneIntersections)
                                {// no coincident intersection was found
                                    if(roots.size()!=2)
                                    {
                                        if(posCointainer.size())
                                        {
                                            std::cout<<"Plane0 bounding box:"<<std::endl;
                                            std::cout<<glidePlane0.meshIntersections<<std::endl;
                                            
                                            std::cout<<"Plane1 bounding box:"<<std::endl;
                                            std::cout<<glidePlane1.meshIntersections<<std::endl;
                                            
                                            
                                            std::cout<<"BOUNDARY POINTS OF INCIDENT PLANES ARE:"<<std::endl;
                                            for(const auto& root : roots)
                                            {
                                                std::cout<<root.transpose()<<std::endl;
                                            }
                                            assert(false && "THERE MUST BE 2 INTERSECTION POINTS BETWEEN GLIDEPLANE(s) and GRAIN-BOUNDARY PERIMETER");
                                        }
                                    }
                                    else
                                    {// two roots
                                        _glidePlaneIntersections.reset(new FiniteLineSegment<dim>(*roots.begin(),*roots.rbegin()));
                                    }
                                    
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
                            {// parallel planes. _glidePlaneIntersections remians null and boundingBoxSegments is empty
                                if(posCointainer.size())
                                {
                                    assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                                }
                            }
                            
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
                                            //                                    std::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                                            std::cout<<"MeshPlanes are:"<<std::endl;
                                            for(const auto& plane : glidePlanes())
                                            {
                                                std::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
                                                std::cout<<"Plane bounding box:"<<std::endl;
                                                std::cout<<plane->meshIntersections<<std::endl;
                                                
                                            }
                                            
                                            std::cout<<"lastGlidePlane is:"<<std::endl;
                                            std::cout<<std::setprecision(15)<<std::scientific<<"  P="<<lastGlidePlane->P.transpose()<<", n="<<lastGlidePlane->unitNormal.transpose()<<std::endl;
                                            
                                            std::cout<<"MeshPlane intersection is:"<<std::endl;
                                            std::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P1="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
                                            
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
        template <int dim>
        void ConfinedDislocationObject<dim>::updateConfinement()
        {
            for(const auto& glidePlane : glidePlanes())
            {
                
                for(const auto& pos : posCointainer)
                {
                    if (!glidePlane->contains(pos))
                    {
                        std::cout<<"GlidePlane size is  "<<glidePlanes().size()<<std::endl;
                        std::cout<<"pos is "<<pos.transpose()<<std::endl;
                        std::cout<<"Position different is "<<(pos-(glidePlane->snapToPlane(pos))).squaredNorm()<<"for glide plane size is "<<glidePlanes().size()<<std::endl;
                    }
                    assert(glidePlane->contains(pos) && "glidePlane MUST CONTAIN POSITION");
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
                            cointained=(cointained && face.second->asPlane().contains(pos));
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
                    if(!face->asPlane().contains(pos))
                    {
                        throw std::runtime_error("FACE MUS CONTAIN POSITION");
                    }
//                    assert(face->asPlane().contains(pos) && "FACE MUS CONTAIN POSITION");
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
        template <int dim>
        typename ConfinedDislocationObject<dim>::GrainContainerType ConfinedDislocationObject<dim>::grains() const
        {
            GrainContainerType temp;
            for(const auto& glidePlane : glidePlanes())
            {
                temp.insert(&glidePlane->grain);
            }
            return temp;
        }
        
        /**********************************************************************/
        template <int dim>
        std::vector<std::pair<const GlidePlane<dim>* const,const GlidePlane<dim>* const>> ConfinedDislocationObject<dim>::parallelAndCoincidentGlidePlanes(const GlidePlaneContainerType& other) const
        {
            std::vector<std::pair<const GlidePlane<dim>* const,const GlidePlane<dim>* const>> pp;
            
            for(const auto& plane : glidePlanes())
            {
                for(const auto& otherPlane : other)
                {
//                    if(plane!=otherPlane && plane->n.cross(otherPlane->n).squaredNorm()==0)
                    if(plane->n.cross(otherPlane->n).squaredNorm()==0)
                    {// parallel planes
                        pp.emplace_back(plane,otherPlane);
                    }
                }
            }
            return pp;
        }
        template struct ConfinedDislocationObject<3>;
}
#endif
