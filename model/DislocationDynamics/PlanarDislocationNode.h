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

#ifndef model_PLANARDISLOCATIONNODE_H_
#define model_PLANARDISLOCATIONNODE_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <tuple>

#include <DislocationNetworkTraits.h>
#include <SplineNode.h>
#include <GramSchmidt.h>
#include <Simplex.h>
#include <LatticeMath.h>
#include <SimplexBndNormal.h>
#include <Grain.h>
#include <GlidePlane.h>
#include <PlanePlaneIntersection.h>
#include <PlaneLineIntersection.h>
#include <LineSegment.h>
#include <DislocationNodeIO.h>
#include <BoundingLineSegments.h>
#include <DefectiveCrystalParameters.h>


#ifndef NDEBUG
#define VerbosePlanarDislocationNode(N,x) if(verbosePlanarDislocationNode>=N){std::cout<<x;}
#else
#define VerbosePlanarDislocationNode(N,x)
#endif

namespace model
{
    
    template <typename Derived,typename InterpolationType>
    class PlanarDislocationNode :
    /*          */ public SplineNode<Derived,TypeTraits<Derived>::dim,TypeTraits<Derived>::corder,InterpolationType>
    /*          */,private std::set<const GrainBoundary<TypeTraits<Derived>::dim>*>
    /*          */,private std::set<const Grain<TypeTraits<Derived>::dim>*>
    /*          */,private std::set<const MeshPlane<TypeTraits<Derived>::dim>*>
    /*          */,private BoundingLineSegments<TypeTraits<Derived>::dim>
    {
        
    public:
        
        constexpr static int dim=TypeTraits<Derived>::dim; // make dim available outside class
        constexpr static int corder=TypeTraits<Derived>::corder; // make dim available outside class
        typedef Derived NodeType;
        typedef typename TypeTraits<Derived>::LinkType LinkType;
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
        static int verbosePlanarDislocationNode;
        
    private:
        
        /**********************************************************************/
        void updateMeshPlaneIntersections(const MeshPlaneType& lastGlidePlane)
        {
//            BoundingLineSegments<dim> temp;
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections"<<std::endl;);
            VerbosePlanarDislocationNode(2,"  lastGlidePlane.P="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.P.transpose()<<std::endl;);
            VerbosePlanarDislocationNode(2,"  lastGlidePlane.unitNormal="<<std::setprecision(15)<<std::scientific<<lastGlidePlane.unitNormal.transpose()<<std::endl;);
            
            switch (meshPlanes().size())
            {
                case 0:
                {// there must be at least one glide plane
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                    VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 1"<<std::endl;);
                    _glidePlaneIntersections.reset(nullptr);
                    break;
                }
                    
                case 2:
                {// a second plane is being added, so we must have no _glidePlaneIntersections
                    VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2"<<std::endl;);
//                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                    assert(!_glidePlaneIntersections && "_glidePlaneIntersections must be empty");

                    // Grab the infinite line of intersection between the two planes
                    const PlanePlaneIntersection<dim>& ppi(this->network().glidePlaneIntersection(&meshPlane(0),&meshPlane(1)));
                    
                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                    {/* Two distinct glide planes can be coincident only if they belong to different grains
                      * In that case, the intersection of their bounding boxes should be one line segment
                      */
                        VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2a"<<std::endl;);
                        if(boundingBoxSegments().size()!=1)
                        {
                            model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                            model::cout<<"glidePlane(0) is "<<meshPlane(0).P.transpose()<<","<<meshPlane(0).unitNormal.transpose()<<std::endl;
                            model::cout<<"glidePlane(1) is "<<meshPlane(1).P.transpose()<<","<<meshPlane(1).unitNormal.transpose()<<std::endl;
                            assert(false && "There should be only one line in boundingBoxSegments()");
                        }
                        //assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
//                        _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
//                        _glidePlaneIntersections.push_back(boundingBoxSegments()[0]);
                        _glidePlaneIntersections.reset(new LineSegment<dim>(boundingBoxSegments().begin()->second));
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
                                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2b"<<std::endl;);
//                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1);
                                _glidePlaneIntersections.reset(new LineSegment<dim>(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().begin()->second.P1));
                                VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P0="<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;);
                                VerbosePlanarDislocationNode(4,"  boundingBoxSegments().begin()->second.P1="<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;);
                                
                                break;
                            }
                                
                            case 2:
                            {// The two intersections must be degenerate (2 boundary points)
                                //                                std::cout<<boundingBoxSegments().begin()->second.P0.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().begin()->second.P1.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P0.transpose()<<std::endl;
                                //                                std::cout<<boundingBoxSegments().rbegin()->second.P1.transpose()<<std::endl;
                                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" updateMeshPlaneIntersections, case 2c"<<std::endl;);
                                assert((boundingBoxSegments(). begin()->second.P0-boundingBoxSegments(). begin()->second.P1).squaredNorm()<FLT_EPSILON);
                                assert((boundingBoxSegments().rbegin()->second.P0-boundingBoxSegments().rbegin()->second.P1).squaredNorm()<FLT_EPSILON);
                                _glidePlaneIntersections.reset(new LineSegment<dim>(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().rbegin()->second.P0));
//                                _glidePlaneIntersections.emplace_back(boundingBoxSegments().begin()->second.P0,boundingBoxSegments().rbegin()->second.P0);
                                break;
                            }
                                
                            default:
                            {
                                model::cout<<"PlanarDislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
                                model::cout<<boundingBoxSegments();
//                                for(const auto& pair : boundingBoxSegments())
//                                {
//                                    model::cout<<"("<<pair.second.P0.transpose()<<","<<pair.second.P1.transpose()<<")"<<std::endl;
//                                }
                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
                            }
                        }
                    }
                    else
                    {
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                    // Now we must have exactly one _glidePlaneIntersections
//                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must have size 1");
                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");

                    break;
                }
                    
                default:
                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
//                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must exist");
                    assert(_glidePlaneIntersections && "_glidePlaneIntersections must exist");

                    // intersect the _glidePlaneIntersections with the new plane
//                    PlaneLineIntersection<dim> pli(lastGlidePlane.P,
//                                                   lastGlidePlane.unitNormal,
//                                                   _glidePlaneIntersections->P0, // origin of line
//                                                   _glidePlaneIntersections->P1-_glidePlaneIntersections->P0 // line direction
//                                                   );

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
                            for(const auto& plane : meshPlanes())
                            {
                                model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
                            }
                            model::cout<<"MeshPlane intersection is:"<<std::endl;
                            model::cout<<std::setprecision(15)<<std::scientific<<"  P0="<<_glidePlaneIntersections->P0.transpose()<<", P2="<<_glidePlaneIntersections->P1.transpose()<<std::endl;
                            
                            assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                            break;
                        }
                    }
                    
//                    if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
//                    {// nothing to do, _glidePlaneIntersections remains unchanged
//
//                    }
//                    else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
//                    {// _glidePlaneIntersections becomes a singular point
//                        _glidePlaneIntersections[0].P0 =pli.P;
//                        _glidePlaneIntersections[0].P1=pli.P;
//                    }
//                    else
//                    {
//                        model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
//                        model::cout<<"MeshPlanes are:"<<std::endl;
//                        for(const auto& plane : meshPlanes())
//                        {
//                            model::cout<<std::setprecision(15)<<std::scientific<<"  P="<<plane->P.transpose()<<", n="<<plane->unitNormal.transpose()<<std::endl;
//                        }
//                        model::cout<<"MeshPlane intersection is:"<<std::endl;
//                        model::cout<<std::setprecision(15)<<std::scientific<<"  P1="<<_glidePlaneIntersections[0].P0.transpose()<<", P2="<<_glidePlaneIntersections[0].P1.transpose()<<std::endl;
//
//                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                    }
                    
                }
                    
            }
            
            
            VerbosePlanarDislocationNode(2,"  _glidePlaneIntersections are: "<<_glidePlaneIntersections->P0.transpose()<<", "<<_glidePlaneIntersections->P1.transpose()<<std::endl;);
//            for(const auto& pair : _glidePlaneIntersections)
//            {
//                VerbosePlanarDislocationNode(2,"P1="<<std::setprecision(15)<<std::scientific<<pair.P0.transpose()<<", P2="<<pair.P1.transpose()<<std::endl;);
//
//            }
            
//            assert(_glidePlaneIntersections.size()<=1 && "_glidePlaneIntersections can have at the most size 1");
        }
        
        /**********************************************************************/
        bool addMeshPlane(const MeshPlaneType& gp)
        {
            const bool success=meshPlanes().insert(&gp).second;
            if(success)
            {
                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" addGlidePlane. meshPlanes().size()="<<meshPlanes().size()<<std::endl;);
                VerbosePlanarDislocationNode(4,"Current bounding box"<<std::endl;);
                VerbosePlanarDislocationNode(4,boundingBoxSegments()<<std::endl;);
                VerbosePlanarDislocationNode(4,"Plane bounding box"<<std::endl;);
                VerbosePlanarDislocationNode(4,BoundingLineSegments<dim>(gp)<<std::endl;);
                VerbosePlanarDislocationNode(4,"Plane meshIntersections: "<<gp.meshIntersections.size()<<std::endl;);
                //VerbosePlanarDislocationNode(4,BoundingLineSegments<dim>(gp)<<std::endl;);
                
                
                
                assert(gp.contains(this->get_P()) && "Glide Plane does not contain PlanarDislocationNode");
                boundingBoxSegments().updateWithMeshPlane(gp); // Update boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                assert((boundingBoxSegments().size() || !_isOnBoundingBox) && "EMPTY boundingBoxSegments");
                updateMeshPlaneIntersections(gp);
                
                VerbosePlanarDislocationNode(4,"new bounding box"<<std::endl;);
                VerbosePlanarDislocationNode(4,boundingBoxSegments()<<std::endl;);
                
                
                //                grains().insert(&this->network().poly.grain(gp.regionIDs.first));    // Insert new grain in grainSet
                //                grains().insert(&this->network().poly.grain(gp.regionIDs.second));   // Insert new grain in grainSet
            }
            return success;
        }
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes()
        {
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" adding GrainBoundaryPlanes"<<std::endl;);
            
            size_t addedGp=0;
            // Check if node is on a GB
            for(const auto& grain : grains())
            {
                for(const auto& gb : grain->grainBoundaries())
                {
                    VerbosePlanarDislocationNode(4,"GB "<<gb.second->tag()<<", d="<<gb.second->distanceTo(this->get_P())<<", contained="<<gb.second->contains(this->get_P())<<std::endl;);
                    if(gb.second->contains(this->get_P()))
                    {
                        grainBoundaries().insert(gb.second);
                        addedGp+=addMeshPlane(*gb.second);
                    }
                }
            }
            
            if(isGrainBoundaryNode())
            {
                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" added "<<addedGp<<" GrainBoundaryPlanes"<<std::endl;);
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
            
            VerbosePlanarDislocationNode(4,"snapping P="<<P.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            std::multimap<double,VectorDim> snapMap;
            
            // Collect possible snap points, sorted by distance to P
//            for(size_t k=0;k<boundingBoxSegments().size();++k)
            for(const auto& seg : boundingBoxSegments())
            {
                
                snapMap.emplace((P-seg.second.P0).squaredNorm(),seg.second.P0);
                snapMap.emplace((P-seg.second.P1).squaredNorm(),seg.second.P1);
                const VectorDim x(seg.second.snap(P));
                snapMap.emplace((P-x).squaredNorm(),x);

                
//                const auto& vertexPair(boundingBoxSegments()[k]);
//                const VectorDim segm(vertexPair.P1-vertexPair.P0);
//                const double segmNorm2(segm.squaredNorm());
//                if(segmNorm2>FLT_EPSILON)
//                {
//                    snapMap.emplace((P-vertexPair.P0).squaredNorm(),vertexPair.P0);
//                    snapMap.emplace((P-vertexPair.P1).squaredNorm(),vertexPair.P1);
//
//                    double u((P-vertexPair.P0).dot(segm)/segmNorm2);
//                    if(u>0.0 && u<1.0)
//                    {
//                        const VectorDim x(vertexPair.P0+u*segm);
//                        snapMap.emplace((P-x).squaredNorm(),x);
//                    }
//                }
//                else
//                {
//                    const VectorDim x(0.5*(vertexPair.P1+vertexPair.P0));
//                    snapMap.emplace((P-x).squaredNorm(),x);
//                }
            }
            
            VerbosePlanarDislocationNode(4,"there are "<<snapMap.size()<<" possible snap points."<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
            
            
            // Return the first point to which we can snap
            for(const auto& pair : snapMap)
            {
                VerbosePlanarDislocationNode(4,"Checking snap point "<<pair.second.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
                if(isMovableTo(pair.second))
                {
                    VerbosePlanarDislocationNode(4,"Snapping to "<<pair.second.transpose()<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
                    return pair.second;
                }
            }
            
            assert(false && "snapToBoundingBox FAILED.");
        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const VectorDim& X,const Simplex<dim,dim>* const guess) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (guess==NULL)
            {
                temp=this->network().mesh.search(X);
            }
            else
            {
                if(grains().size()==1)
                {// node only in one region
                    if((*grains().begin())->grainID!=guess->region->regionID)
                    {
                        temp=this->network().mesh.searchRegion((*grains().begin())->grainID,X);
                    }
                    else
                    {
                        temp=this->network().mesh.searchRegionWithGuess(X,guess);
                    }
                }
                else
                {
                    temp=this->network().mesh.searchWithGuess(X,guess);
                }
            }
            if(!temp.first) // PlanarDislocationNode not found inside mesh
            {// Detect if the PlanarDislocationNode is sligtly outside the boundary
                int faceID;
                const double baryMin(temp.second->pos2bary(X).minCoeff(&faceID));
                const bool isApproxOnBoundary(std::fabs(baryMin)<1.0e3*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
                if(!isApproxOnBoundary)
                {
                    model::cout<<"PlanarDislocationNode "<<this->sID<<" @ "<<X.transpose()<<std::endl;
                    model::cout<<"Simplex "<<temp.second->xID<<std::endl;
                    model::cout<<"bary "<<temp.second->pos2bary(X)<<std::endl;
                    model::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                    model::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                    model::cout<<"face of barymin is region-boundary Simplex? "<<temp.second->child(faceID).isRegionBoundarySimplex()<<std::endl;
                    assert(0 && "DISLOCATION NODE OUTSIDE MESH.");
                }
            }
            return temp.second;
        }
        
        /**********************************************************************/
//        BoundingLineSegments<dim> _glidePlaneIntersections; //
        std::unique_ptr<LineSegment<dim>> _glidePlaneIntersections; //

        
        bool _isGlissile;
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        //! The current velocity vector of *this PlanarDislocationNode
        VectorDofType velocity;
        //! The previous velocity vector of *this PlanarDislocationNode
        VectorDofType vOld;
        double velocityReductionCoeff;
        //! The normal unit vector of the boundary on which *this PlanarDislocationNode is moving on
        VectorDim boundaryNormal;
        bool _isOnBoundingBox;
        std::shared_ptr<NodeType> virtualNode;
        
//        BoundingLineSegments<TypeTraits<Derived>::dim> _boundingLineSegments;
        
    public:
        
        const NodeType* const masterNode;
        const bool isVirtualBoundaryNode;
        
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            use_velocityFilter=TextFileParser(fileName).readScalar<double>("use_velocityFilter",true);
            velocityReductionFactor=TextFileParser(fileName).readScalar<double>("velocityReductionFactor",true);
            assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
            //            EDR.readScalarInFile(fullName.str(),"verbosePlanarDislocationNode",NodeType::verbosePlanarDislocationNode);
            verbosePlanarDislocationNode=TextFileParser(fileName).readScalar<int>("verbosePlanarDislocationNode",true);
        }
        
        
        /**********************************************************************/
        PlanarDislocationNode(LoopNetworkType* const ln,
                              const VectorDim& Pin,
                              const VectorDofType& Vin,
                              const double& vrc) :
        /* base constructor */ NodeBaseType(ln,Pin)
        /* init */,_isGlissile(true)
        /* init */,p_Simplex(get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL))
        /* init */,velocity(Vin)
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(vrc)
        /* init */,boundaryNormal(SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol))
        /* init */,_isOnBoundingBox(boundaryNormal.squaredNorm()>FLT_EPSILON)
        /* init */,masterNode(nullptr)
        /* init */,isVirtualBoundaryNode(false)
        {/*! Constructor from DOF
          */
            VerbosePlanarDislocationNode(1,"Creating PlanarDislocationNode "<<this->sID<<" from position"<<std::endl;);
        }
        
        /**********************************************************************/
        PlanarDislocationNode(const LinkType& pL,
                              const double& u) :
        /* init */ NodeBaseType(pL.loopNetwork,pL.get_r(u))
        /* init */,_isGlissile(true)
        /* init */,p_Simplex(get_includingSimplex(this->get_P(),pL.source->includingSimplex()))
        /* init */,velocity((pL.source->velocity+pL.sink->velocity)*0.5) // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(std::min(pL.source->velocityReduction(),pL.sink->velocityReduction()))
        /* init */,boundaryNormal(SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) )
        /* init */,_isOnBoundingBox(pL.isBoundarySegment())
        /* init */,masterNode(nullptr)
        /* init */,isVirtualBoundaryNode(false)
        {/*! Constructor from ExpandingEdge and DOF
          */
            VerbosePlanarDislocationNode(1,"Creating PlanarDislocationNode "<<this->sID<<" from expanding "<<pL.source->sID<<"->"<<pL.sink->sID<<std::endl;);
            if(pL.isBoundarySegment())
            {
                assert(isBoundaryNode());
            }
            
        }
        
        /**********************************************************************/
        PlanarDislocationNode(LoopNetworkType* const ln,
                              const VectorDim& Pin,
                              const NodeType* const master) :
        /* base constructor */ NodeBaseType(ln,Pin)
        /* init */,_isGlissile(false)
        /* init */,p_Simplex(this->network().simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC? get_includingSimplex(this->get_P(),(const Simplex<dim,dim>*) NULL) : NULL)
        /* init */,velocity(VectorDim::Zero())
        /* init */,vOld(velocity)
        /* init */,velocityReductionCoeff(1.0)
        /* init */,boundaryNormal(VectorDim::Zero())
        /* init */,_isOnBoundingBox(false)
        /* init */,masterNode(master)
        /* init */,isVirtualBoundaryNode(true)
        {/*! Constructor from DOF
          */
            VerbosePlanarDislocationNode(1,"Creating VirtualPlanarDislocationNode "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        ~PlanarDislocationNode()
        {
            VerbosePlanarDislocationNode(1,"Destroying PlanarDislocationNode "<<this->sID<<" ("<<this<<")"<<std::endl;);
            //            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<", virtual node count="<<virtualNode.use_count()<<std::endl;);
            //            if(virtualNode)
            //            {
            //                this->network().remove(virtualNode->sID);
            //            }
        }
        
        /**********************************************************************/
        const std::shared_ptr<NodeType>& virtualBoundaryNode() const
        {
            return virtualNode;
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
//            return _boundingLineSegments;
            return *this;
        }
        
        /**********************************************************************/
        BoundingLineSegments<dim>& boundingBoxSegments()
        {
//            return _boundingLineSegments;
            return *this;
        }
        
//        /**********************************************************************/
//        const BoundingLineSegments<dim>& glidePlaneIntersections() const
//        {
//            return _glidePlaneIntersections;
//        }

        /**********************************************************************/
        const std::unique_ptr<LineSegment<dim>>& glidePlaneIntersections() const
        {
            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        VectorDim snapToMeshPlaneIntersection(const VectorDim& P)
        {
            
            if(_glidePlaneIntersections)
            {
                return _glidePlaneIntersections->snap(P);
            }
            else
            {
                VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 0"<<std::endl;);
                assert(meshPlanes().size()==1);
                return meshPlane(0).snapToPlane(P);
            }
            
//            switch (_glidePlaneIntersections.size())
//            {
//                case 0:
//                {
//                    //                    assert(glidePlanes().size()>0);
//                    VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 0"<<std::endl;);
//                    assert(meshPlanes().size()==1);
//                    return meshPlane(0).snapToPlane(P);
//                    break;
//                }
//
//                case 1:
//                {
//                    const VectorDim D(_glidePlaneIntersections[0].P1-_glidePlaneIntersections[0].P0);
//                    const double normD2(D.squaredNorm());
//                    if(normD2>FLT_EPSILON)
//                    {
//                        const double u=(P-_glidePlaneIntersections[0].P0).dot(D)/normD2;
//                        VerbosePlanarDislocationNode(3,"u="<<u<<std::endl;);
//
//                        if(u<0.0)
//                        {
//                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1a"<<std::endl;);
//                            return _glidePlaneIntersections[0].P0;
//                        }
//                        else if(u>1.0)
//                        {
//                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1b"<<std::endl;);
//                            return _glidePlaneIntersections[0].P1;
//                        }
//                        else
//                        {
//                            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1c"<<std::endl;);
//                            return _glidePlaneIntersections[0].P0+u*D;
//                        }
//                    }
//                    else
//                    {
//                        VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 1d"<<std::endl;);
//                        return _glidePlaneIntersections[0].P0;
//                    }
//
//                    break;
//                }
//
//                default:
//                {
//                    VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" snapToMeshPlaneIntersection, case 2"<<std::endl;);
//                    assert(0 && "THERE CAN BE AT MOST ONE LINE OF INTERSECTION");
//                    return VectorDim::Zero();
//                    break;
//                }
//            }
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" addLoopLink "<<pL->tag()<<std::endl;);
            
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
                
                if(boundingBoxSegments().contains(this->get_P()))
                {
                    _isOnBoundingBox=true;
                    setToBoundary(this->get_P());
                }
                
            }
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" finished addLoopLink "<<pL->tag()<<std::endl;);
            
            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" removeLoopLink "<<pL->tag()<<std::endl;);
            
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            if(!pL->loop()->isVirtualBoundaryLoop())
            {// Re-construct nodeConfinement
                
                _isGlissile=true;
                meshPlanes().clear();
                boundingBoxSegments().clear();
                _glidePlaneIntersections.reset(nullptr);
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
                
                if(boundingBoxSegments().contains(this->get_P()))
                {
                    _isOnBoundingBox=true;
                    setToBoundary(this->get_P());
                }
            }
            
            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" finished removeLoopLink "<<pL->tag()<<std::endl;);
            
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
        
        /**********************************************************************/
        const Simplex<dim,dim>* includingSimplex() const
        {/*!\returns A pointer to the const Simplex imcluding *this PlanarDislocationNode
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
                    model::cout<<"PlanarDislocationNode "<<this->sID<<", P="<<this->get_P().transpose()<<std::endl;
                    model::cout<<"_isOnBoundingBox="<<_isOnBoundingBox<<std::endl;
                    model::cout<<"isGrainBoundaryNode()="<<isGrainBoundaryNode()<<std::endl;
                    model::cout<<"NODE ON BoundingBox MUST BE EITHER A BOUNDARY NODE OR A GB NODE"<<std::endl;
                    assert(0 && "NODE ON BoundingBox MUST BE EITHER A BOUNDARY NODE OR A GB NODE");
                }
            }
            return isBndNode ;
        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return grainBoundaries().size();
        }
        
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
        bool isZeroBurgersNode() const
        {
            
            bool temp=true;
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=std::get<1>(neighborIter.second)->hasZeroBurgers();
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }

        /**********************************************************************/
        bool isSimpleZeroBurgersNode() const
        {
            bool temp(this->isSimple() && isZeroBurgersNode());
            if(temp)
            {// make sure attached sessile segments are aligned
                const LinkType* firstLink(std::get<1>(this->neighbors().begin()->second));
                const LinkType* secondLink(std::get<1>(this->neighbors().rbegin()->second));
                temp*=(firstLink->chord().normalized().cross(secondLink->chord().normalized()).norm()<FLT_EPSILON);
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
                      && (   isSimpleBoundaryNode()
                          || isSimpleGrainBoundaryNode()
                          || isSimpleSessileNode()
                          || isSimpleZeroBurgersNode()
                          || isGeometricallyRemovable(Lmin,cosRemove)
                          )
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
            if(!this->isBoundaryNode() && !this->isGrainBoundaryNode())
            {
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
        void resetVirtualBoundaryNode(const VectorDim& X)
        {
            if(isBoundaryNode() && !isVirtualBoundaryNode)
            {
                switch (this->network().simulationParameters.simulationType)
                {
                    case DefectiveCrystalParameters::FINITE_FEM:
                    {
                        if(virtualNode)
                        {
                            static_cast<NodeBaseType*>(virtualNode.get())->set_P(X+this->network().simulationParameters.virtualSegmentDistance*boundaryNormal);
                        }
                        else
                        {
                            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
                            virtualNode.reset(new NodeType(&this->network(),X+this->network().simulationParameters.virtualSegmentDistance*boundaryNormal,this->p_derived()));
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
                            VerbosePlanarDislocationNode(2,"PlanarDislocationNode "<<this->sID<<" resetting virtualBoundaryNode"<<std::endl;);
                            virtualNode.reset(new NodeType(&this->network(),X-((this->network().mesh.xMax()-this->network().mesh.xMin()).cwiseProduct(boundaryNormal)),this->p_derived()));
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
        }
        
        /**********************************************************************/
        void setToBoundary(const VectorDim& X)
        {
            _isOnBoundingBox=true;
            p_Simplex=get_includingSimplex(X,p_Simplex);
//            boundaryNormal=SimplexBndNormal::get_boundaryNormal(X,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
            boundaryNormal=boundingBoxSegments().boundaryNormal(X);
            assert(boundaryNormal.squaredNorm()>FLT_EPSILON); // automatically checks that bounding box contains X
            resetVirtualBoundaryNode(X);
            NodeBaseType::set_P(X); // in turn this calls PlanarDislocationSegment::updateGeometry, so the boundaryNormal must be computed before this line
//            assert(boundingBoxSegments().contains(this->get_P()));
        }
        
        /**********************************************************************/
        bool set_P(const VectorDim& newP)
        {
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" current P="<< this->get_P().transpose()<<"set_P to "<<newP.transpose()<<std::endl;);
            // make sure that node is on glide planes
            bool glidePlanesContained=true;
            for(const auto& gp : meshPlanes())
            {
                glidePlanesContained*=gp->contains(newP);
            }
            
            if(glidePlanesContained)
            {
                if(_isOnBoundingBox || boundingBoxSegments().contains(newP))
                {// node was on bounding box, it must remain on bounding box
                    const VectorDim X(snapToBoundingBox(newP));
                    setToBoundary(X);
                    if(addGrainBoundaryPlanes())
                    {// GB-planes were added, the bounding box has changed, so snap again
                        VerbosePlanarDislocationNode(3,"case 4"<<std::endl;);
                        const VectorDim X1(snapToBoundingBox(this->get_P()));
                        setToBoundary(X1);
                    }
                }
                else
                {// internal node
                    std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                    if(temp.first)
                    {// internal node, and newP is inside current grain
                        if(   isConnectedToBoundaryNodes()
                           && boundingBoxSegments().size()==2
                           && glidePlaneIntersections())
                        {// force special case to boundary to get rid of small debris
                            if((newP-glidePlaneIntersections()->P0).norm()<this->network().surfaceAttractionDistance)
                            {
                                setToBoundary(glidePlaneIntersections()->P0);
                            }
                            else if((newP-glidePlaneIntersections()->P1).norm()<this->network().surfaceAttractionDistance)
                            {
                                setToBoundary(glidePlaneIntersections()->P1);
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
                        setToBoundary(X);
                        if(addGrainBoundaryPlanes())
                        {// GB-planes were added, the bounding box has changed, and it may now be a set of degenerate lines
                            if(boundingBoxSegments().contains(this->get_P()))
                            {// new bounding box contains node
                                VerbosePlanarDislocationNode(3,"case 5"<<std::endl;);
                            }
                            else
                            {// new bounding box does not contain node
                                VerbosePlanarDislocationNode(3,"case 6"<<std::endl;);
                                const VectorDim X1(snapToMeshPlaneIntersection(this->get_P()));
                                _isOnBoundingBox=false;
                                boundaryNormal=boundingBoxSegments().boundaryNormal(X1);
                                //boundaryNormal=SimplexBndNormal::get_boundaryNormal(X1,*p_Simplex,bndTol); // must be updated before NodeBaseType::set_P
                                NodeBaseType::set_P(X1); // kill numerical errors
                            }
                        }
                    }
                }
                
                p_Simplex=get_includingSimplex(this->get_P(),p_Simplex); // update including simplex
                
                if(_isOnBoundingBox)
                {
//                    assert(boundingBoxSegments().contains(this->get_P()));
//                    boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // check if node is now on a boundary
                    boundaryNormal=boundingBoxSegments().boundaryNormal(this->get_P());
                    if(boundaryNormal.squaredNorm()<FLT_EPSILON && !isGrainBoundaryNode())
                    {
                        std::cout<<"PlanarDislocationNode "<<this->sID<<", @"<<std::setprecision(15)<<std::scientific<<this->get_P().transpose()<<std::endl;
                        std::cout<<"BoundingBox Lines:"<<std::endl;
                        std::cout<<boundingBoxSegments();
                        std::cout<<"BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL"<<std::endl;
                        //assert(false && "BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL");
                    }
                }
                else
                {
                    boundaryNormal.setZero();
                }
            }
            else
            {
                model::cout<<"PlanarDislocationNode "<<this->sID<<std::endl;
                assert(0 && "new position outside glide planes.");
            }
            
            
            
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" _isOnBoundingBox="<<_isOnBoundingBox<<std::endl;);
            const double posDelta((this->get_P()-newP).norm());
            VerbosePlanarDislocationNode(3,"PlanarDislocationNode "<<this->sID<<" posDelta="<<posDelta<<std::endl;);

            return posDelta<FLT_EPSILON;
        }
        
        /**********************************************************************/
        bool isMovableTo(const VectorDim& X) const
        {
            bool isMovable=true;
            
            VerbosePlanarDislocationNode(4,"checking if PlanarDislocationNode "<<this->sID<< " isMovable:"<<std::endl;);
            
            for(const auto& gp : meshPlanes())
            {// X must be contained by all glidePlanes
                isMovable*=gp->contains(X);
            }
            VerbosePlanarDislocationNode(4,"  meshPlanes contains X? "<<isMovable<<std::endl;);
            
            if(isMovable)
            {
                if(isOnBoundingBox())
                {
                    isMovable*=boundingBoxSegments().contains(X);
                }
            }
            
            if(isMovable)
            {
                
                
                
                for(const auto& pair : this->neighbors())
                {
                    if(std::get<1>(pair.second)->isBoundarySegment())
                    {// boundary segments other than must remain boundary if this node is moved
                        const bool bndNeighborMovable=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X));
                        VerbosePlanarDislocationNode(4,"  boundaryNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<bndNeighborMovable<<std::endl;);
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
                            VerbosePlanarDislocationNode(4,"  gbNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<gbNeighborMovable<<std::endl;);
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
                        const double currentNorm((std::get<0>(pair.second)->get_P()-this->get_P()).norm());
                        const double newNorm((std::get<0>(pair.second)->get_P()-X).norm());
                        if(currentNorm>FLT_EPSILON && newNorm>FLT_EPSILON)
                        {
                            const bool sessileNeighborMovable=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*currentNorm*newNorm);
                            VerbosePlanarDislocationNode(4,"  sessileNeighbor "<<std::get<1>(pair.second)->tag()<< " movable?"<<sessileNeighborMovable<<std::endl;);
                            isMovable*=sessileNeighborMovable;
                            if(!isMovable)
                            {
                                break;
                            }
                        }
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
        
        //        /**********************************************************************/
        //        template <class T>
        //        friend T& operator << (T& os, const NodeType& ds)
        //        {
        //            os<< DislocationNodeIO<dim>(ds);
        //            return os;
        //        }
        
    };
    
    // static data
    template <typename Derived,typename InterpolationType>
    bool PlanarDislocationNode<Derived,InterpolationType>::use_velocityFilter=true;
    
    template <typename Derived,typename InterpolationType>
    double PlanarDislocationNode<Derived,InterpolationType>::velocityReductionFactor=0.75;
    
    template <typename Derived,typename InterpolationType>
    const double PlanarDislocationNode<Derived,InterpolationType>::bndTol=1.0e-4;
    
    template <typename Derived,typename InterpolationType>
    int PlanarDislocationNode<Derived,InterpolationType>::verbosePlanarDislocationNode=0;
    
}
#endif
