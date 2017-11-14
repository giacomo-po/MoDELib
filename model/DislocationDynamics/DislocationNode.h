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
    /*          */ private std::set<const GrainBoundary<typename TypeTraits<DislocationNode<_dim,corder,InterpolationType>>::LoopNetworkType>*>,
    /*          */ private std::set<const Grain<typename TypeTraits<DislocationNode<_dim,corder,InterpolationType>>::LoopNetworkType>*>,
    /*          */ private std::set<const GlidePlane<typename TypeTraits<DislocationNode<_dim,corder,InterpolationType>>::LoopNetworkType>*>,
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
        typedef GlidePlane<LoopNetworkType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef std::set<const Grain<LoopNetworkType>*> GrainContainerType;
        typedef std::set<const GrainBoundary<LoopNetworkType>*> GrainBoundaryContainerType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        static const double bndTol;
        static int verboseDislocationNode;
        
    private:
        
        /**********************************************************************/
        void updateGlidePlaneIntersections(const GlidePlaneType& lastGlidePlane)
        {
            BoundingLineSegments<dim> temp;
            
            switch (glidePlanes().size())
            {
                case 0:
                {// there must be at least one glide plane
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
                    _glidePlaneIntersections.clear();
                    break;
                }
                    
                case 2:
                {// a second plane is being added, so we must have no _glidePlaneIntersections
                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
                    
                    // Grab the infinite line of intersection between the two planes
                    GlidePlaneObserver<LoopNetworkType>* const gpo(glidePlane(0).glidePlaneObserver);
                    const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&glidePlane(0),&glidePlane(1)));
                    
                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
                    {/* Two distinct glide planes can be coincident only if they belong to different grains
                      * In that case, the intersection of their bounding boxes should be one line segment
                      */
                        assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
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
                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
                                break;
                            }
                                
                            case 2:
                            {// The two intersections must be degenerate (2 boundary points)
//                                std::cout<<boundingBoxSegments()[0].first.transpose()<<std::endl;
//                                std::cout<<boundingBoxSegments()[0].second.transpose()<<std::endl;
//                                std::cout<<boundingBoxSegments()[1].first.transpose()<<std::endl;
//                                std::cout<<boundingBoxSegments()[1].second.transpose()<<std::endl;
                                
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
                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
                    }
                    
                }
                    
            }
            assert(_glidePlaneIntersections.size()<=1 && "_glidePlaneIntersections can have at the most size 1");
        }
        
        
        
//        /**********************************************************************/
//        void clear()
//        {
//            glidePlanes().clear();
//            boundingBoxSegments().clear();
//            _glidePlaneIntersections.clear();
//            grains().clear();
//        }
        
        /**********************************************************************/
        bool addGlidePlane(const GlidePlaneType& gp)
        {
            const bool success=glidePlanes().insert(&gp).second;
            if(success)
            {
                VerboseDislocationNode(1,"DislocationNode "<<this->sID<<" addGlidePlane. glidePlanes().size()="<<glidePlanes().size()<<std::endl;);
                assert(gp.contains(this->get_P()) && "Glide Plane does not contain DislocationNode");
                //                _isGlissile*=pL->loop()->isGlissile;
                boundingBoxSegments().updateWithGlidePlane(gp); // Update _boundingBoxSegments. This must be called before updateGlidePlaneIntersections
                updateGlidePlaneIntersections(gp);
                //                grains().insert(&(gp.grain)); // Insert new grain in grainSet
                //                if(grains().size()>1)
                //                {
                //                    std::cout<<"WARNING: CHECK THAT NODE IS ON REGION BND"<<std::endl;
                //                }
                //
                //                const VectorDim bbP(_boundingBoxSegments.snap(this->get_P()));
                //                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                //                {
                //                    set_P(bbP);
                //                    _isOnBoundingBox=true;
                //                }
            }
            return success;
        }
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes() __attribute__ ((deprecated)) // HERE glidePlanes().begin() IS TEMPORARY, UNTIL WE STORE THE GLIDE PLANE OF THE CSL AND DSCL
        {
            size_t addedGp=0;
            // Check if node is on a GB
            for(const auto& gb : this->network().poly.grainBoundaries())
            {
                
                const auto& gp(gb.second.glidePlanes().begin()->second); // HERE glidePlanes().begin() IS TEMPORARY, UNTIL WE STORE THE GLIDE PLANE OF THE CSL AND DSCL
                if(gp.contains(this->get_P()))
                {
                    grainBoundaries().insert(&gb.second);
                    addedGp+=addGlidePlane(gp);
                }
                
//                for(const auto& gp : gb.second.glidePlanes())
//                {
//                    if(gp.second.contains(this->get_P()))
//                    {
//                        grainBoundaries().insert(&gb.second);
//                        addedGp+=addGlidePlane(gp.second);
//                    }
//                }
            }
            
            if(addedGp)
            {
                VerboseDislocationNode(1,"DislocationNode "<<this->sID<<" adding "<<addedGp<<" GrainBoundaryPlanes"<<std::endl;);
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
            
            
            const VectorDim pL=std::get<0>(boundingBoxSegments().snap(P));
            const VectorDim pV=boundingBoxSegments().snapToVertex(P).second;
            
            bool pLcontained=true;
            bool pVcontained=true;
            
            for(const auto& pair : this->neighbors())
            {
                VerboseDislocationNode(2,"checking "<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<std::endl;);
                
                if(std::get<1>(pair.second)->isBoundarySegment())
                {// boundary segments must not become internal
                    VerboseDislocationNode(2,std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<" is boundary"<<std::endl;);
                    
                    pLcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pL+std::get<0>(pair.second)->get_P())).first;
                    pVcontained*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pV+std::get<0>(pair.second)->get_P())).first;
                    VerboseDislocationNode(2,"pLcontained="<<pLcontained<<std::endl;);//<<", lineID="<<pLcontained.second<<std::endl;
                    VerboseDislocationNode(2,"pVcontained="<<pVcontained<<std::endl;);//", lineID="<<pVcontained.second<<std::endl;
                    
                }
                
                if(std::get<1>(pair.second)->isGrainBoundarySegment())
                {// grainBoundary segments must not become internal
                    for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
                    {
                        pLcontained*=gb->glidePlanes().begin()->second.contains(pL);
                        pVcontained*=gb->glidePlanes().begin()->second.contains(pV);
                    }
                }
            }
            
            if(pLcontained)
            {
                return pL;
            }
            else
            {
                if(pVcontained)
                {
                    return pV;
                }
                else
                {
                    model::cout<<"DislocationNode "<<this->sID<<" snapToBoundingBox FAILED."<<std::endl;
                    assert(false && "snapToBoundingBox FAILED.");
                    return VectorDim::Zero();
                }
            }
            
        }
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
            //std::cout<<"DislocationNode "<<this->sID<<" get_includingSimplex "<<std::flush;
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (this->network().use_boundary)
            {
                //std::cout<<" 1 "<<std::flush;
                
                if (guess==NULL)
                {
                    //std::cout<<" 2 "<<std::flush;
                    
                    temp=this->network().mesh.search(this->get_P());
                }
                else
                {
                    //std::cout<<" 3 "<<std::flush;
                    
                    if(grains().size()==1)
                    {// node only in one region
                        //std::cout<<" 4 "<<std::flush;
                        
                        if((*grains().begin())->grainID!=guess->region->regionID)
                        {
                            //std::cout<<" 5 "<<std::flush;
                            
                            temp=this->network().mesh.searchRegion((*grains().begin())->grainID,this->get_P());
                        }
                        else
                        {
                            //std::cout<<" 6 "<<std::flush;
                            
                            temp=this->network().mesh.searchRegionWithGuess(this->get_P(),guess);
                        }
                    }
                    else
                    {
                        //std::cout<<" 7 "<<std::flush;
                        
                        //                        std::cout<<"WARNING: CHECK THAT NODE IS ON THE REGION BOUNDARY"<<std::endl;
                        temp=this->network().mesh.searchWithGuess(this->get_P(),guess);
                    }
                }
                //std::cout<<" 8 "<<std::flush;
                
                
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
            }
            
            //std::cout<<" done"<<std::endl;
            
            
            return temp.second;
        }
        
        
        /**********************************************************************/
        void make_projectionMatrix()
        {
            
            Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
            VectorOfNormalsType  CN;
            for(const auto& plane : glidePlanes())
            {
                CN.push_back(plane->n.cartesian().normalized());
            }
            
            if(meshLocation()==MeshLocation::onMeshBoundary)
            {
                CN.push_back(boundaryNormal);
                
            }
            
            // Add normal to region boundary
            //            CN.push_back(regionBndNormal);
            
            // Find independent vectors
            GramSchmidt::orthoNormalize(CN);
            
            // Assemble projection matrix (prjM)
            this->prjM.setIdentity();
            for (size_t k=0;k<CN.size();++k)
            {
                this->prjM*=( I-CN[k]*CN[k].transpose() );
            }
            
        }
        
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
        VectorDim C;
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
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
        /* init list        */ boundaryNormal(this->network().use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) : VectorDim::Zero()),
        /* init list        */ C(Pin)
        {/*! Constructor from DOF
          */
            
            VerboseDislocationNode(1,"Creating DislocationNode "<<this->sID<<" from position"<<std::endl;);
            
            addGrainBoundaryPlanes();
            //            set_P(this->get_P()); // trigger _isOnBoundingBox and _isGrainBoundaryNode
            std::cout<<"WARNING INITIALIZE _isOnBoundingBox"<<std::endl;
            //            std::cout<<"WARNING INITIALIZE _isGrainBoundaryNode"<<std::endl;
            std::cout<<"WARNING INITIALIZE C FROM INPUT FILE"<<std::endl;
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const VectorDim& Pin) :
        /* base constructor */ NodeBaseType(pL.loopNetwork,Pin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ _isOnBoundingBox(pL.isBoundarySegment()),
        /* init list        */ boundaryNormal(this->network().use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) : VectorDim::Zero()),
        /* init list        */ C(this->get_P())
        {/*! Constructor from ExpandingEdge and DOF
          */
            
            VerboseDislocationNode(1,"Creating DislocationNode "<<this->sID<<" from expansion"<<std::endl;);
            
            
            addGrainBoundaryPlanes();
            //            std::cout<<"WARNING INITIALIZE _isGrainBoundaryNode"<<std::endl;
            std::cout<<"WARNING INITIALIZE C FROM INPUT FILE"<<std::endl;
            //            set_P(this->get_P()); // trigger _isOnBoundingBox and _isGrainBoundaryNode
        }
        
        /**********************************************************************/
        ~DislocationNode()
        {
            VerboseDislocationNode(1,"Destroying DislocationNode "<<this->sID<<std::endl;);
        }
        
        /**********************************************************************/
        const GlidePlaneContainerType& glidePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        GlidePlaneContainerType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GlidePlaneType& glidePlane(const size_t& n) const
        {
            assert(n<glidePlanes().size());
            auto iter=glidePlanes().begin();
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
        VectorDim snapToGlidePlaneIntersection(const VectorDim& P)
        {
            
            switch (_glidePlaneIntersections.size())
            {
                case 0:
                {
                    //                    assert(glidePlanes().size()>0);
                    assert(glidePlanes().size()==1);
                    return glidePlane(0).snapToPlane(P);
                    break;
                }
                    
                case 1:
                {
                    const VectorDim D=_glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first;
                    const double normD2(D.squaredNorm());
                    if(normD2>FLT_EPSILON)
                    {
                        const double u=(P-_glidePlaneIntersections[0].first).dot(D)/normD2;
                        if(u<0.0)
                        {
                            return _glidePlaneIntersections[0].first;
                        }
                        else if(u>1.0)
                        {
                            return _glidePlaneIntersections[0].second;
                        }
                        else
                        {
                            return _glidePlaneIntersections[0].first+u*D;
                        }
                    }
                    else
                    {
                        return _glidePlaneIntersections[0].first;
                    }
                    
                    break;
                }
                    
                default:
                {
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
            
            VerboseDislocationNode(1,"DislocationNode "<<this->sID<<" addLoopLink"<<std::endl;);
            
            NodeBaseType::addLoopLink(pL); // forward to base class
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            const bool success = addGlidePlane(pL->loop()->glidePlane);
            if(success)
            {
                _isGlissile*=pL->loop()->isGlissile;
                
                if(boundingBoxSegments().contains(this->get_P()).first)
                {
                    _isOnBoundingBox=true;
                }
                
                
                //                if(boundingBoxSegments().size())
                //                {
                //                    const VectorDim bbP(std::get<0>(boundingBoxSegments().snap(this->get_P())));
                //                    if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                //                    {
                //                        set_P(bbP);
                //                        _isOnBoundingBox=true;
                //                    }
                //                }
                //                else
                //                {
                //                    _isGlissile=false;
                //                }
            }
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            
            VerboseDislocationNode(1,"DislocationNode "<<this->sID<<" removeLoopLink"<<std::endl;);

            
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            
            
            // Re-construct nodeConfinement
            _isGlissile=true;
            glidePlanes().clear();
            boundingBoxSegments().clear();
            _glidePlaneIntersections.clear();
            grains().clear();
            
            for(const auto& loopLink : this->loopLinks())
            {
                const bool success = addGlidePlane(loopLink->loop()->glidePlane);
                
                if(success)
                {
                    _isGlissile*=loopLink->loop()->isGlissile;
                }
            }
            
            addGrainBoundaryPlanes();
            
            if(boundingBoxSegments().contains(this->get_P()).first)
            {
                //                    const VectorDim bbP(std::get<0>(boundingBoxSegments().snap(this->get_P())));
                //                    if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                //                    {
                //                        set_P(bbP);
                _isOnBoundingBox=true;
                //                    }
            }
            
            //            if(boundingBoxSegments().size())
            //            {// not the last segment being removed
            //                const VectorDim bbP(std::get<0>(boundingBoxSegments().snap(this->get_P())));
            //                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
            //                {
            //                    set_P(bbP);
            //                    _isOnBoundingBox=true;
            //                }
            //            }
            //            else
            //            {
            //                _isGlissile=false;
            //            }
            
            //            std::cout<<"DislocationNode "<<this->sID<<" removeLoopLink. _isOnBoundingBox="<<_isOnBoundingBox<<std::endl;
            
        }
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const
        {
            VectorOfNormalsType temp;
            
            if(_isGlissile)
            {
                for(const auto& plane : glidePlanes())
                {
                    temp.push_back(plane->n.cartesian().normalized());
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
        
        
        /**********************************************************************/
        void set_V(const VectorDofType& vNew)
        {
            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
        }
        
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
        
        /**********************************************************************/
        const VectorDofType& get_V() const
        {/*! The nodal velocity vector
          */
            return velocity;
        }
        
        /**********************************************************************/
        void applyVelocityFilter(double vMaxGood)
        {
            if(use_velocityFilter)
            {
                const double Rc=3.0*10.0;
                const double filterThreshold=0.05*vMaxGood*vMaxGood;
                
                if((this->get_P()-C).norm()<Rc)
                {
                    if(velocity.dot(vOld)<-filterThreshold)
                    {
                        velocityReductionCoeff*=velocityReductionFactor;
                    }
                }
                else
                {
                    if(velocity.dot(vOld)<-filterThreshold)
                    {
                        velocityReductionCoeff*=velocityReductionFactor;
                        C=this->get_P();
                    }
                    else if(velocity.dot(vOld)>0.0)
                    {
                        velocityReductionCoeff/=velocityReductionFactor;
                    }
                    else
                    {
                        // don't change velocityReductionCoeff
                    }
                }
                
                //                const double filterThreshold=0.05*velocity.norm()*vOld.norm();
                
                //                const double VdotVold=
                //                if(velocity.dot(vOld)<-filterThreshold)
                //                {
                //                    velocityReductionCoeff*=velocityReductionFactor;
                //                }
                //                else if(velocity.dot(vOld)>filterThreshold)
                //                {
                //                    velocityReductionCoeff/=velocityReductionFactor;
                //                }
                //                else
                //                {
                //                    // don't change velocityReductionCoeff
                //                }
                if(velocityReductionCoeff>1.0)
                {
                    velocityReductionCoeff=1.0;
                }
                if(velocityReductionCoeff<0.0001)
                {
                    velocityReductionCoeff=0.0001;
                }
                
                
                
                if( isOscillating()
                   //&& velocity.norm()>vMaxGood /*velocity.norm()>0.0*/
                   )
                { // oscillating node
                    std::cout<<"node "<<this->sID<<" BAD, ";
                    //                    velocity*=(velocityReductionCoeff*vMaxGood/velocity.norm());
                    if(velocity.norm()>vMaxGood)
                    {
                        velocity=velocityReductionCoeff*vMaxGood*velocity.normalized();
                    }
                    //                    std::cout<<velocity.norm()<<std::endl;
                }
                else
                { // good node
                    std::cout<<"node "<<this->sID<<" GOOD, ";
                    velocity*=velocityReductionCoeff;
                }
                
                //                velocity*=velocityReductionCoeff;
                std::cout<<"velocityReductionCoeff="<<velocityReductionCoeff<<", vMaxGood="<<vMaxGood<<", velocity.norm()="<<velocity.norm()<<", newVelocity="<<velocity.norm()<<std::endl;
                
            }
        }
        
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
        
        
        
        /**********************************************************************/
        bool isOscillating() const
        {
            return velocityReductionCoeff<std::pow(velocityReductionFactor,3);
        }
        
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
            
            
            if(_isOnBoundingBox)
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
            
            //            if(boundaryNormal.squaredNorm())
            //            {
            //                temp=onMeshBoundary;
            //            }
            //            else
            //            {
            //                if(_isOnBoundingBox)
            //                {
            //                    temp=onRegionBoundary;
            //
            //                }
            //                else
            //                {
            //                    temp=insideMesh;
            //                }
            //            }
            
            return temp;
        }
        
        /**********************************************************************/
        const bool& isOnBoundingBox() const
        {
            return _isOnBoundingBox;
        }
        
        /**********************************************************************/
        const bool& isBoundaryNode() const
        {
            return _isOnBoundingBox;
        }
        
        //        /**********************************************************************/
        //        const bool& isGrainBoundaryNode() const
        //        {
        //            return _isGrainBoundaryNode;
        //            //            return meshLocation()==onRegionBoundary;
        //        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return grainBoundaries().size();
            //            return meshLocation()==onRegionBoundary;
        }
        
        /**********************************************************************/
        bool isPureBoundaryNode() const
        {
            return isBoundaryNode() && isConnectedToBoundaryNodes();
        }
        //
        /**********************************************************************/
        bool isConnectedToBoundaryNodes() const
        {
            bool temp(true);
            for (const auto& neighborIter : this->neighbors())
            {
                temp*=std::get<0>(neighborIter.second)->isBoundaryNode();
            }
            
            return temp;
        }
        
        
        /**********************************************************************/
        bool isSimpleBndNode() const
        {
            bool temp=false;
            if(isOnBoundingBox())
            {
                if(this->isSimple())
                {
                    temp=true;
                    
                    //                    std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> normalDeq;
                    std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> chordDeq;
                    
                    for (const auto& neighborIter : this->neighbors())
                    {
                        //                        if (!std::get<2>(neighborIter.second)==0)
                        //                        {
                        if (
                            //std::get<1>(neighborIter.second)->isBoundarySegment()
                            !std::get<1>(neighborIter.second)->hasZeroBurgers()
                            )
                        {  // neighbor not searched
                            //temp*=std::get<0>(neighborIter.second)->isOnBoundingBox();
                            temp*=std::get<1>(neighborIter.second)->isBoundarySegment();
                            //                            normalDeq.push_back(std::get<0>(neighborIter.second)->bndNormal());
                            chordDeq.push_back(std::get<1>(neighborIter.second)->chord());
                        }
                        //                        }
                    }
                    
                    
                    //                    if(normalDeq.size())
                    //                    {
                    //                        for(const auto& bndN : normalDeq)
                    //                        {
                    //                            temp*=((bndN-normalDeq[0]).squaredNorm()<FLT_EPSILON);
                    //                        }
                    //                    }
                    //                    else
                    //                    {
                    //                        temp=false;
                    //                    }
                    
                    if(chordDeq.size())
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
                
            }
            
            return temp;
        }
        
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
        
        /******************************************************************************/
        void neighborsAt(const LatticeVectorType& L0, std::set<size_t>& temp) const
        {/*!\param[in] P0 position to be serached
          * \param[out]temp set of IDs of neighbors of this which are located at P0 (possibly including *this)
          * \param[in] tol tolerance used to detect position overlap
          */
            for (const auto& nIiter : this->neighbors())
            { // loop over neighborhood
                if((std::get<0>(nIiter.second)->get_L()-L0).squaredNorm()==0)
                { // a neighbor of I exists at P0
                    temp.insert(std::get<0>(nIiter.second)->sID);
                }
            }
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
        void set_P(const VectorDim& newP)
        {
            VerboseDislocationNode(1,"DislocationNode "<<this->sID<<" set_P"<<std::endl;);
            // make sure that node is on glide planes
            bool glidePlanesContained=true;
            for(const auto& gp : glidePlanes())
            {
                glidePlanesContained*=gp->contains(newP);
            }
            
            
            if(glidePlanesContained)
            {
                //                NodeBaseType::set_P(newP);
                if(this->network().use_boundary) // using confining mesh
                {
                    
                    std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                    if(temp.first)
                    {// newP is inside current grain
                        if(_isOnBoundingBox)
                        {// if the node is already on the the bounding box, keep it there
                            VerboseDislocationNode(2,"case 1"<<std::endl;);
                            NodeBaseType::set_P(snapToBoundingBox(newP));
                        }
                        else
                        {// node was internal to the grain and remains internal
                            VerboseDislocationNode(2,"case 2"<<std::endl;);
                            NodeBaseType::set_P(newP);
                        }
                    }
                    else
                    {// newP is outside current grain (or on grain boundary)
                        VerboseDislocationNode(2,"case 3"<<std::endl;);
                        NodeBaseType::set_P(snapToBoundingBox(newP)); // put node on the bouding box
                        if(_isOnBoundingBox)
                        {// node was on bounding box, and exited the grain
                            if(addGrainBoundaryPlanes())
                            {// GB-planes were added, the bounding box has changed, so snap again
                                VerboseDislocationNode(2,"case 4"<<std::endl;);
                                NodeBaseType::set_P(snapToBoundingBox(newP)); // put node back on bounding box
                            }
                        }
                        else
                        {// node was internal and exited the grain
                            if(addGrainBoundaryPlanes())
                            {// GB-planes were added, the bounding box has changed
                                if(boundingBoxSegments().contains(this->get_P()).first)
                                {// new bounding box contains node
                                    VerboseDislocationNode(2,"case 5"<<std::endl;);
                                    NodeBaseType::set_P(snapToBoundingBox(newP)); // kill numerical errors
                                    _isOnBoundingBox=true;
                                }
                                else
                                {// new bounding box does not contain node
                                    VerboseDislocationNode(2,"case 6"<<std::endl;);
                                    NodeBaseType::set_P(snapToGlidePlaneIntersection(newP)); // kill numerical errors
                                    _isOnBoundingBox=false;
                                }
                            }
                            else
                            {// node is now on bounding box, and no GB planes were added
                                _isOnBoundingBox=true;
                            }
                        }
                        
                    }
                    
                    p_Simplex=get_includingSimplex(p_Simplex); // update including simplex
                    
                    if(_isOnBoundingBox)
                    {
                        assert(boundingBoxSegments().contains(this->get_P()).first);
                        
                        boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // check if node is now on a boundary
                        if(boundaryNormal.squaredNorm()<FLT_EPSILON)
                        {
                            model::cout<<"DislocationNode "<<this->sID<<", @"<<this->get_P().transpose()<<std::endl;
                            assert(false && "BOUNDARY NODES MUST HAVE A NON-ZERO NORMAL");
                        }
                        
                        
                    }
                    
                }
                else
                {
                    VerboseDislocationNode(2,"case 7"<<std::endl;);
                    NodeBaseType::set_P(newP);
                }
                
                
                make_projectionMatrix();
            }
            else
            {
                model::cout<<"DislocationNode "<<this->sID<<std::endl;
                assert(0 && "new position outside glide planes.");
            }
            
        }
        
        /**********************************************************************/
        bool isMovableTo(const VectorDim& X) const
        {
            bool isMovable=true;
            
            
            for(const auto& gp : glidePlanes())
            {// X must be contained by all glidePlanes
                isMovable*=gp->contains(X);
            }
            
            //            std::cout<<"A "<<isMovable<<std::endl;
            
            if(isMovable)
            {
                
                //                HERE CHECK THAT X IS IN REGION
                
                for(const auto& pair : this->neighbors())
                {
                    //                            std::cout<<"neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
                    
                    if(std::get<1>(pair.second)->isSessile())
                    {// sessile segments cannot change direction
                        
                        isMovable*=((std::get<0>(pair.second)->get_P()-X).cross(std::get<0>(pair.second)->get_P()-this->get_P()).norm()<FLT_EPSILON*(std::get<0>(pair.second)->get_P()-X).norm()*(std::get<0>(pair.second)->get_P()-this->get_P()).norm());
                        
                        //                        isMovable*=LineSegment<dim>(std::get<0>(pair.second)->get_P(),X).contains(this->get_P());
                        //                        isMovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
                    }
                    
                    //                                std::cout<<"B "<<isMovable<<std::endl;
                    
                    if(std::get<1>(pair.second)->isBoundarySegment())
                    {// boundary segments other than A->B must remain boundary
                        isMovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
                    }
                    
                    //                                std::cout<<"C "<<isMovable<<std::endl;
                    
                    if(std::get<1>(pair.second)->isGrainBoundarySegment())
                    {// boundary segments other than A->B must remain boundary
                        for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
                        {
                            //                                        std::cout<<"grainBoundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nB->get_P())<<std::endl;
                            isMovable*=gb->glidePlanes().begin()->second.contains(X);
                        }
                    }
                    
                    //                                std::cout<<"D "<<isMovable<<std::endl;
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
                const VectorDim newP=snapToGlidePlaneIntersection(this->get_P()+dX);
                
                set_P(newP);
                
                //                if(this->network().use_boundary)
                //                {// using confining mesh
                //
                //                    if(_isOnBoundingBox
                //                       //|| grains().size()>1
                //                       )
                //                    {// if the node is already on the the bounding box, keep it there
                //                        snapToBoundingBox(newP);
                //                    }
                //                    else
                //                    {// check if node is ou
                //
                //
                //                        //                        std::set<const Simplex<dim,dim>*> path;
                //                        //                        const bool searchAllRegions=false;
                //                        std::pair<bool,const Simplex<dim,dim>*> temp(this->network().mesh.searchRegionWithGuess(newP,p_Simplex));
                //                        THIS SEARCH SHOULD BE DONE IN set_P. IF REGIONS BEFORE AND AFTER DON'T MATHC, TRIGGER addGrainBoundaryPlanes
                //
                //                        if(temp.first)
                //                        {// newP is inside box
                //                            set_P(newP);
                //                        }
                //                        else
                //                        {
                //                            snapToBoundingBox(newP);
                //                            addGrainBoundaryPlanes(); // bounding box changes
                //
                //                        }
                //
                //
                //                    }
                //
                //
                //                }
                //                else
                //                {// no confining mesh, move freely
                //                    set_P(newP);
                //                }
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
            
            vOld=velocity; // store current value of velocity before updating
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
