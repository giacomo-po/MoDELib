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
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/Math/GramSchmidt.h>
#include <model/Mesh/Simplex.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/SegmentSegmentIntersection.h>
#include <model/Geometry/PlanePlaneIntersection.h>
#include <model/Geometry/PlaneLineIntersection.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/DislocationNodeConfinement.h>

namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationNode :
    /*          */ public SplineNode<DislocationNode<_dim,corder,InterpolationType>,_dim,corder,InterpolationType>,
    /*          */ public DislocationNodeConfinement<DislocationNode<_dim,corder,InterpolationType>>
    
    {
        
    public:
        
        enum NodeMeshLocation{outsideMesh=-1, insideMesh=0, onMeshBoundary=1, onRegionBoundary=2};
        //enum BoundaryType {noBoundary=0, softBoundary=1, hardBoundary=2};
        
        
        constexpr static int dim=_dim; // make dim available outside class
        typedef DislocationNode       <dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment    <dim,corder,InterpolationType> LinkType;
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        typedef DislocationNodeConfinement<NodeType>  DislocationNodeConfinementType;
        
        typedef typename NodeBaseType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        //        typedef std::deque<LatticePlane> SpecialLatticePlaneContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        //        typedef 	std::set<VectorDim,
        //        /*                */ CompareVectorsByComponent<double,dim,float>,
        //        /*                */ Eigen::aligned_allocator<VectorDim> > ConfiningPlaneIntersectionContainerType;
        //        typedef std::deque<std::pair<VectorDim,VectorDim>,Eigen::aligned_allocator<std::pair<VectorDim,VectorDim>>> LineSegmentContainerType;
        //        typedef typename BoundingLineSegments<dim>:: LineSegmentContainerType;
        
        //        typedef typename BoundingLineSegments<dim>::LineSegmentContainerType LineSegmentContainerType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        //static constexpr double bndTol=static_cast<double>(FLT_EPSILON);
        static const double bndTol;
        
    private:
        
        bool _isGlissile;
        
        
        DislocationSharedObjects<dim> shared;
        
//        std::set<const Grain<dim>*> grainSet; // this must be defined before p_Simplex
        
        
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        //        GlidePlaneContainerType _confiningPlanes;
        
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        
        //! The previous velocity vector of *this DislocationNode
        VectorDofType vOld;
        
        double velocityReductionCoeff;
        
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        bool _isOnBoundingBox;
        VectorDim boundaryNormal;
        
        VectorDim C;
        
        //LineSegmentContainerType _boundingBoxSegments; // this is the intersection of the bounding boxes of each glide plane
//        BoundingLineSegments<dim> _boundingBoxSegments;
        
//        BoundingLineSegments<dim> _glidePlaneIntersections; //
        
        /**********************************************************************/
        void snapToBoundingBox(const VectorDim& P)
        {
            
            set_P(boundingBoxSegments().snap(P));
            _isOnBoundingBox=true;
            
        }
        
//        /**********************************************************************/
//        VectorDim snapToGlidePlaneIntersection(const VectorDim& P)
//        {
//            
//            switch (_glidePlaneIntersections.size())
//            {
//                case 0:
//                {
////                    assert(nodeConfinement().glidePlanes().size()>0);
//                    assert(nodeConfinement().glidePlanes().size()==1);
//                    return glidePlane(0).snapToPlane(P);
//                    break;
//                }
//                    
//                case 1:
//                {
//                    const VectorDim D=_glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first;
//                    const double normD2(D.squaredNorm());
//                    return normD2>FLT_EPSILON? _glidePlaneIntersections[0].first+(P-_glidePlaneIntersections[0].first).dot(D)*D/normD2 : _glidePlaneIntersections[0].first;
//                    break;
//                }
//                    
//                default:
//                {
//                    assert(0 && "THERE CAN BE AT MOST ONE LINE OF INTERSECTION");
//                    return VectorDim::Zero();
//                    break;
//                }
//            }
//            
//        }
        
//        /**********************************************************************/
//        void updateGlidePlaneIntersections(const GlidePlaneType& lastGlidePlane)
//        {
//            BoundingLineSegments<dim> temp;
//            
//            switch (nodeConfinement().glidePlanes().size())
//            {
//                case 0:
//                {// there must be at least one glide plane
//                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
//                    break;
//                }
//                    
//                case 1:
//                {// if there is only one glide plane, then _glidePlaneIntersections must be empty
//                    _glidePlaneIntersections.clear();
//                    break;
//                }
//                    
//                case 2:
//                {// a second plane is being added, so we must have no _glidePlaneIntersections
//                    assert(_glidePlaneIntersections.size()==0 && "_glidePlaneIntersections must be empty");
//                    
//                    // Grab the infinite line of intersection between the two planes
//                    GlidePlaneObserver<LoopType>* const gpo(glidePlane(0).glidePlaneObserver);
//                    const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(&glidePlane(0),&glidePlane(1)));
//                    
//                    if(ppi.type==PlanePlaneIntersection<dim>::COINCIDENT)
//                    {/* Two distinct glide planes can be coincident only if they belong to different grains
//                      * In that case, the intersection of their bounding boxes should be one line segment
//                      */
//                        assert(boundingBoxSegments().size()==1 && "There should be only one line in boundingBoxSegments()");
//                        _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
//                    }
//                    else if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
//                    {/* If the two planes are incident then the intersection of
//                      * their bounding boxes is either a pair of singluar segments (2 points)
//                      * or a line segment on the boundary
//                      */
//                        switch (boundingBoxSegments().size())
//                        {
//                            case 1:
//                            {// the bounding boxes of the two planes intersect on a boundary segment. Add end points to _glidePlaneIntersections
//                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[0].second);
//                                break;
//                            }
//                                
//                            case 2:
//                            {// The two intersections must be degenerate (2 boundary points)
//                                assert((boundingBoxSegments()[0].first-boundingBoxSegments()[0].second).squaredNorm()<FLT_EPSILON);
//                                assert((boundingBoxSegments()[1].first-boundingBoxSegments()[1].second).squaredNorm()<FLT_EPSILON);
//                                _glidePlaneIntersections.emplace_back(boundingBoxSegments()[0].first,boundingBoxSegments()[1].first);
//                                break;
//                            }
//                                
//                            default:
//                            {
//                                model::cout<<"DislocationNode "<<this->sID<<" boundingBoxSegments() are:"<<std::endl;
//                                for(const auto& pair : boundingBoxSegments())
//                                {
//                                    model::cout<<"("<<pair.first.transpose()<<","<<pair.second.transpose()<<")"<<std::endl;
//                                }
//                                assert(0 && "Bounding boxes of two incident planes must intersect on a boundary segment or on two boundary points.");
//                            }
//                        }
//                    }
//                    else
//                    {
//                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                    }
//                    
//                    // Now we must have exactly one _glidePlaneIntersections
//                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must have size 1");
//                    
//                    break;
//                }
//                    
//                default:
//                {// Case of more that 2 planes. A _glidePlaneIntersections must exist
//                    assert(_glidePlaneIntersections.size()==1 && "_glidePlaneIntersections must exist");
//                    
//                    // intersect the _glidePlaneIntersections with the new plane
//                    PlaneLineIntersection<dim> pli(lastGlidePlane.P.cartesian(),
//                                                   lastGlidePlane.n.cartesian(),
//                                                   _glidePlaneIntersections[0].first, // origin of line
//                                                   _glidePlaneIntersections[0].second-_glidePlaneIntersections[0].first // line direction
//                                                   );
//                    
//                    if(pli.type==PlaneLineIntersection<dim>::COINCIDENT)
//                    {// nothing to do, _glidePlaneIntersections remains unchanged
//                        
//                    }
//                    else if(pli.type==PlaneLineIntersection<dim>::INCIDENT)
//                    {// _glidePlaneIntersections becomes a singular point
//                        _glidePlaneIntersections[0].first =pli.P;
//                        _glidePlaneIntersections[0].second=pli.P;
//                    }
//                    else
//                    {
//                        assert(0 && "Intersection must be COINCIDENT or INCIDENT.");
//                    }
//                    
//                }
//                    
//            }
//            
//        }
        
        
        
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
            //std::cout<<"DislocationNode "<<this->sID<<" get_includingSimplex "<<std::flush;
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (DislocationSharedObjects<dim>::use_boundary)
            {
                //std::cout<<" 1 "<<std::flush;
                
                if (guess==NULL)
                {
                    //std::cout<<" 2 "<<std::flush;
                    
                    temp=DislocationSharedObjects<dim>::mesh.search(this->get_P());
                }
                else
                {
                    //std::cout<<" 3 "<<std::flush;
                    
                    if(nodeConfinement().grains().size()==1)
                    {// node only in one region
                        //std::cout<<" 4 "<<std::flush;
                        
                        if((*nodeConfinement().grains().begin())->grainID!=guess->region->regionID)
                        {
                            //std::cout<<" 5 "<<std::flush;
                            
                            temp=DislocationSharedObjects<dim>::mesh.searchRegion((*nodeConfinement().grains().begin())->grainID,this->get_P());
                        }
                        else
                        {
                            //std::cout<<" 6 "<<std::flush;
                            
                            temp=DislocationSharedObjects<dim>::mesh.searchRegionWithGuess(this->get_P(),guess);
                        }
                    }
                    else
                    {
                        //std::cout<<" 7 "<<std::flush;
                        
                        std::cout<<"WARNING: CHECK THAT NODE IS ON THE REGION BOUNDARY"<<std::endl;
                        temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(this->get_P(),guess);
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
            for(const auto& plane : nodeConfinement().glidePlanes())
            {
                CN.push_back(plane->n.cartesian().normalized());
            }
            
            if(meshLocation()==onMeshBoundary)
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
        
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**********************************************************************/
        DislocationNode(const VectorDim& Pin,
                        //                        const int& grainID,
                        const VectorDofType& Vin,
                        const double& vrc) :
        /* base constructor */ NodeBaseType(Pin),
        /* base constructor */ DislocationNodeConfinementType(this),
        //        /* init list        */ grain(shared.poly.grain(grainID)),
        //        /* init list        */ L(grain.latticeVector(Pin)),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex((const Simplex<dim,dim>*) NULL)),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(vrc),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) : VectorDim::Zero()),
        /* init list        */ C(Pin)
        {/*! Constructor from DOF
          */
            std::cout<<"WARNING INITIALIZE C FROM INPUT FILE"<<std::endl;
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const VectorDim& Pin) :
        /* base constructor */ NodeBaseType(Pin),
        /* base constructor */ DislocationNodeConfinementType(this),
        //        /* init list        */ grain(pL.grain),
        //        /* init list        */ L(Lin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol) : VectorDim::Zero()),
        /* init list        */ C(this->get_P())
        {/*! Constructor from ExpandingEdge and DOF
          */
        }
        
        /**********************************************************************/
        const DislocationNodeConfinementType& nodeConfinement() const
        {
            return *this;
        }
        
        /**********************************************************************/
        DislocationNodeConfinementType& nodeConfinement()
        {
            return *this;
        }
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            NodeBaseType::addLoopLink(pL); // forward to base class
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            const bool success = nodeConfinement().addGlidePlane(pL->loop()->glidePlane);
            if(success)
            {
                _isGlissile*=pL->loop()->isGlissile;
//                updateGlidePlaneIntersections(pL->loop()->glidePlane);
                
                const VectorDim bbP(boundingBoxSegments().snap(this->get_P()));
                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                {
                    set_P(bbP);
                    _isOnBoundingBox=true;
                }
            }
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            // Re-construct nodeConfinement
            _isGlissile=true;
            nodeConfinement().clear();
//            _glidePlaneIntersections.clear();

            for(const auto& loopLink : this->loopLinks())
            {
                const bool success = nodeConfinement().addGlidePlane(loopLink->loop()->glidePlane);
                
                if(success)
                {
                    _isGlissile*=loopLink->loop()->isGlissile;
//                    updateGlidePlaneIntersections(loopLink->loop()->glidePlane);
                }
            }
            
            if(boundingBoxSegments().size())
            {// not the last segment being removed
                const VectorDim bbP(boundingBoxSegments().snap(this->get_P()));
                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                {
                    set_P(bbP);
                    _isOnBoundingBox=true;
                }
            }
        }
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const
        {
            VectorOfNormalsType temp;
            
            if(_isGlissile)
            {
                for(const auto& plane : nodeConfinement().glidePlanes())
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
        
        /**********************************************************************/
        bool is_simple() const
        {
            size_t nonZeroLink=0;
            for (const auto& neighborIter : this->neighbors())
            {
                if (!std::get<2>(neighborIter.second)==0)
                {
                    if (std::get<1>(neighborIter.second)->burgers().squaredNorm())
                    {  // neighbor not searched
                        nonZeroLink++;
                    }
                }
            }
            return (nonZeroLink==2);
        }
        
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
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& glidePlaneIntersections() const
        {
            return nodeConfinement().glidePlaneIntersections();

//            return _glidePlaneIntersections;
        }
        
        /**********************************************************************/
        const BoundingLineSegments<dim>& boundingBoxSegments() const
        {
            return nodeConfinement().boundingBoxSegments();
//            return _boundingBoxSegments;
        }
        
        //        /**********************************************************************/
        //        const GlidePlaneContainerType& confiningPlanes() const
        //        {
        //            return _confiningPlanes;
        //        }
        
        /**********************************************************************/
        bool isOscillating() const
        {
            return velocityReductionCoeff<std::pow(velocityReductionFactor,3);
        }
        
//        /**********************************************************************/
//        const GlidePlaneType& glidePlane(const size_t& n) const
//        {
//            assert(n<nodeConfinement().glidePlanes().size());
//            auto iter=nodeConfinement().glidePlanes().begin();
//            std::advance(iter,n);
//            return **iter;
//        }
        
        /**********************************************************************/
        void set_P(const VectorDim& P_in)
        {
            
            // make sure that node is on glide planes
            for(const auto& gp : nodeConfinement().glidePlanes())
            {
                assert(gp->contains(P_in) && "NodePosition outside GlidePlane");
            }
            
            NodeBaseType::set_P(P_in);
            if(shared.use_boundary) // using confining mesh
            {
                p_Simplex=get_includingSimplex(p_Simplex);
                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndTol); // check if node is now on a boundary
            }
            
            make_projectionMatrix();
            
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
        NodeMeshLocation meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            
            NodeMeshLocation temp = outsideMesh;
            
            if(boundaryNormal.squaredNorm())
            {
                temp=onMeshBoundary;
            }
            else
            {
                if(_isOnBoundingBox)
                {
                    temp=onRegionBoundary;
                    
                }
                else
                {
                    temp=insideMesh;
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        bool isBoundaryNode() const
        {
            return meshLocation()==onMeshBoundary;
        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return meshLocation()==onRegionBoundary;
        }
        
        /**********************************************************************/
        bool isPureBoundaryNode() const
        {
            return isBoundaryNode() && isConnectedToBoundaryNodes();
        }
        
        /**********************************************************************/
        bool isConnectedToBoundaryNodes() const
        {
            bool temp(true);
            for (const auto& neighborIter : this->neighbors())
            {
                if (std::get<2>(neighborIter.second)) // not self
                {
                    temp*=std::get<0>(neighborIter.second)->isBoundaryNode();
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        bool isPureGBNode() const
        {
            bool temp(!this->is_isolated());
            for (const auto& neighborIter : this->neighbors())
            {
                if (std::get<2>(neighborIter.second)) // not self
                {
                    temp*=std::get<0>(neighborIter.second)->isGrainBoundaryNode();
                }
            }
            
            return (isGrainBoundaryNode()&&temp);
        }
        
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
        const bool& isOnBoundingBox() const
        {
            return _isOnBoundingBox;
        }
        
        /**********************************************************************/
        void move(const double & dt,const double& dxMax)
        {
            
            //velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
            
            
            
            
            
            //const VectorDim P_old(this->get_P());
            
            //VectorDim dX=2.0*A2-A1-this->get_P();
            
            VectorDim dX=velocity.template segment<dim>(0)*dt;
            //VectorDim dX=velocity.template segment<dim>(0)*dt*velocityReductionCoeff;
            //VectorDim dX=(1.5*velocity.template segment<dim>(0)-0.5*vOld)*dt; // Adams–Bashforth
            //VectorDim dX=(0.5*velocity.template segment<dim>(0)+0.5*vOld)*dt; // Adams–Bashforth
            
            
            //Limit dX for boundaryNodes bec
            
            const double dXnorm(dX.norm());
            if((isBoundaryNode() || isConnectedToBoundaryNodes()) && dXnorm>dxMax)
            {
                dX*=dxMax/dXnorm;
            }
            
            if (dX.squaredNorm()>0.0 && _isGlissile) // move a node only if |v|!=0
            {
                
                // Make sure that new position is at intersection of glidePlanes
                const VectorDim newP=nodeConfinement().snapToGlidePlaneIntersection(this->get_P()+dX);
                
                if(shared.use_boundary)
                {// using confining mesh
                    
                    if(_isOnBoundingBox || nodeConfinement().grains().size()>1)
                    {// if the node is already on the the bounding box, keep it there
                        snapToBoundingBox(newP);
                    }
                    else
                    {// check if node is ou
                        
                        
                        //                        std::set<const Simplex<dim,dim>*> path;
                        //                        const bool searchAllRegions=false;
                        std::pair<bool,const Simplex<dim,dim>*> temp(DislocationSharedObjects<dim>::mesh.searchRegionWithGuess(newP,p_Simplex));
                        if(temp.first)
                        {// newP is inside box
                            set_P(newP);
                        }
                        else
                        {
                            snapToBoundingBox(newP);
                        }
                        
                        
                    }
                    
                    
                }
                else
                {// no confining mesh, move freely
                    set_P(newP);
                }
            }
            else
            {
                velocity.setZero();
            }
            
            //            // Store actual velocity
            //            if(dt>0.0)
            //            {
            //                velocity=(this->get_P()-P_old)/dt;
            //            }
            
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
    const double DislocationNode<_dim,corder,InterpolationType>::bndTol=FLT_EPSILON;
    
    
}
#endif


//        /**********************************************************************/
//        VectorDim boundingBoxProjection(const VectorDim& P) const
//        {
//
//            std::map<double,VectorDim,std::less<double>,Eigen::aligned_allocator<std::pair<double,VectorDim>>> snapMap;
//
//            for(const auto& vertexPair : _boundingBoxSegments)
//            {
//                const VectorDim segm(vertexPair.second-vertexPair.first);
//                const double segmNorm2(segm.squaredNorm());
//                if(segmNorm2>FLT_EPSILON)
//                {
//                    double u((P-vertexPair.first).dot(segm)/segmNorm2);
//                    if(u<0.0)
//                    {
//                        u=0.0;
//                    }
//                    if(u>1.0)
//                    {
//                        u=1.0;
//                    }
//                    const VectorDim x(vertexPair.first+u*segm);
//                    snapMap.emplace((P-x).squaredNorm(),x);
//                }
//                else
//                {
//                    const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
//                    snapMap.emplace((P-x).squaredNorm(),x);
//                }
//            }
//
//            return snapMap.begin()->second;
//
//        }

//        /**********************************************************************/
//        static LineSegmentContainerType updateBoundingSegments(const LineSegmentContainerType& old,
//                                                           const GlidePlaneType& gp)
//        {
//            //model::cout<<"DislocationNode "<<this->sID<<" adding GlidePlane "<<gp.sID<<std::endl;
//
//
//            LineSegmentContainerType temp;
//
//            if(old.size())
//            {
//                for(const auto& oldPair : old)
//                {
//                    const LineSegmentContainerType psi=GlidePlaneObserver<LoopType>::planeSegmentIntersection(gp.P.cartesian(),
//                                                                                                          gp.n.cartesian(),
//                                                                                                          oldPair.first,
//                                                                                                          oldPair.second);
//                    if(psi.size())
//                    {// plane and current segment intersect
//                        for(size_t k=0;k<gp.meshIntersections.size();++k)
//                        {
//                            const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
//
//                            SegmentSegmentIntersection<dim> ssi(gp.meshIntersections[k].second,
//                                                                gp.meshIntersections[k1].second,
//                                                                oldPair.first,
//                                                                oldPair.second);
//
//                            if(ssi.size)
//                            {
//                                temp.emplace_back(ssi.x0,ssi.x1);
//                            }
//                        }
//                    }
//                }
//            }
//            else
//            {
//                for(size_t k=0;k<gp.meshIntersections.size();++k)
//                {
//                    const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
//                    temp.emplace_back(gp.meshIntersections[k].second,gp.meshIntersections[k1].second);
//                }
//            }
//
//            //            assert(temp.size()==1 || temp.size()>=3 && "updateBoundingSegments failed");
//            return temp;
//
//        }

//    /**********************************************************************/
//    template <int _dim, short unsigned int corder, typename InterpolationType,
//    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
//    VectorDim DislocationNode<_dim,corder,Interpolation,QuadratureRule>::boundingBoxProjection(const VectorDim& P) const
//    {
//
//        std::map<double,VectorDim> snapMap;
//
//        for(const auto& vertexPair : _boundingBoxSegments)
//        {
//            const VectorDim segm(vertexPair.second-vertexPair.first);
//            const double segmNorm2(segm.squaredNorm());
//            if(segmNorm2>FLT_EPSILON)
//            {
//                double u((P-vertexPair.first).dot(segm)/segmNorm2);
//                if(u<0.0)
//                {
//                    u=0.0;
//                }
//                if(u>1.0)
//                {
//                    u=1.0;
//                }
//                const VectorDim x(vertexPair.first+u*segm);
//                snapMap.emplace((P-x).squaredNorm(),x);
//            }
//            else
//            {
//                const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
//                snapMap.emplace((P-x).squaredNorm(),x);
//            }
//        }
//
//        return snapMap.begin()->second;
//
//    }



