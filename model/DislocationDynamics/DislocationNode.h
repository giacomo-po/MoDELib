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


namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationNode : public SplineNode<DislocationNode<_dim,corder,InterpolationType>,
    /*                                         */ _dim,corder,InterpolationType>
    
    {
        
    public:
        
        enum NodeMeshLocation{outsideMesh=-1, insideMesh=0, onMeshBoundary=1, onRegionBoundary=2};
        //enum BoundaryType {noBoundary=0, softBoundary=1, hardBoundary=2};
        
        
        constexpr static int dim=_dim; // make dim available outside class
        typedef DislocationNode       <dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment    <dim,corder,InterpolationType> LinkType;
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        typedef typename NodeBaseType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef GlidePlane<LoopType> GlidePlaneType;
        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        typedef std::deque<LatticePlane> SpecialLatticePlaneContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        //        typedef 	std::set<VectorDim,
        //        /*                */ CompareVectorsByComponent<double,dim,float>,
        //        /*                */ Eigen::aligned_allocator<VectorDim> > ConfiningPlaneIntersectionContainerType;
        typedef std::deque<std::pair<VectorDim,VectorDim>,Eigen::aligned_allocator<std::pair<VectorDim,VectorDim>>> BoundingSegmentsType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        //static constexpr double bndDistance=static_cast<double>(FLT_EPSILON);
        static const double bndDistance;
        
    private:
        
        bool _isGlissile;

        
        DislocationSharedObjects<dim> shared;
        
        std::set<const Grain<dim>*> grainSet; // this must be defined before p_Simplex

        
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        GlidePlaneContainerType _confiningPlanes;
        
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        
        //! The previous velocity vector of *this DislocationNode
        VectorDofType vOld;
        
        double velocityReductionCoeff;
        
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        bool _isOnBoundingBox;
        VectorDim boundaryNormal;
        
        VectorDim C;
        
        BoundingSegmentsType _boundingBoxSegments;
        BoundingSegmentsType _glidePlaneIntersections;

//        /**********************************************************************/
//        VectorDim boundingBoxProjection(const VectorDim& P) const;
        
        /**********************************************************************/
        VectorDim boundingBoxProjection(const VectorDim& P) const
        {
            
            std::map<double,VectorDim> snapMap;
            
            for(const auto& vertexPair : _boundingBoxSegments)
            {
                const VectorDim segm(vertexPair.second-vertexPair.first);
                const double segmNorm2(segm.squaredNorm());
                if(segmNorm2>FLT_EPSILON)
                {
                    double u((P-vertexPair.first).dot(segm)/segmNorm2);
                    if(u<0.0)
                    {
                        u=0.0;
                    }
                    if(u>1.0)
                    {
                        u=1.0;
                    }
                    const VectorDim x(vertexPair.first+u*segm);
                    snapMap.emplace((P-x).squaredNorm(),x);
                }
                else
                {
                    const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
                    snapMap.emplace((P-x).squaredNorm(),x);
                }
            }
            
            return snapMap.begin()->second;
            
        }
        

        
        /**********************************************************************/
        VectorDim snapToGlidePlaneIntersection(const VectorDim& P)
        {
            
            switch (_glidePlaneIntersections.size())
            {
                case 0:
                {
                    assert(_confiningPlanes.size()>0);
                    return glidePlane(0).snapToPlane(P);
                    break;
                }
                    
                case 1:
                {
                    return _glidePlaneIntersections[0].first+(P-_glidePlaneIntersections[0].first).dot(_glidePlaneIntersections[0].second)*_glidePlaneIntersections[0].second;
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
        void updateGlidePlaneIntersections(const GlidePlaneType& lastGlidePlane)
        {
            BoundingSegmentsType temp;
            
            switch (_confiningPlanes.size())
            {
                case 0:
                {
                    assert(0 && "AT LEAST ONE GLIDE PLANE MUST EXIST");
                    break;
                }
                    
                case 1:
                {
                    _glidePlaneIntersections.clear();
                    break;
                }
                    
                default:
                {
                    //const GlidePlaneType* const  lastGlidePlane(_confiningPlanes[_confiningPlanes.size()-1]);
                    if(_glidePlaneIntersections.size())
                    {/* some line of intersection already exist (maybe a degenerate point)
                      * intersect the last plane with the existing line and overwrite the line
                      */
                        assert(_glidePlaneIntersections.size()==1 && "THERE CANNOT BE MORE THAN ONE LINE OF INTERSECTION");
                        PlaneLineIntersection<dim> pli(lastGlidePlane.P.cartesian(),
                                                       lastGlidePlane.n.cartesian(),
                                                       _glidePlaneIntersections[0].first,
                                                       _glidePlaneIntersections[0].second);
                        
                        assert(pli.type==PlaneLineIntersection<dim>::COINCIDENT || pli.type==PlaneLineIntersection<dim>::INCIDENT);
                        _glidePlaneIntersections[0].first=pli.P;
                        _glidePlaneIntersections[0].second=pli.d;
                    }
                    else
                    {/* No line has been found for the previous planes.
                      * This means either _confiningPlanes.size()==1, or
                      * all previous planes are parallel (e.g. parallel planes in multiple grains). 
                      * Intersect last plane with first.
                      */
                        GlidePlaneObserver<LoopType>* const gpo(glidePlane(0).glidePlaneObserver);
                        for(const auto& otherGlidePlane : _confiningPlanes)
                        {
                            if(otherGlidePlane!=&lastGlidePlane)
                            {/* compute line of intersection between lastGlidePlane and any
                              * plane different from lastGlidePlane.
                              */
                                const PlanePlaneIntersection<dim>& ppi(gpo->glidePlaneIntersection(otherGlidePlane,&lastGlidePlane));
                                
                                if(ppi.type==PlanePlaneIntersection<dim>::INCIDENT)
                                {
                                    _glidePlaneIntersections.emplace_back(ppi.P,ppi.d);
                                }
                                break;
                            }
                        }
                        
                        
                    }
                    
                    break;
                }
                    
            }
            
        }
        
        
        /**********************************************************************/
        static BoundingSegmentsType updateBoundingSegments(const BoundingSegmentsType& old,
                                                           const GlidePlaneType& gp)
        {
            //model::cout<<"DislocationNode "<<this->sID<<" adding GlidePlane "<<gp.sID<<std::endl;
            
            
            BoundingSegmentsType temp;
            
            if(old.size())
            {
                for(const auto& oldPair : old)
                {
                    const BoundingSegmentsType psi=GlidePlaneObserver<LoopType>::planeSegmentIntersection(gp.P.cartesian(),
                                                                                                          gp.n.cartesian(),
                                                                                                          oldPair.first,
                                                                                                          oldPair.second);
                    if(psi.size())
                    {// plane and current segment intersect
                        for(size_t k=0;k<gp.meshIntersections.size();++k)
                        {
                            const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
                            
                            SegmentSegmentIntersection<dim> ssi(gp.meshIntersections[k].second,
                                                                gp.meshIntersections[k1].second,
                                                                oldPair.first,
                                                                oldPair.second);
                            
                            if(ssi.size)
                            {
                                temp.emplace_back(ssi.x0,ssi.x1);
                            }
                        }
                    }
                }
            }
            else
            {
                for(size_t k=0;k<gp.meshIntersections.size();++k)
                {
                    const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
                    temp.emplace_back(gp.meshIntersections[k].second,gp.meshIntersections[k1].second);
                }
            }
            
            //            assert(temp.size()==1 || temp.size()>=3 && "updateBoundingSegments failed");
            return temp;
            
        }
        
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

                    if(grainSet.size()==1)
                    {// node only in one region
                        //std::cout<<" 4 "<<std::flush;

                        if((*grainSet.begin())->grainID!=guess->region->regionID)
                        {
                            //std::cout<<" 5 "<<std::flush;

                            temp=DislocationSharedObjects<dim>::mesh.searchRegion((*grainSet.begin())->grainID,this->get_P());
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
            for(const auto& plane : _confiningPlanes)
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
        //                        const Simplex<dim,dim>* guess=(const Simplex<dim,dim>*) NULL) :
        /* base constructor */ NodeBaseType(Pin),
        //        /* init list        */ grain(shared.poly.grain(grainID)),
        //        /* init list        */ L(grain.latticeVector(Pin)),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex((const Simplex<dim,dim>*) NULL)),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(vrc),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero()),
        C(Pin)
        //        /* init list        */ grainBoundary_rID2(-1)
        //        oldP(this->get_P()),
        //        A1(this->get_P()),
        //                A2(this->get_P())
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        {/*! Constructor from DOF
          */
            std::cout<<"WARNING INITIALIZE C FROM INPUT FILE"<<std::endl;
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const VectorDim& Pin) :
        //        /* base constructor */ NodeBaseType(pL,Lin.cartesian()),
        /* base constructor */ NodeBaseType(Pin),
        //        /* init list        */ grain(pL.grain),
        //        /* init list        */ L(Lin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero()),
        //        otherGrains
        //        /* init list        */ grainBoundary_rID2((pL.source->grainBoundary_rID2==pL.sink->grainBoundary_rID2 && pL.sink->grainBoundary_rID2>0)? pL.sink->grainBoundary_rID2 : -1) // TO DO: CHANGE THIS FOR TRIPLE JUNCTIONS
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        //        oldP(this->get_P()),
        //        A1(this->get_P()),
        //        A2(this->get_P())
        C(this->get_P())
        {/*! Constructor from ExpandingEdge and DOF
          */
                        //std::cout<<"DislocationNode from ExpadingLink "<<this->sID<<std::endl;
            //forceBoundaryNode(pL);
            //            assert(0 && "Initialize C");
        }
        
        
        
        const GlidePlaneContainerType& confiningPlanes() const
        {
            return _confiningPlanes;
        }
        
        
        
        

        
        
        
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
                        //std::cout<<"DislocationNode "<<this->sID<<" addLoopLink"<<std::flush;
            NodeBaseType::addLoopLink(pL); // forward to base class
            _isGlissile*=pL->loop()->isGlissile;
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            const bool success=_confiningPlanes.insert(&(pL->loop()->glidePlane)).second;
            if(success)
            {
                assert(pL->loop()->glidePlane.contains(this->get_P()) && "Glide Plane does not contain DislocationNode");
                updateGlidePlaneIntersections(pL->loop()->glidePlane);
                _boundingBoxSegments=updateBoundingSegments(_boundingBoxSegments,pL->loop()->glidePlane); // Update _boundingBoxSegments
                grainSet.insert(&(pL->loop()->grain)); // Insert new grain in grainSet
                if(grainSet.size()>1)
                {
                    std::cout<<"WARNING: CHECK THAT NODE IS ON REGION BND"<<std::endl;
                }
                
                const VectorDim bbP(boundingBoxProjection(this->get_P()));
                if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                {
                    set_P(bbP);
                    _isOnBoundingBox=true;
                }

                
                //                boxCenter.setZero();
                //                for(const auto& posPair : _boundingBoxSegments)
                //                {
                //                    boxCenter+=posPair.first;
                //                    boxCenter+=posPair.second;
                //                }
                //                boxCenter/=(2*_boundingBoxSegments.size());
                

            }
            
            //std::cout<<" done"<<std::endl;

            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
                        //std::cout<<"DislocationNode "<<this->sID<<" removeLoopLink"<<std::flush;
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            // Re-construct _confiningPlanes and grainSet
            _isGlissile=true;
            _confiningPlanes.clear();
            _boundingBoxSegments.clear();
            grainSet.clear();
            for(const auto& loopLink : this->loopLinks())
            {
                //                std::cout<<"loopLink "<<loopLink->source()->sID<<"->"<<loopLink->sink()->sID<<std::endl;
                _isGlissile*=loopLink->loop()->isGlissile;
                const bool success=_confiningPlanes.insert(&(loopLink->loop()->glidePlane)).second;
                
                if(success)
                {
                    assert(loopLink->loop()->glidePlane.contains(this->get_P()) && "Glide Plane does not contain DislocationNode");
                    updateGlidePlaneIntersections(loopLink->loop()->glidePlane);
                    _boundingBoxSegments=updateBoundingSegments(_boundingBoxSegments,loopLink->loop()->glidePlane);
                    grainSet.insert(&(loopLink->loop()->grain));

                    
                }
            }
            
            if(grainSet.size()>1)
            {
                std::cout<<"WARNING: CHECK THAT NODE IS ON REGION BND"<<std::endl;
            }
            
            const VectorDim bbP(boundingBoxProjection(this->get_P()));
            if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
            {
                set_P(bbP);
                _isOnBoundingBox=true;
            }
            
            //std::cout<<" done"<<std::endl;
            //            boxCenter.setZero();
            //            for(const auto& posPair : _boundingBoxSegments)
            //            {
            //                boxCenter+=posPair.first;
            //                boxCenter+=posPair.second;
            //            }
            //            boxCenter/=(2*_boundingBoxSegments.size());
        }
        
        
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const
        {
            VectorOfNormalsType temp;
            
            if(_isGlissile)
            {
                for(const auto& plane : _confiningPlanes)
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
        
        bool isOscillating() const
        {
            return velocityReductionCoeff<std::pow(velocityReductionFactor,3);
        }
        
        /**********************************************************************/
        const GlidePlaneType& glidePlane(const size_t& n) const
        {
            assert(n<_confiningPlanes.size());
            auto iter=_confiningPlanes.begin();
            std::advance(iter,n);
            return **iter;
        }
        
        /**********************************************************************/
        void set_P(const VectorDim& P_in)
        {
            
            // make sure that node is on glide planes
            for(const auto& gp : _confiningPlanes)
            {
                assert(gp->contains(P_in) && "NodePosition outside GlidePlane");
            }
            
            NodeBaseType::set_P(P_in);
            if(shared.use_boundary) // using confining mesh
            {
                p_Simplex=get_includingSimplex(p_Simplex);
                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance); // check if node is now on a boundary
            }
            
            make_projectionMatrix();
            
        }
        
        /**********************************************************************/
        void snapToBoundingBox(const VectorDim& P)
        {
            
            set_P(boundingBoxProjection(P));
            _isOnBoundingBox=true;
            
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
                const VectorDim newP=snapToGlidePlaneIntersection(this->get_P()+dX);
                
                if(shared.use_boundary)
                {// using confining mesh
                    
                    if(_isOnBoundingBox || grainSet.size()>1)
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
            //            os  << ds.sID<<"\t"
            //            /**/<< std::setprecision(15)<<std::scientific<<ds.get_P().transpose()<<"\t"
            //            /**/<< ds.get_V().transpose()<<"\t"
            //            /**/<< ds.velocityReduction()<<"\t"
            //            /**/<< ds.pSN()->sID<<"\t"
            //            /**/<< (ds.meshLocation()==onMeshBoundary);
            return os;
        }
        
    };
    
    
    // static data
    template <int _dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<_dim,corder,InterpolationType>::use_velocityFilter=true;
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    double DislocationNode<_dim,corder,InterpolationType>::velocityReductionFactor=0.75;
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    const double DislocationNode<_dim,corder,InterpolationType>::bndDistance=FLT_EPSILON;
    
    
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
    
}
#endif





