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
//#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/Math/GramSchmidt.h>
#include <model/Mesh/Simplex.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/PlanePlaneIntersection.h>
#include <model/Geometry/PlaneLineIntersection.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/DislocationNodeConfinement.h>

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
    /*          */ public DislocationNodeConfinement<DislocationNode<_dim,corder,InterpolationType>>
    
    {
        
    public:
        
//        enum NodeMeshLocation{outsideMesh=-1, insideMesh=0, onMeshBoundary=1, onRegionBoundary=2};
        //enum BoundaryType {noBoundary=0, softBoundary=1, hardBoundary=2};
        
        
        constexpr static int dim=_dim; // make dim available outside class
        typedef DislocationNode       <dim,corder,InterpolationType> NodeType;
        typedef DislocationSegment    <dim,corder,InterpolationType> LinkType;
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        typedef DislocationNodeConfinement<NodeType>  DislocationNodeConfinementType;
        
        typedef typename NodeBaseType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<NodeType>::LoopType LoopType;
        typedef typename TypeTraits<NodeType>::LoopNetworkType LoopNetworkType;
        
        constexpr static int NdofXnode=NodeBaseType::NdofXnode;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef GlidePlane<LoopNetworkType> GlidePlaneType;
        //        typedef std::set<const GlidePlaneType*> GlidePlaneContainerType;
        //        typedef std::deque<LatticePlane> SpecialLatticePlaneContainerType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        
        typedef typename TypeTraits<NodeType>::MeshLocation MeshLocation;
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
        static int verboseDislocationNode;
    private:
        
        
        
        /**********************************************************************/
        size_t addGrainBoundaryPlanes()
        {
            size_t addedGp=0;
            // Check if node is on a GB
            for(const auto& gb : this->network().poly.grainBoundaries())
            {
                const GlidePlaneType& gp(gb.second.glidePlanes().begin()->second);// HERE BEGIN IS TEMPORARY, UNTIL WE STORE THE GLIDE PLANE OF THE CSL AND DSCL
                
                if(gp.contains(this->get_P()))
                {
                    addedGp+=nodeConfinement().addGlidePlane(gp);
                    //                    std::cout<<"DIslocationNode "<<this->sID<<" addGrainBoundaryPlanes"<<std::endl;
                    //                    if(nodeConfinement().addGlidePlane(gp))
                    //                    {// adding the gb-plane was successful
                    //                        _isGrainBoundaryNode=true;
                    //
                    //                        // The bounding box had changed, so see if _isOnBoundingBox applies
                    //                        if(nodeConfinement().contains(this->get_P()))
                    //                        {
                    //                            _isOnBoundingBox=true;
                    //                        }
                    //                        else
                    //                        {
                    //                            _isOnBoundingBox=false;
                    //                        }
                    //                    }
                }
            }
            
            if(addedGp)
            {
                VerboseDislocationNode(2,"DIslocationNode "<<this->sID<<" addGrainBoundaryPlanes"<<std::endl;);
                _isGrainBoundaryNode=true;
                
                
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
            
            
            const VectorDim pL=std::get<0>(nodeConfinement().boundingBoxSegments().snap(P));
            const VectorDim pV=nodeConfinement().boundingBoxSegments().snapToVertex(P).second;
            
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
                    //                    for(const auto& pair : this->neighbors())
                    //                    {
                    //                        std::cout<<"checking "<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<std::endl;
                    //                        if(std::get<1>(pair.second)->isBoundarySegment() ||
                    //                           std::get<1>(pair.second)->isGrainBoundarySegment() )
                    //                        {
                    //                            std::cout<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<" is boundary"<<std::endl;
                    //                            const auto pLtemp=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pL+std::get<0>(pair.second)->get_P()));
                    //                            const auto pVtemp=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pV+std::get<0>(pair.second)->get_P()));
                    //                            std::cout<<"pLcontained="<<pLtemp.first<<", lineID="<<pLtemp.second<<std::endl;
                    //                            std::cout<<"pVcontained="<<pVtemp.first<<", lineID="<<pVtemp.second<<std::endl;
                    //                        }
                    //                    }
                    //
                    //                    for(const auto& pair : this->neighbors())
                    //                    {
                    //                        std::cout<<"checking "<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<std::endl;
                    //                        if(std::get<1>(pair.second)->isBoundarySegment())
                    //                        {// boundary and grainBoundary segments must not become internal
                    //                            const auto pLtemp=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pL+std::get<0>(pair.second)->get_P()));
                    //                            const auto pVtemp=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(pV+std::get<0>(pair.second)->get_P()));
                    //                            std::cout<<"pLcontained="<<pLtemp.first<<", lineID="<<pLtemp.second<<std::endl;
                    //                            std::cout<<"pVcontained="<<pVtemp.first<<", lineID="<<pVtemp.second<<std::endl;
                    //                        }
                    //
                    //                        if(std::get<1>(pair.second)->isGrainBoundarySegment())
                    //                        {// boundary and grainBoundary segments must not become internal
                    //                            std::cout<<std::get<1>(pair.second)->source->sID<<"->"<<std::get<1>(pair.second)->sink->sID<<" is grainBoundary"<<std::endl;
                    //
                    //                            bool pLtemp=true;
                    //                            bool pVtemp=true;
                    //
                    //
                    //                            for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
                    //                            {
                    //                                pLtemp*=gb->glidePlanes().begin()->second.contains(pL);
                    //                                pVtemp*=gb->glidePlanes().begin()->second.contains(pV);
                    //                            }
                    //                            std::cout<<"pLcontained="<<pLtemp<<std::endl;
                    //                            std::cout<<"pVcontained="<<pVtemp<<std::endl;
                    //                        }
                    //                    }
                    //
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
                    
                    if(nodeConfinement().grains().size()==1)
                    {// node only in one region
                        //std::cout<<" 4 "<<std::flush;
                        
                        if((*nodeConfinement().grains().begin())->grainID!=guess->region->regionID)
                        {
                            //std::cout<<" 5 "<<std::flush;
                            
                            temp=this->network().mesh.searchRegion((*nodeConfinement().grains().begin())->grainID,this->get_P());
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
                        
                        std::cout<<"WARNING: CHECK THAT NODE IS ON THE REGION BOUNDARY"<<std::endl;
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
            for(const auto& plane : nodeConfinement().glidePlanes())
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
        bool _isGrainBoundaryNode;
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
        /* base constructor */ DislocationNodeConfinementType(this),
        //        /* init list        */ grain(this->network().poly.grain(grainID)),
        //        /* init list        */ L(grain.latticeVector(Pin)),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex((const Simplex<dim,dim>*) NULL)),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(vrc),
        /* init list        */ _isOnBoundingBox(false),
        /* init list        */ _isGrainBoundaryNode(false),
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
        /* base constructor */ DislocationNodeConfinementType(this),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ _isOnBoundingBox(pL.isBoundarySegment()),
        /* init list        */ _isGrainBoundaryNode(false),
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
            
            
            //            std::cout<<"DislocationNode "<<this->sID<<" addLoopLink"<<std::endl;
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            const bool success = nodeConfinement().addGlidePlane(pL->loop()->glidePlane);
            if(success)
            {
                _isGlissile*=pL->loop()->isGlissile;
                
                if(nodeConfinement().boundingBoxSegments().contains(this->get_P()).first)
                {
                    //                    const VectorDim bbP(std::get<0>(nodeConfinement().boundingBoxSegments().snap(this->get_P())));
                    //                    if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                    //                    {
                    //                        set_P(bbP);
                    _isOnBoundingBox=true;
                    //                    }
                }
                
                
                //                if(nodeConfinement().boundingBoxSegments().size())
                //                {
                //                    const VectorDim bbP(std::get<0>(nodeConfinement().boundingBoxSegments().snap(this->get_P())));
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
            
            //                        std::cout<<"DislocationNode "<<this->sID<<" addLoopLink. _isOnBoundingBox="<<_isOnBoundingBox<<std::endl;
            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            
            //                        std::cout<<"DislocationNode "<<this->sID<<" removeLoopLink"<<std::endl;
            
            // Re-construct nodeConfinement
            _isGlissile=true;
            nodeConfinement().clear();
            
            for(const auto& loopLink : this->loopLinks())
            {
                const bool success = nodeConfinement().addGlidePlane(loopLink->loop()->glidePlane);
                
                if(success)
                {
                    _isGlissile*=loopLink->loop()->isGlissile;
                }
            }
            
            addGrainBoundaryPlanes();
            
            if(nodeConfinement().boundingBoxSegments().contains(this->get_P()).first)
            {
                //                    const VectorDim bbP(std::get<0>(nodeConfinement().boundingBoxSegments().snap(this->get_P())));
                //                    if((this->get_P()-bbP).squaredNorm()<FLT_EPSILON)
                //                    {
                //                        set_P(bbP);
                _isOnBoundingBox=true;
                //                    }
            }
            
            //            if(nodeConfinement().boundingBoxSegments().size())
            //            {// not the last segment being removed
            //                const VectorDim bbP(std::get<0>(nodeConfinement().boundingBoxSegments().snap(this->get_P())));
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
        //            return nodeConfinement().glidePlaneIntersections();
        //        }
        
        //        /**********************************************************************/
        //        const BoundingLineSegments<dim>& boundingBoxSegments() const
        //        {
        //            return nodeConfinement().boundingBoxSegments();
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
                if(_isGrainBoundaryNode)
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
        
        /**********************************************************************/
        const bool& isGrainBoundaryNode() const
        {
            return _isGrainBoundaryNode;
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
            for(const auto& gp : nodeConfinement().glidePlanes())
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
                                if(nodeConfinement().boundingBoxSegments().contains(this->get_P()).first)
                                {// new bounding box contains node
                                    VerboseDislocationNode(2,"case 5"<<std::endl;);
                                    NodeBaseType::set_P(snapToBoundingBox(newP)); // kill numerical errors
                                    _isOnBoundingBox=true;
                                }
                                else
                                {// new bounding box does not contain node
                                    VerboseDislocationNode(2,"case 6"<<std::endl;);
                                    NodeBaseType::set_P(nodeConfinement().snapToGlidePlaneIntersection(newP)); // kill numerical errors
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
                        assert(nodeConfinement().boundingBoxSegments().contains(this->get_P()).first);

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
            
            
            for(const auto& gp : nodeConfinement().glidePlanes())
            {// X must be contained by all glidePlanes
                isMovable*=gp->contains(X);
            }
            
            if(isMovable)
            {
                
//                HERE CHECK THAT X IS IN REGION
                
                for(const auto& pair : this->neighbors())
                {
                    //                            std::cout<<"neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<temp.contains(0.5*(std::get<0>(pair.second)->get_P()+nB->get_P())).first<<std::endl;
                    if(std::get<1>(pair.second)->isBoundarySegment())
                    {// boundary segments other than A->B must remain boundary
                        isMovable*=std::get<1>(pair.second)->boundingBoxSegments().contains(0.5*(std::get<0>(pair.second)->get_P()+X)).first;
                    }
                    
                    if(std::get<1>(pair.second)->isGrainBoundarySegment())
                    {// boundary segments other than A->B must remain boundary
                        for(const auto& gb : std::get<1>(pair.second)->grainBoundaries())
                        {
                            //                                        std::cout<<"grainBoundary neighbor "<<nA->sID<<" "<<pair.first<<" ("<<std::get<0>(pair.second)->sID<<") "<<gb->glidePlanes().begin()->second.contains(nB->get_P())<<std::endl;
                            isMovable*=gb->glidePlanes().begin()->second.contains(X);
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
                const VectorDim newP=nodeConfinement().snapToGlidePlaneIntersection(this->get_P()+dX);
                
                set_P(newP);
                
                //                if(this->network().use_boundary)
                //                {// using confining mesh
                //
                //                    if(_isOnBoundingBox
                //                       //|| nodeConfinement().grains().size()>1
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
