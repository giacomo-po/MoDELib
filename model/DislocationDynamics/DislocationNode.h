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

//#include <model/Network/Operations/EdgeExpansion.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/Geometry/Splines/SplineNode.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/Math/GramSchmidt.h>
//#include <model/DislocationDynamics/DislocationEnergyRules.h>
#include <model/Mesh/Simplex.h>
//#include <model/DislocationDynamics/Junctions/EdgePermutation.h>
#include <model/LatticeMath/LatticeMath.h>
//#include <model/LatticeMath/LineMeshIntersection.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Geometry/SegmentSegmentIntersection.h>


namespace model
{
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationNode : public SplineNode<DislocationNode<_dim,corder,InterpolationType,QuadratureRule>,
    /*                                         */ _dim,corder,InterpolationType>
    
    {
        
    public:
        
        constexpr static int dim=_dim; // make dim available outside class
        typedef DislocationNode       <dim,corder,InterpolationType,QuadratureRule> NodeType;
        typedef DislocationSegment    <dim,corder,InterpolationType,QuadratureRule> LinkType;
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
        typedef 	std::set<VectorDim,
        /*                */ CompareVectorsByComponent<double,dim,float>,
        /*                */ Eigen::aligned_allocator<VectorDim> > ConfiningPlaneIntersectionContainerType;
        typedef std::deque<std::pair<VectorDim,VectorDim>,Eigen::aligned_allocator<std::pair<VectorDim,VectorDim>>> BoundingSegmentsType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        static double bndDistance;
        bool _isGlissile;
        
    private:
        
        DislocationSharedObjects<dim> shared;
        
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        GlidePlaneContainerType _confiningPlanes;
        std::set<const Grain<dim>*> grainSet;
        
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        
        //! The previous velocity vector of *this DislocationNode
        VectorDofType vOld;
        
        double velocityReductionCoeff;
        
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        VectorDim boundaryNormal;
        
        VectorDim C;
        
        BoundingSegmentsType _boundingSegments;
        
        
        /**********************************************************************/
        VectorDim snapToBoundingBox(const VectorDim& P)
        {
            
            std::map<double,VectorDim> snapMap;
            
            for(const auto& vertexPair : _boundingSegments)
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
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
            //            if(guess->region->regionID!=grain.grainID)
            //            {
            //                model::cout<<"DislocationNode "<<this->sID<<std::endl;
            //                model::cout<<"grainID "<<grain.grainID<<std::endl;
            //                model::cout<<"guess= "<<guess->xID<<std::endl;
            //                model::cout<<"guess regionID "<<guess->region->regionID<<std::endl;
            //                assert(0 && "guess does not belong to grain.");
            //            }
            
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (DislocationSharedObjects<dim>::use_boundary)
            {
                if (guess==NULL)
                {
                    temp=DislocationSharedObjects<dim>::mesh.search(this->get_P());
                }
                else
                {
                    temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(this->get_P(),guess);
                }
                
                //                temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(this->get_P(),guess);
                
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
        
        //        /**********************************************************************/
        //        DislocationNode(const LatticeVectorType& Lin,
        //                        const int& grainID,
        //                        const Simplex<dim,dim>* guess=(const Simplex<dim,dim>*) NULL) :
        //        /* base constructor */ NodeBaseType(Lin.cartesian()),
        //        /* init list        */ L(Lin),
        //        /* init list        */ p_Simplex(get_includingSimplex(guess)),
        //        /* init list        */ velocity(VectorDofType::Zero()),
        //        /* init list        */ vOld(velocity),
        //        /* init list        */ velocityReductionCoeff(1.0),
        //        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero())
        //        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        //        {/*! Constructor from DOF
        //          */
        //        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const LatticeVectorType& Lin) :
        //        /* base constructor */ NodeBaseType(pL,Lin.cartesian()),
        /* base constructor */ NodeBaseType(Lin.cartesian()),
        //        /* init list        */ grain(pL.grain),
        //        /* init list        */ L(Lin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity((pL.source->velocity+pL.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
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
            //            std::cout<<"DislocationNode from ExpadingLink A "<<this->sID<<std::endl;
            forceBoundaryNode(pL);
            //            assert(0 && "Initialize C");
        }
        
        /**********************************************************************/
        DislocationNode(const LinkType& pL,
                        const LatticeVectorType& Lin,
                        const VectorDofType& Vin) :
        //        /* base constructor */ NodeBaseType(pL,Lin.cartesian()),
        /* base constructor */ NodeBaseType(Lin.cartesian()),
        //        /* init list        */ grain(pL.grain),
        //        /* init list        */ L(Lin),
        /* init list        */ _isGlissile(true),
        /* init list        */ p_Simplex(get_includingSimplex(pL.source->includingSimplex())),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.source->velocityReduction()+pL.sink->velocityReduction())),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero()),
        //        /* init list        */ grainBoundary_rID2((pL.source->grainBoundary_rID2==pL.sink->grainBoundary_rID2 && pL.sink->grainBoundary_rID2>0)? pL.sink->grainBoundary_rID2 : -1) // TO DO: CHANGE THIS FOR TRIPLE JUNCTIONS
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        //        oldP(this->get_P()),
        //        A1(this->get_P()),
        //        A2(this->get_P())
        C(this->get_P())
        {
            //            std::cout<<"DislocationNode from ExpadingLink B "<<this->sID<<std::endl;
            forceBoundaryNode(pL);
        }
        
        //        /**********************************************************************/
        //        DislocationNode(const ContractingVertices<NodeType,LinkType>& cv,
        //                        const LatticeVectorType& Lin) :
        //        /* base constructor */ NodeBaseType(Lin.cartesian()),
        //        /* init list        */ grain(cv.v0.grain),
        //        /* init list        */ L(Lin),
        //        /* init list        */ p_Simplex(get_includingSimplex(cv.v0.includingSimplex())),
        //        /* init list        */ velocity(0.5*(cv.v0.get_V()+cv.v1.get_V())),
        //        /* init list        */ vOld(velocity),
        //        /* init list        */ velocityReductionCoeff(0.5*(cv.v0.velocityReduction()+cv.v1.velocityReduction())),
        //        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero()),
        //        /* init list        */ grainBoundary_rID2((cv.v0.grainBoundary_rID2==cv.v1.grainBoundary_rID2 && cv.v1.grainBoundary_rID2>0)? cv.v1.grainBoundary_rID2 : -1) // TO DO: CHANGE THIS FOR TRIPLE JUNCTIONS
        //
        //        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        //        {/*! Constructor from VertexContraction
        //          */
        //
        //        }
        
        //        /**********************************************************************/
        //        void make_bndNormal()
        //        {
        //            boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance); // check if node is now on a boundary
        //        }
        
        //        /**********************************************************************/
        //        const LatticeVectorType& get_L() const
        //        {
        //            return L;
        //        }
        
        //        /**********************************************************************/
        //        void set(const LatticeVectorType& Lin)
        //        {
        //            L=Lin;
        //            NodeBaseType::set_P(L.cartesian());
        //        }
        
        //        /**********************************************************************/
        //		void make_planeNormals()
        //        {
        //            //! 1- Clear and re-builds the std::vector planenormals
        //			planenormals.clear();
        //			for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
        //            {
        //				if (std::get<2>(neighborIter->second))
        //                {
        //					LinkType* pL(std::get<1>(neighborIter->second));
        //					if (std::find(planenormals.begin(),planenormals.end(), pL->glidePlaneNormal )==planenormals.end() &&
        //						std::find(planenormals.begin(),planenormals.end(),-pL->glidePlaneNormal )==planenormals.end()   )
        //                    {
        //						planenormals.push_back(pL->glidePlaneNormal );
        //					}
        //					if(pL->sessilePlaneNormal.norm()>FLT_EPSILON)
        //                    {
        //						planenormals.push_back(pL->sessilePlaneNormal);
        //					}
        //
        //				}
        //			}
        //
        //			//! 2- Compute projectionMatrix
        //			make_projectionMatrix();
        //        }
        
        const GlidePlaneContainerType& confiningPlanes() const
        {
            return _confiningPlanes;
        }
        
        
        
        
        
        
        
        
        
        /**********************************************************************/
        BoundingSegmentsType updateBoundingSegments(const BoundingSegmentsType& old,
                                                    const GlidePlaneType& gp) const
        {
            //model::cout<<"DislocationNode "<<this->sID<<" adding GlidePlane "<<gp.sID<<std::endl;
            
            
            BoundingSegmentsType temp;
            
            if(old.size())
            {
                
                ConfiningPlaneIntersectionContainerType uniquePts;
                
                
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
                                uniquePts.insert(ssi.x0);
                                uniquePts.insert(ssi.x1);
                                
                            }
                            
                            
                            
                        }
                    }
                    
                }
                
                switch (uniquePts.size())
                {
                    case 1:
                    {
                        temp.emplace_back(*uniquePts.begin(),*uniquePts.begin());
                        break;
                    }
                    case 2:
                    {
                        temp.emplace_back(*uniquePts.begin(),*uniquePts.rbegin());
                        break;
                    }
                        
                    default:
                    {
                        model::cout<<"DislocationNode "<<this->sID<<" uniquePts.size()="<<uniquePts.size()<<std::endl;
                        for(const auto& pt : uniquePts)
                        {
                            model::cout<<pt.transpose()<<std::endl;
                        }
                        assert(0);
                        break;
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
            
            assert(temp.size()==1 || temp.size()>=3 && "updateBoundingSegments failed");
            return temp;
            
        }
        
        
        /**********************************************************************/
        void addLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::addLoopLink
          */
            NodeBaseType::addLoopLink(pL); // forward to base class
            
            // Insert new plane in _confiningPlanes. If plane already exists nothing will happen
            const bool success=_confiningPlanes.insert(&(pL->loop()->glidePlane)).second;
            if(success)
            {
                _isGlissile*=pL->loop()->isGlissile;
                // Update _boundingSegments
                _boundingSegments=updateBoundingSegments(_boundingSegments,pL->loop()->glidePlane);
                
                // Insert new grain in grainSet
                grainSet.insert(&(pL->loop()->grain));
            }
            
        }
        
        /**********************************************************************/
        void removeLoopLink(LoopLinkType* const pL)
        {/*@param[in] pL LoopLink pointer
          *
          * This functin overrides LoopNode::removeLoopLink
          */
            NodeBaseType::removeLoopLink(pL); // forward to base class
            
            // Re-construct _confiningPlanes and grainSet
            _isGlissile=true;
            _confiningPlanes.clear();
            _boundingSegments.clear();
            grainSet.clear();
            for(const auto& loopLink : this->loopLinks())
            {
                _isGlissile*=loopLink->loop()->isGlissile;
                _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
                grainSet.insert(&(loopLink->loop()->grain));
                _boundingSegments=updateBoundingSegments(_boundingSegments,loopLink->loop()->glidePlane);
            }
            
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
                //
                //                make_bndNormal();
            }
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
        int meshLocation() const
        {/*!\returns the position of *this relative to the bonudary:
          * 1 = inside mesh
          * 2 = on mesh boundary
          */
            return (boundaryNormal.squaredNorm()>FLT_EPSILON? onMeshBoundary : insideMesh);
        }
        
        /**********************************************************************/
        bool isBoundaryNode() const
        {
            return meshLocation()==onMeshBoundary;
        }
        
        /**********************************************************************/
        bool isGrainBoundaryNode() const
        {
            return grainSet.size()>1;
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
            
//            switch (_confiningPlanes.size())
//            {
//                case 1:
//                    //                    dX=_confiningPlanes[0]->n.snapToLattice(dX).cartesian();
//                    dX=(*_confiningPlanes.begin())->n.snapToPlane(dX);
//                    break;
//                    
//                case 2:
//                {
//                    //                    PlanePlaneIntersection ppi(glidePlane(0),glidePlane(1));
//                    const LatticeDirection<dim> d(glidePlane(0).n.cross(glidePlane(1).n));
//                    //LatticeLine line(ppi.P,ppi.d);
//                    //                    dX=line.d.snapToDirection(dX).cartesian();
//                    dX=d.snapToDirection(dX);
//                    break;
//                }
//                    
//                default:
//                    dX.setZero();
//                    break;
//            }
            
            
            if (dX.squaredNorm()>0.0 && _isGlissile) // move a node only if |v|!=0
            {
                
                
                
                if(shared.use_boundary) // using confining mesh
                {
                    
                    if(_boundingSegments.size()>1)
                    {
                        
                        set_P(this->get_P()+dX);

                        
                        //                    // See if the new position is inside mesh
                        //                    std::set<const Simplex<dim,dim>*> path;
                        //                    const bool searchAllRegions=true; // CHANGE THIS FOR MULTIPLE REGIONS
                        //                    std::pair<bool,const Simplex<dim,dim>*> temp(DislocationSharedObjects<dim>::mesh.searchWithGuess(searchAllRegions,
                        //                                                                                                                     this->get_P()+dX,
                        //                                                                                                                     p_Simplex,
                        //                                                                                                                     path));
                        
                    }
                    else
                    {
                        set_P(snapToBoundingBox(this->get_P()+dX));
                    }
                    
                    
                    make_projectionMatrix();
                    
                }
                else // move node freely
                {
                    set_P(this->get_P()+dX);
                    //set(L+grain.latticeVector(dX));
                    //                    L+=LatticeVectorType(dX);
                    //                    this->set(this->get_nodeDof()+dX);
                }
                
                
                //                REMEBER TO UPDATE p_Simplex and boundaryNormal
                
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
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.get_P().transpose()<<"\t"
            /**/<< ds.get_V().transpose()<<"\t"
            /**/<< ds.velocityReduction()<<"\t"
            /**/<< ds.pSN()->sID<<"\t"
            /**/<< (ds.meshLocation()==onMeshBoundary);
            return os;
        }
        
    };
    
    
    // static data
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    bool DislocationNode<_dim,corder,InterpolationType,QuadratureRule>::use_velocityFilter=true;
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    double DislocationNode<_dim,corder,InterpolationType,QuadratureRule>::velocityReductionFactor=0.75;
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    double DislocationNode<_dim,corder,InterpolationType,QuadratureRule>::bndDistance=FLT_EPSILON;
    
    
}
#endif


//        /**********************************************************************/
//        void forceBoundaryNode(const LinkType& pL)
//        {
//
//            assert(0 && "RE_ENABLE forceBoundaryNode");
//            //            if(shared.use_boundary)
//            //            {
//            //                if(pL.is_boundarySegment() && !isBoundaryNode())
//            //                {
//            //                    const VectorDim nb=(pL.source->bndNormal()+pL.sink->bndNormal()).normalized();
//            //                    const VectorDim np=pL.glidePlane->n.cartesian().normalized();
//            //                    VectorDim outDir=nb-nb.dot(np)*np;
//            //                    if(outDir.squaredNorm()>FLT_EPSILON)
//            //                    {
//            //                        //                        std::cout<<"DislocationNode "<<this->sID<<", outDir="<<outDir.transpose()<<std::endl;
//            //                        //                        std::cout<<pL.source->get_P().transpose()<<std::endl;
//            //                        //                        std::cout<<pL.sink->get_P().transpose()<<std::endl;
//            //                        //                        std::cout<<this->get_P().transpose()<<std::endl;
//            //
//            //                        outDir.normalize();
//            //                        LatticeVectorType dL(pL.glidePlane->n.snapToLattice(outDir));
//            //                        assert(dL.squaredNorm()>0.0);
//            //
//            //                        LatticeVectorType L0=L;
//            //                        if(!DislocationSharedObjects<dim>::mesh.searchWithGuess(L0.cartesian(),p_Simplex).first)
//            //                        {
//            //                            L0=LatticeVectorType(pL.glidePlane->snapToLattice(0.5*(pL.source->get_P()+pL.sink->get_P())));
//            //                        }
//            //                        assert(DislocationSharedObjects<dim>::mesh.searchWithGuess(L0.cartesian(),p_Simplex).first && "L0 outside mesh");
//            //
//            //                        LatticeLine line(L0,dL);
//            //                        LineMeshIntersection lmi(line,L+dL,shared.mesh,p_Simplex);
//            //                        if(lmi.search.first)
//            //                        {
//            //                            if(lmi.search.second->region->regionID==grain.grainID)
//            //                            {
//            //                                p_Simplex=lmi.search.second;
//            //                                set(lmi.L);
//            //                                make_bndNormal();
//            //                                //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance);
//            //                                if(!isBoundaryNode())
//            //                                {
//            //                                    //                                std::cout<<"DislocaitonNode "<<this->sID<<" not on mesh boundary"<<std::endl;
//            //                                    //                                std::cout<<"dL="<<dL.transpose()<<std::endl;
//            //                                }
//            //                                assert(isBoundaryNode());
//            //                            }
//            //                            else
//            //                            {
//            //                                assert(0);
//            //                            }
//            //
//            //                        }
//            //                        else
//            //                        {
//            //                            assert(0 && "lmi failed.");
//            //                        }
//            //                    }
//            //                    else
//            //                    {
//            //                        model::cout<<"DislocationNode "<<this->sID<<std::endl;
//            //                        assert(0 && "outDir has zero norm");
//            //                    }
//            //                }
//            //            }
//        }
//
//        /**********************************************************************/
//        void moveToBoundary(const VectorDim& outDir,const std::pair<bool,
//                            const Simplex<dim,dim>*>& temp,
//                            const VectorDim& dX)
//        {
//
//            assert(0 && "RE-ENABLE moveToBoundary");
//
//            //            if(outDir.squaredNorm()<FLT_EPSILON)
//            //            {
//            //                model::cout<<"DislocationNode "<<this->sID<<" outDir has zero norm"<<std::endl;
//            //                assert(0 && "outDir has zero norm");
//            //            }
//            //
//            //            LatticeVectorType dL(grain.latticeVector(VectorDim::Zero()));
//            //            switch (_confiningPlanes.size())
//            //            {
//            //                case 1:
//            //                {
//            //                    //std::cout<<"boundary motion case 1, DislocationNode "<<this->sID<<std::endl;
//            //                    const VectorDim planeN(_confiningPlanes[0]->n.cartesian().normalized());
//            //                    VectorDim dD=outDir-outDir.dot(planeN)*planeN;
//            //                    const double dDnorm=dD.norm();
//            //                    if(dDnorm>0.0)
//            //                    {
//            //                        dD.normalize();
//            //                        dL=LatticeVectorType(_confiningPlanes[0]->n.snapToLattice(dD)); // a lattice vector on the plane pointing ouside mesh
//            //                        assert(dL.squaredNorm()>0);
//            //                    }
//            //                    break;
//            //                }
//            //
//            //                case 2:
//            //                {
//            //                    //std::cout<<"boundary motion case 2, DislocationNode "<<this->sID<<std::endl;
//            //                    const PlanePlaneIntersection ppi(glidePlane(0) ,glidePlane(1) );
//            //                    //dL=LatticeVectorType(ppi.d.snapToDirection(10.0*outDir));
//            //                    dL=ppi.d.snapToDirection(10.0*outDir);
//            //                    break;
//            //                }
//            //
//            //            }
//            //
//            //            if(dL.squaredNorm()>FLT_EPSILON) // a line direction could be found
//            //            {
//            //                const LatticeVectorType L0(grain.latticeVector((this->get_P()+dX*temp.first).eval()));
//            //                const LatticeLine line(L0,dL);
//            //                LineMeshIntersection lmi(line,L0+dL,shared.mesh,temp.second);
//            //                assert(lmi.search.first);
//            //                if(lmi.search.second->region->regionID==grain.grainID)
//            //                {
//            //                    p_Simplex=lmi.search.second;
//            //                    set(lmi.L);
//            //                    make_bndNormal();
//            //                    //                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance);
//            //                    if(meshLocation()!=onMeshBoundary)
//            //                    {
//            //                        model::cout<<"DislocaitonNode "<<this->sID<<std::endl;
//            //                        assert(0 && "NODE MUST BE ON MESH-BOUNDARY");
//            //                    }
//            //                }
//            //                else
//            //                {
//            //                    iterateToBoundary(dX,lmi.search);
//            //                }
//            //            }
//
//        }
//
//
//
//
//        /**********************************************************************/
//        void iterateToBoundary(	const VectorDim& dX,
//                               const std::pair<bool,const Simplex<dim,dim>*>& temp)
//        {
//
//            assert(0 && "RE-ENABLE iterateToBoundary");
//
//            //            LatticeDirectionType dXprimitive(grain.latticeDirection(dX));
//            //            grainBoundary_rID2=temp.second->region->regionID;
//            //            assert(temp.second->region->regionID!=grain.grainID);
//            //            for(int i=1;(dX-dXprimitive.cartesian()*i).norm()>0;i++)
//            //            {
//            //                std::pair<bool,const Simplex<dim,dim>*>  temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(this->get_P()+dX-dXprimitive.cartesian(),p_Simplex);
//            //                if(temp.first&&temp.second->region->regionID==grain.grainID)
//            //                {
//            //                    model::cout<<"A boundary node was crossing the region boundary. Instead, the node is moved as close as possible to the interface"<<std::endl;
//            //                    p_Simplex=temp.second;
//            //                    set(L+grain.latticeVector(dX-dXprimitive.cartesian()));
//            //
//            //                    if(shared.poly.grainBoundary(grain.grainID,grainBoundary_rID2).latticePlane(grain.grainID).contains(this->get_P()))
//            //                    {
//            //                        assert(grainBoundaryPlane().contains(this->get_P()));
//            //                    }
//            //                    else
//            //                    {
//            //                        grainBoundary_rID2=-1;
//            //                    }
//            //                    break;
//            //                }
//            //                if(i>25)
//            //                {
//            //                    break;
//            //                }
//            //
//            //            }
//        }



//move()
//{
//                    // See if the new position is inside mesh
//                    std::set<const Simplex<dim,dim>*> path;
//                    const bool searchAllRegions=true; // CHANGE THIS FOR MULTIPLE REGIONS
//                    std::pair<bool,const Simplex<dim,dim>*> temp(DislocationSharedObjects<dim>::mesh.searchWithGuess(searchAllRegions,
//                                                                                                                     this->get_P()+dX,
//                                                                                                                     p_Simplex,
//                                                                                                                     path));
//                    std::cout<<"IS THIS NECESSARY WITH THE NEW CONFINING PERIMETER?"<<std::endl;
//                    //p_Simplex=temp.second;
//
//
//                    if(isBoundaryNode()) // node already a boundary node, it must remain on boundary
//                    {
//                        assert(0 && "RE_ENABLE THIS 1");
//                        //                        const VectorDim bndNrml=SimplexBndNormal::get_boundaryNormal(this->get_P()+dX,*temp.second,bndDistance); // boundary normal at new position
//                        //                        if(temp.first && bndNrml.squaredNorm()>0.0) // new position is on boundary
//                        //                        {
//                        //                            if(temp.second->region->regionID!=grain.grainID)
//                        //                            {
//                        //                                iterateToBoundary(dX,temp);
//                        //                            }
//                        //                            else
//                        //                            {
//                        //                                p_Simplex=temp.second;
//                        //                                set(L+grain.latticeVector(dX));
//                        //                                //boundaryNormal=bndNrml;
//                        //                                make_bndNormal();
//                        //                                assert(meshLocation()==onMeshBoundary);
//                        //                           	}
//                        //                        }
//                        //                        else // new position is not on boundary
//                        //                        {
//                        //                            moveToBoundary(boundaryNormal,temp,dX);
//                        //                        }
//
//
//                    }
//                    else // not a boundary node
//                    {
//                        if(!temp.first) // node moved outside or already on boundary
//                        {
//                            assert(0 && "RE_ENABLE THIS 2");
//
//                            //                            VectorDim outDir=boundaryNormal;
//                            //                            if(outDir.squaredNorm()==0.0)
//                            //                            {// node is exiting for the first time, we need a tentative boundary normal to identify the "outside direction"
//                            //                                for(const auto& simplex : path) // loop over the pathof simplices  connecting P to P+dX
//                            //                                {
//                            //                                    outDir=SimplexBndNormal::get_boundaryNormal(this->get_P()+dX,*simplex,dX.norm()+FLT_EPSILON);
//                            //
//                            //                                    if(outDir.squaredNorm()>0.0)
//                            //                                    {
//                            //                                        break;
//                            //                                    }
//                            //
//                            //                                }
//                            //
//                            //                            }
//                            //                            assert(outDir.squaredNorm()>0 && "COULD NOT DETERMINE OUTDIR");
//                            //                            moveToBoundary(outDir,temp,dX);
//                        }
//                        else // node is internal and remains internal
//                        {
//                            if(false /*temp.second->region->regionID!=grain.grainID*/)
//                            {// node is crossing regions
//                                assert(0 && "RE_ENABLE THIS 2");
//                                //
//                                //                                grainBoundary_rID2=temp.second->region->regionID;
//                                //
//                                //                                const LatticeDirectionType dL(grain.latticeVector(dX));
//                                //                                LatticeLine trajectory(L,dL);
//                                //                                PlaneLineIntersection pli(grainBoundaryPlane(),trajectory);
//                                //
//                                //                                switch (pli.intersectionType)
//                                //                                {
//                                //                                    case PlaneLineIntersection::coincident:
//                                //                                    {
//                                //                                        set(L+grain.latticeVector(dX));
//                                //                                        break;
//                                //                                    }
//                                //                                    case PlaneLineIntersection::intersecting:
//                                //                                    {
//                                //                                        set(pli.P);
//                                //                                        break;
//                                //                                    }
//                                //                                    case PlaneLineIntersection::offLattice:
//                                //                                    {
//                                //                                        switch (_confiningPlanes.size())
//                                //                                        {
//                                //                                            case 1:
//                                //                                            {
//                                //                                                PlanePlaneIntersection ppi(grainBoundaryPlane(),glidePlane(0) );
//                                //                                                //model::cout<<redColor<<"About to check for offlattice!!"<<Color<<std::endl;
//                                //                                                if(ppi.intersectionType==3)
//                                //                                                {
//                                //                                                    std::cout<<"DETECTED NON-LATTICE INTERSECTION WITH GB PLANE. ITERATIVELY MOVING NODE "<<this->sID<<" TOWARDS GB"<<std::endl;
//                                //                                                    iterateToBoundary(dX,temp);
//                                //                                                    //assert(0 && "ppi intersection found, but off lattice!");
//                                //                                                    }
//                                //                                                    LatticeLine gbLine(ppi.P,ppi.d);
//                                //                                                    set(gbLine.snapToLattice(pli.P));
//                                //                                                    break;
//                                //                                                    }
//                                //                                                case 2:
//                                //                                                    {
//                                //                                                        std::cout<<"DETECTED NON-LATTICE INTERSECTION WITH 2 CONFINING PLANES"<<std::endl;
//                                //                                                        break;
//                                //                                                    }
//                                //                                                default:
//                                //                                                    break;
//                                //                                                    }
//                                //                                                    break;
//                                //                                                    }
//                                //
//                                //                                                default:
//                                //                                                    assert(0 && "parallel intersection detected");
//                                //                                                    break;
//                                //                                                    }
//                                //
//                                //                                                    temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(false,
//                                //                                                                                                             this->get_P(),
//                                //                                                                                                             p_Simplex,
//                                //                                                                                                             path);
//                                //
//                                //                                                    if(!(temp.first && temp.second->region->regionID==grain.grainID))
//                                //                                                    {
//                                //                                                        model::cout<<"DislocationNode "<<this->sID<<std::endl;
//                                //                                                        model::cout<<"GrainID= "<<grain.grainID<<std::endl;
//                                //                                                        model::cout<<"search found region= "<<temp.second->region->regionID<<std::endl;
//                                //                                                        model::cout<<"search start="<<p_Simplex->xID<<std::endl;
//                                //                                                        model::cout<<"search end="<<temp.second->xID<<std::endl;
//                                //                                                        model::cout<<"bary coords="<<temp.second->pos2bary(this->get_P())<<std::endl;
//                                //                                                        assert(0 && "Intersection point of GB-node not found in correct region");
//                                //                                                    }
//                                //                                                    p_Simplex=temp.second;
//                                //                                                    make_confiningPlanes();
//                                //
//                                //
//                                //
//                                //                                                    //assert(0 && "RE-ENABLE THIS");
//                                //                                                    //                            //// ////std::cout<<"DislocationNode "<<this->sID<<" crossing region. Path size="<<path.size()<<std::endl;
//                                //                                                    //
//                                //                                                    //                            int faceID=-1;
//                                //                                                    ////                            const Simplex<dim,dim>* regBndSimplex=(const Simplex<dim,dim>*) NULL;
//                                //                                                    //                            const Simplex<dim,dim>* regBndSimplex(NULL);
//                                //                                                    //
//                                //                                                    //                            Eigen::Matrix<double,dim+1,1> faceInt(Eigen::Matrix<double,dim+1,1>::Zero());
//                                //                                                    //
//                                //                                                    //                            for(const auto& simplex : path) // loop over the pathof simplices  connecting P to P+dX
//                                //                                                    //                            {
//                                //                                                    //                                const Eigen::Matrix<double,dim+1,1> baryOld(simplex->pos2bary(this->get_P()));
//                                //                                                    //                                const Eigen::Matrix<double,dim+1,1> baryNew(simplex->pos2bary(this->get_P()+dX));
//                                //                                                    //                                for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
//                                //                                                    //                                {
//                                //                                                    //                                    if(simplex->child(f).isRegionBoundarySimplex())
//                                //                                                    //                                    {
//                                //                                                    //                                        faceInt=simplex->faceLineIntersection(baryOld,baryNew,f);
//                                //                                                    //                                        //                                    //////std::cout<<"DislocationNode "<<this->sID<<", baryMin="<<faceInt.minCoeff()<<std::endl;
//                                //                                                    //
//                                //                                                    //                                        if(faceInt.minCoeff()>=-FLT_EPSILON) // faceInt belongs to triangle
//                                //                                                    //                                        {
//                                //                                                    //                                            regBndSimplex=simplex; // current simplex is the region boundary simplex wanted
//                                //                                                    //                                            faceID=f; // intersection face is f
//                                //                                                    //                                            break;
//                                //                                                    //                                        }
//                                //                                                    //                                    }
//                                //                                                    //                                }
//                                //                                                    //
//                                //                                                    //                                if(faceID>=0)
//                                //                                                    //                                {
//                                //                                                    //                                    break;
//                                //                                                    //                                }
//                                //                                                    //                            }
//                                //                                                    //
//                                //                                                    //                            assert(faceID>=0 && "FACE INTERSECTION NOT FOUND");
//                                //                                                    //                            this->set(regBndSimplex->bary2pos(faceInt)); // move node to intersction with region boundary
//                                //                                                    //                            regionBndNormal=regBndSimplex->nda.col(faceID).normalized();
//                            }
//                            else // node not crossing regions
//                            {
//                                //p_Simplex=temp.second;
//                                set_P(this->get_P()+dX);
//                                //set(L+grain.latticeVector(dX));
//                                //make_bndNormal();
//                                //                                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance); // check if node is now on a boundary
//
//                                //                                L+=LatticeVectorType(dX);
//                                //                                this->set(this->get_P()+dX); // move node
//                                //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,10.0); // check if node is now on a boundary
//                                //
//                                //                            PROBLEM HERE, BOUNDSARY NODES ARE ASSIGNED ZERO NORMAL
//
//                            }
//
//                        }
//                    }
//}

//        void storeP()
//        {
//            oldP=this->get_P();
//        }
//
//        void stepA1(const double& dt)
//        {
//            oldP=this->get_P();
//            A1=this->get_P()+velocity*dt;
//            set_P(this->get_P()+velocity*dt*0.5);
//        }
//
//        double stepA2(const double& dt)
//        {
//            A2=this->get_P()+velocity*dt*0.5;
////            set_P(this->get_P()+velocity*dt*0.5);
//            set_P(oldP);
//
//            return (A1-A2).norm();
//        }

//        /**********************************************************************/
//        void make_confiningPlanes()
//        {
//
////            std::cout<<"make_confiningPlanes 0"<<std::endl;
//            _isGlissile=true;
//            _confiningPlanes.clear();
//            for(const auto& loopLink : this->loopLinks())
//            {
//
//                if(_confiningPlanes.find(&loopLink->loop()->glidePlane)==_confiningPlanes.end())
//                {// new GlidePlane
//                    if(_confiningPlanes.size()>1)
//                    {// Remove redunant planes
//
//                        if(glidePlane(0).grain.grainID==glidePlane(1).grain.grainID)
//                        {// first two planes in same lattice
//                            const LatticeDirection<dim> d=glidePlane(0).n.cross(glidePlane(1).n);
//                            assert(d.squaredNorm()>0 && "Planes are parallel or the same");
//
//                            if(fabs(loopLink->loop()->glidePlane.n.cartesian().normalized().dot(d.cartesian().normalized()))>FLT_EPSILON)
//                            {
//                                _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
//                            }
//
//                        }
//                        else
//                        {// first two planes in different lattices
//                            const VectorDim d=glidePlane(0).n.cartesian().normalized().cross(glidePlane(1).n.cartesian().normalized());
//                            const double dNorm(d.norm());
//                            if(dNorm>FLT_EPSILON)
//                            {
//                                if(fabs(loopLink->loop()->glidePlane.n.cartesian().normalized().dot(d/dNorm))>FLT_EPSILON)
//                                {
//                                    _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
//                                }
//                            }
//                            else
//                            {// glidePlane(0) and glidePlane(1) are coincident and in different grains
//                                if(loopLink->loop()->glidePlane.n.cartesian().normalized().cross(glidePlane(0).n.cartesian().normalized()).norm()<FLT_EPSILON)
//                                {
//                                    assert(loopLink->loop()->grain.grainID!=glidePlane(0).grain.grainID && loopLink->loop()->grain.grainID!=glidePlane(1).grain.grainID);
//                                    _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
//
//                                }
//                                else
//                                {
//                                    const int erased=_confiningPlanes.erase(&glidePlane(0));
//                                    assert(erased==1);
//                                    _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
//                                }
//
//                            }
//
//
//                        }
//
//                    }
//                    else
//                    {
//                        _confiningPlanes.insert(&(loopLink->loop()->glidePlane));
//                    }
//
//                }
//
//                if(loopLink->loop()->isSessile || _confiningPlanes.size()>2)
//                {
//                    _isGlissile=false;
//                }
//
//            }
//
//            assert(_confiningPlanes.size());
//
//
//            confiningPlaneIntersections.clear();
//            if(_confiningPlanes.size()==2)
//            {
//
////                            std::cout<<"make_confiningPlanes 3"<<std::endl;
//                const auto& mi0 = glidePlane(0).meshIntersections;
//                const auto& mi1 = glidePlane(1).meshIntersections;
//
//                for(size_t i0=0;i0<mi0.size();++i0)
//                {
//                    const size_t j0= (i0==mi0.size()-1)? 0 : i0+1;
//
//                    const VectorDim& A0(mi0[i0].second);
//                    const VectorDim& B0(mi0[j0].second);
//
//                    for(size_t i1=0;i1<mi1.size();++i1)
//                    {
//
////                                    std::cout<<"make_confiningPlanes 4"<<std::endl;
//                        const size_t j1= (i1==mi1.size()-1)? 0 : i1+1;
//
//                        const VectorDim& A1(mi1[i1].second);
//                        const VectorDim& B1(mi1[j1].second);
//
////                        const std::deque<VectorDim,Eigen::aligned_allocator<VectorDim>> temp=segmentSegmentIntersection(A0, B0, A1, B1);
//                        for(const auto& x : segmentSegmentIntersection(A0, B0, A1, B1))
//                        {
////                                        std::cout<<"make_confiningPlanes 5"<<std::endl;
////                            std::cout<<x.transpose()<<std::endl;
//                            confiningPlaneIntersections.insert(x);
////                            std::cout<<"make_confiningPlanes 5a"<<std::endl;
//
//                        }
//
//                    }
//                }
////                            std::cout<<"make_confiningPlanes 6"<<std::endl;
//            }
//            else
//            {
//
//            }
//
//            //            // add to _confiningPlanes the special planes of this node
//            //            for(const auto& plane : specialConfiningPlanes)
//            //            {
//            //                _confiningPlanes.push_back(&plane);
//            //            }
//
////                        std::cout<<"make_confiningPlanes 7"<<std::endl;
//
//            //std::cout<<"a"<<std::endl;
//
//
//            //std::cout<<"b"<<std::endl;
//
//            //! 2- Compute projectionMatrix
//            make_projectionMatrix();
//
////            std::cout<<"make_confiningPlanes 8"<<std::endl;
//
//
//            //std::cout<<"done"<<std::endl;
//
//        }


//        /**********************************************************************/
//        void make_grains()
//        {
//            grainSet.clear();
//            for(const auto& loopLink : this->loopLinks())
//            {
//                grainSet.insert(&(loopLink->loop()->grain));
//                //                if(loopLink->loop()->isSessile)
//                //                {
//                //                    assert(0 && "FINISH HERE, ADD MORE PLANES TO FULLY CONSTRAIN NODE");
//                //                }
//            }
//        }


//        /**********************************************************************/
//        bool invertedMotion() const
//        {/*! The nodal velocity vector
//          */
//            return velocity.template segment<dim>(0).dot( vOld.template segment<dim>(0) ) < 0.0;
//        }

//        /**********************************************************************/
//        const LatticePlane& grainBoundaryPlane() const
//        {
//            assert(grainBoundary_rID2>0 && "Node not on Grain Boundary");
//            return shared.poly.grainBoundary(grain.grainID,grainBoundary_rID2).latticePlane(grain.grainID);
//        }


//        /**********************************************************************/
//        std::pair<std::deque<std::pair<size_t,size_t> >,std::deque<std::pair<size_t,size_t> >> edgeDecomposition() const
//        {
//
//            std::deque<std::pair<size_t,size_t> > firstLinkDeq;
//            std::deque<std::pair<size_t,size_t> > secondLinkDeq;
//
//            if(this->Neighborhood.size()>2 && meshLocation()==insideMesh)
//            {
//                std::deque<std::pair<int,int> > linkDeq;
//
//
//                const size_t neighSize=this->Neighborhood.size()-1;
//
//                Eigen::MatrixXd temp(neighSize,dim);
//
//                int k=0;
//                for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
//                {
//                    if (std::get<2>(neighborIter->second)) // not self
//                    {
//                        LinkType* pL(std::get<1>(neighborIter->second));
//                        temp.row(k)=pL->Burgers;
//                        linkDeq.emplace_back(pL->source->sID,pL->sink->sID);
//                        k++;
//                    }
//                }
//
//                //                if(temp.rows()>7)
//                //                {
//                //                    model::cout<<"Dislocation Node "<<this->sID<<std::endl;
//                //                }
//                auto vecVec=EdgePermutations::edgeStats(temp);
//                if(vecVec.size()==2)
//                {
//                    if((vecVec[0].array()*vecVec[1].array()).matrix().squaredNorm()<FLT_EPSILON) // decompisition is unique
//                    {
//                        //                       //std::cout<<"DislocationNode "<<this->sID<<std::endl;
//                        //                       //std::cout<<"Matrix of Burgers= "<<std::endl<<temp<<std::endl;
//                        //                       //std::cout<<"null configurations= "<<std::endl;
//                        //                        for(auto v : vecVec)
//                        //                        {
//                        //                           //std::cout<<v<<std::endl;
//                        //                        }
//                        for(size_t n=0;n<neighSize;++n)
//                        {
//                            if(vecVec[0](n)!=0)
//                            {
//                                firstLinkDeq.push_back(linkDeq[n]);
//                            }
//                            else
//                            {
//                                secondLinkDeq.push_back(linkDeq[n]);
//                            }
//                        }
//
//                    }
//                }
//            }
//
//            return make_pair(firstLinkDeq,secondLinkDeq);
//        }


//        /**********************************************************************/
//        LatticeVectorType zeroFlow() const
//        {
//            return LatticeVectorType(grain.lattice());
//        }

//        /**********************************************************************/
//        const int& rID2() const
//        {
//            return grainBoundary_rID2;
//        }

//        /**********************************************************************/
//        void removeFromNeighborhood(LinkType* const pL)
//        {/*!@param[in] pL A pointer to the DislocationSegment being disconnected
//          * from this node.
//          *
//          *  Overwrites NetworkNode::removeFromNeighborhood in order to modify
//          *  planeNormals and tangent after the DislocationSegment is disconnected
//          */
//            NodeBaseType::removeFromNeighborhood(pL);
//            //            make_planeNormals();
//            make_confiningPlanes();
//            //            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
////            NodeBaseType::make_T();
//        }
//
//        /**********************************************************************/
//        void addToNeighborhood(LinkType* const pL)
//        {/*!@param[in] pL A pointer to the DislocationSegment being connected
//          * to this node.
//          *
//          *  Overwrites NetworkNode::addToNeighborhood in order to modify
//          *  planeNormals and tangent after the DislocationSegment is disconnected
//          */
//            NodeBaseType::addToNeighborhood(pL);
//            //            make_planeNormals();
//            make_confiningPlanes();
//            //            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
////            NodeBaseType::make_T();
//        }


//        /**********************************************************************/
//		VectorOfNormalsType constraintNormals() const
//        {
//			GramSchmidt<dim> GS;
//			if (meshLocation()==insideMesh)
//            { // DislocationNode is inside mesh
//				if (this->is_balanced())
//                { // DislocationNode is balanced
//					GS.orthoNormalize(planenormals); //  constrain by planenormals
//				}
//				else
//                { // DislocationNode is unbalanced and is inside mesh: fix
//					GS.push_back((VectorDim()<<1.0,0.0,0.0).finished());
//					GS.push_back((VectorDim()<<0.0,1.0,0.0).finished());
//					GS.push_back((VectorDim()<<0.0,0.0,1.0).finished());
//				}
//			}
//			else if (meshLocation()==onMeshBoundary)
//            { // DislocationNode is on mesh boundary:
//				VectorOfNormalsType PN(planenormals); // constrain by planenormals
//				PN.push_back(boundaryNormal); // constrain by boundaryNormal
//				GS.orthoNormalize(PN);
//			}
//			else
//            {
//               //std::cout<<"DislocationNode "<<this->sID<< " at "<<this->get_P().transpose()<<" is outside mesh."<<std::endl;
//                assert(0 && "DISLOCATION NODE FOUND OUTSIDE MESH."); //RE-ENABLE THIS
//			}
//			assert(GS.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
//			return GS;
//		}

//        void //showTemp(const VectorOfNormalsType& temp) const
//        {
//            for(int k=0;k<temp.size();++k)
//            {
//               //std::cout<<temp[k].transpose()<<std::endl;
//            }
//        }

//                            const LatticeVectorType dL(outdX(outDir));

//                        LatticeVectorType dL(LatticeVectorType::Zero());
//                        switch (_confiningPlanes.size())
//                        {
//                            case 1:
//                            {
//                                //std::cout<<"boundary motion case 1, DislocationNode "<<this->sID<<std::endl;
//                                const VectorDim planeN(_confiningPlanes[0]->n.cartesian().normalized());
//                                VectorDim dD=outDir-outDir.dot(planeN)*planeN;
//                                const double dDnorm=dD.norm();
//                                if(dDnorm>0.0)
//                                {
//                                    dD.normalize();
//                                    dL=LatticeVectorType(_confiningPlanes[0]->n.snapToLattice(dD)); // a lattice vector on the plane pointing ouside mesh
//                                    assert(dL.squaredNorm()>0);
//                                }
//                                break;
//                            }
//
//                            case 2:
//                            {
//                                //std::cout<<"boundary motion case 2, DislocationNode "<<this->sID<<std::endl;
//                                const PlanePlaneIntersection ppi(glidePlane(0) ,glidePlane(1) );
//                                dL=LatticeVectorType(ppi.d.snapToDirection(10.0*outDir));
//                                break;
//                            }
//
//                        }

//                            if(dL.squaredNorm()) // a line direction could be found
//                            {
//                                const LatticeLine line(L0,dL);
//                                LineMeshIntersection lmi(line,L0+dL,shared.mesh,temp.second);
//                                assert(lmi.search.first);
//                                p_Simplex=lmi.search.second;
//                                set(lmi.L);
//                                //                                L=lmi.L;
//                                //                                this->set(L.cartesian()); // move node
//                                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,10.0);
//                                assert(meshLocation()==onMeshBoundary);
//                            }
//                            else
//                            {
//                                if(temp.first) // only for boundary nodes
//                                {
//                                    const VectorDim bndNrml=SimplexBndNormal::get_boundaryNormal(this->get_P()+dX,*temp.second,10.0);
//                                    if(bndNrml.squaredNorm()>0.0)
//                                    {
//                                        p_Simplex=temp.second;
//                                        set(L+LatticeVectorType(dX));
//                                        boundaryNormal=bndNrml;
//                                        assert(meshLocation()==onMeshBoundary);
//                                    }
//                                }
//                            }
