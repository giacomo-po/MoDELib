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

#include <model/Network/Operations/EdgeExpansion.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/Geometry/Splines/SplineNode.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/Math/GramSchmidt.h>
#include <model/DislocationDynamics/DislocationEnergyRules.h>
#include <model/Mesh/Simplex.h>
#include <model/DislocationDynamics/Junctions/EdgePermutation.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/LatticeMath/LineMeshIntersection.h>
#include <model/DislocationDynamics/SimplexBndNormal.h>

namespace model
{
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationNode : public SplineNode<DislocationNode<_dim,corder,InterpolationType,QuadratureRule>,
    /*                                         */ _dim,corder,InterpolationType>
    
    {
        
    public:
        
        // make dim available outside class
        constexpr static int dim=_dim;
        
        
        // define Derived to use NetworkTypedefs.h
        typedef DislocationNode       <dim,corder,InterpolationType,QuadratureRule> Derived;
#include <model/Network/NetworkTypedefs.h>
        
        
        typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
        
        using NodeBaseType::NdofXnode;
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
        typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
        typedef std::deque<const LatticePlane*> LatticePlaneContainerType;
        typedef std::deque<LatticePlane> SpecialLatticePlaneContainerType;
        
        typedef LatticeVector<dim> LatticeVectorType;
        typedef LatticeDirection<dim> LatticeDirectionType;
        
        static bool use_velocityFilter;
        static double velocityReductionFactor;
        static double velocityIncreaseFactor;
        static double bndDistance;
        
        
    private:
        
        DislocationSharedObjects<dim> shared;
        
        LatticeVectorType L;
        
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        
        //! The std::vector containing the glidePlaneNormal(s) of the connected DislocationSegment(s)
        //		VectorOfNormalsType planenormals;
        LatticePlaneContainerType _confiningPlanes;
        SpecialLatticePlaneContainerType specialConfiningPlanes;
        
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        
        //! The previous velocity vector of *this DislocationNode
        VectorDofType vOld;
        
        double velocityReductionCoeff;
        
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        VectorDim boundaryNormal;
        
        //! The normal to the region boundary
        //        std::auto_ptr<LatticePlane> regionBndPlane;
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
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
                
                if(!temp.first) // DislocationNode not found inside mesh
                {
                    // Detect if the DislocationNode is sligtly outside the boundary
                    int faceID;
                    const double baryMin(temp.second->pos2bary(this->get_P()).minCoeff(&faceID));
                    const bool isApproxOnBoundary(std::fabs(baryMin)<1000.0*FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
                    if(!isApproxOnBoundary)
                    {
                        model::cout<<"DislocationNode "<<this->sID<<" @ "<<this->get_P().transpose()<<std::endl;
                        model::cout<<"Simplex "<<temp.second->xID<<std::endl;
                        model::cout<<"bary "<<temp.second->pos2bary(this->get_P())<<std::endl;
                        model::cout<<"face of barymin is "<<temp.second->child(faceID).xID<<std::endl;
                        model::cout<<"face of barymin is boundary Simplex? "<<temp.second->child(faceID).isBoundarySimplex()<<std::endl;
                        assert(0 && "DISLOCATION NODE CREATED OUTSIDE MESH.");
                    }
                }
            }
            
            return temp.second;
        }
        
        /**********************************************************************/
        void forceBoundaryNode(const EdgeRef<LinkType>& pL)
        {
            if(shared.use_boundary)
            {
                if(pL.E.is_boundarySegment() && !isBoundaryNode())
                {
                    const VectorDim nb=(pL.E.source->bndNormal()+pL.E.sink->bndNormal()).normalized();
                    const VectorDim np=pL.E.glidePlane.n.cartesian().normalized();
                    VectorDim outDir=nb-nb.dot(np)*np;
                    if(outDir.squaredNorm()>FLT_EPSILON)
                    {
                        //                        std::cout<<"DislocationNode "<<this->sID<<", outDir="<<outDir.transpose()<<std::endl;
                        //                        std::cout<<pL.E.source->get_P().transpose()<<std::endl;
                        //                        std::cout<<pL.E.sink->get_P().transpose()<<std::endl;
                        //                        std::cout<<this->get_P().transpose()<<std::endl;
                        
                        outDir.normalize();
                        LatticeVectorType dL(pL.E.glidePlane.n.snapToLattice(outDir));
                        assert(dL.squaredNorm()>0.0);
                        
                        LatticeVectorType L0=L;
                        if(!DislocationSharedObjects<dim>::mesh.searchWithGuess(L0.cartesian(),p_Simplex).first)
                        {
                            L0=LatticeVectorType(pL.E.glidePlane.snapToLattice(0.5*(pL.E.source->get_P()+pL.E.sink->get_P())));
                        }
                        assert(DislocationSharedObjects<dim>::mesh.searchWithGuess(L0.cartesian(),p_Simplex).first && "L0 outside mesh");
                        
                        LatticeLine line(L0,dL);
                        LineMeshIntersection lmi(line,L+dL,shared.mesh,p_Simplex);
                        if(lmi.search.first)
                        {
                            p_Simplex=lmi.search.second;
                            set(lmi.L);
                            //                            boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance);
                            make_bndNormal();
                            if(!isBoundaryNode())
                            {
                                //                                std::cout<<"DislocaitonNode "<<this->sID<<" not on mesh boundary"<<std::endl;
                                //                                std::cout<<"dL="<<dL.transpose()<<std::endl;
                            }
                            assert(isBoundaryNode());
                        }
                        else
                        {
                            assert(0 && "lmi failed.");
                        }
                    }
                    else
                    {
                        model::cout<<"DislocationNode "<<this->sID<<std::endl;
                        assert(0 && "outDir has zero norm");
                    }
                }
            }
        }
        
        /**********************************************************************/
        void moveToBoundary(const VectorDim& outDir,const std::pair<bool,
                            const Simplex<dim,dim>*>& temp,
                            const VectorDim& dX)
        {
            
            if(outDir.squaredNorm()<FLT_EPSILON)
            {
                model::cout<<"DislocationNode "<<this->sID<<" outDir has zero norm"<<std::endl;
                assert(0 && "outDir has zero norm");
            }
            
            LatticeVectorType dL(LatticeVectorType::Zero());
            switch (_confiningPlanes.size())
            {
                case 1:
                {
                    //std::cout<<"boundary motion case 1, DislocationNode "<<this->sID<<std::endl;
                    const VectorDim planeN(_confiningPlanes[0]->n.cartesian().normalized());
                    VectorDim dD=outDir-outDir.dot(planeN)*planeN;
                    const double dDnorm=dD.norm();
                    if(dDnorm>0.0)
                    {
                        dD.normalize();
                        dL=LatticeVectorType(_confiningPlanes[0]->n.snapToLattice(dD)); // a lattice vector on the plane pointing ouside mesh
                        assert(dL.squaredNorm()>0);
                    }
                    break;
                }
                    
                case 2:
                {
                    //std::cout<<"boundary motion case 2, DislocationNode "<<this->sID<<std::endl;
                    const PlanePlaneIntersection ppi(*_confiningPlanes[0],*_confiningPlanes[1]);
                    dL=LatticeVectorType(ppi.d.snapToDirection(10.0*outDir));
                    break;
                }
                    
            }
            
            if(dL.squaredNorm()>FLT_EPSILON) // a line direction could be found
            {
                const LatticeVectorType L0((this->get_P()+dX*temp.first).eval());
                const LatticeLine line(L0,dL);
                LineMeshIntersection lmi(line,L0+dL,shared.mesh,temp.second);
                assert(lmi.search.first);
                p_Simplex=lmi.search.second;
                set(lmi.L);
                //                boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance);
                make_bndNormal();
                if(meshLocation()!=onMeshBoundary)
                {
                    model::cout<<"DislocaitonNode "<<this->sID<<std::endl;
                    assert(0 && "NODE MUST BE ON MESH-BOUNDARY");
                }
            }
            
        }
        
        
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**********************************************************************/
        DislocationNode(const LatticeVectorType& Lin,
                        const Simplex<dim,dim>* guess=(const Simplex<dim,dim>*) NULL) :
        /* base constructor */ NodeBaseType(Lin.cartesian()),
        /* init list        */ L(Lin),
        /* init list        */ p_Simplex(get_includingSimplex(guess)),
        /* init list        */ velocity(VectorDofType::Zero()),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero())
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        {/*! Constructor from DOF
          */
        }
        
        /**********************************************************************/
        DislocationNode(const EdgeRef<LinkType>& pL,
                        const LatticeVectorType& Lin) :
        /* base constructor */ NodeBaseType(pL,Lin.cartesian()),
        /* init list        */ L(Lin),
        /* init list        */ p_Simplex(get_includingSimplex(pL.E.source->includingSimplex())),
        /* init list        */ velocity((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.E.source->velocityReduction()+pL.E.sink->velocityReduction())),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero())
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        {/*! Constructor from EdgeRef and DOF
          */
            //            std::cout<<"DislocationNode from ExpadingLink A "<<this->sID<<std::endl;
            forceBoundaryNode(pL);
        }
        
        /**********************************************************************/
        DislocationNode(const EdgeRef<LinkType>& pL,
                        const LatticeVectorType& Lin,
                        const VectorDofType& Vin) :
        /* base constructor */ NodeBaseType(pL,Lin.cartesian()),
        /* init list        */ L(Lin),
        /* init list        */ p_Simplex(get_includingSimplex(pL.E.source->includingSimplex())),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(pL.E.source->velocityReduction()+pL.E.sink->velocityReduction())),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero())
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        {
            //            std::cout<<"DislocationNode from ExpadingLink B "<<this->sID<<std::endl;
            forceBoundaryNode(pL);
        }
        
        /**********************************************************************/
        DislocationNode(const ContractingVertices<NodeType,LinkType>& cv,
                        const LatticeVectorType& Lin) :
        /* base constructor */ NodeBaseType(Lin.cartesian()),
        /* init list        */ L(Lin),
        /* init list        */ p_Simplex(get_includingSimplex(cv.v0.includingSimplex())),
        /* init list        */ velocity(0.5*(cv.v0.get_V()+cv.v1.get_V())),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(0.5*(cv.v0.velocityReduction()+cv.v1.velocityReduction())),
        /* init list        */ boundaryNormal(shared.use_boundary? SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance) : VectorDim::Zero())
        //        /* init list        */ regionBndNormal(VectorDim::Zero())
        {/*! Constructor from VertexContraction
          */
            
        }
        
        /**********************************************************************/
        void make_bndNormal()
        {
            boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance); // check if node is now on a boundary
        }
        
        /**********************************************************************/
        const LatticeVectorType& get_L() const
        {
            return L;
        }
        
        /**********************************************************************/
        void set(const LatticeVectorType& Lin)
        {
            L=Lin;
            NodeBaseType::set(L.cartesian());
        }
        
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
        
        const LatticePlaneContainerType& confiningPlanes() const
        {
            return _confiningPlanes;
        }
        
        /**********************************************************************/
        void make_confiningPlanes()
        {
            
            //std::cout<<"DislocationNode "<<this->sID<<"make_confiningPlanes"<<std::endl;
            
            //! 1- Clear and re-builds the std::vector planenormals
            _confiningPlanes.clear();
            specialConfiningPlanes.clear();
            
            if (!this->is_balanced() && meshLocation()!=onMeshBoundary)
            {
                for(const auto& planeBase : CrystalOrientation<dim>::planeNormals())
                {
                    specialConfiningPlanes.emplace_back(L,planeBase);
                }
                
            }
            
            // add to _confiningPlanes the planes of the attached segments
            for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
            {
                if (std::get<2>(neighborIter->second))
                {
                    LinkType* pL(std::get<1>(neighborIter->second));
                    _confiningPlanes.push_back(&(pL->glidePlane));
                    _confiningPlanes.push_back(&(pL->sessilePlane));
                }
            }
            
            // add to _confiningPlanes the special planes of this node
            for(const auto& plane : specialConfiningPlanes)
            {
                _confiningPlanes.push_back(&plane);
            }
            
            //std::cout<<"a"<<std::endl;
            
            GramSchmidt::makeUnique(_confiningPlanes);
            
            //std::cout<<"b"<<std::endl;
            
            //! 2- Compute projectionMatrix
            make_projectionMatrix();
            
            //std::cout<<"done"<<std::endl;
            
        }
        
        /**********************************************************************/
        void removeFromNeighborhood(LinkType* const pL)
        {/*!@param[in] pL A pointer to the DislocationSegment being disconnected
          * from this node.
          *
          *  Overwrites NetworkNode::removeFromNeighborhood in order to modify
          *  planeNormals and tangent after the DislocationSegment is disconnected
          */
            NodeBaseType::removeFromNeighborhood(pL);
            //            make_planeNormals();
            make_confiningPlanes();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
            NodeBaseType::make_T();
        }
        
        /**********************************************************************/
        void addToNeighborhood(LinkType* const pL)
        {/*!@param[in] pL A pointer to the DislocationSegment being connected
          * to this node.
          *
          *  Overwrites NetworkNode::addToNeighborhood in order to modify
          *  planeNormals and tangent after the DislocationSegment is disconnected
          */
            NodeBaseType::addToNeighborhood(pL);
            //            make_planeNormals();
            make_confiningPlanes();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
            NodeBaseType::make_T();
        }
        
        /**********************************************************************/
        VectorOfNormalsType constraintNormals() const
        {
            VectorOfNormalsType temp;
            
            for(const auto& plane : _confiningPlanes)
            {
                temp.push_back(plane->n.cartesian().normalized());
            }
            
            
            
            if (meshLocation()==insideMesh)
            { // DislocationNode is inside mesh
                if (!this->is_balanced())
                { // DislocationNode is not balanced, fix it
                    temp.push_back((VectorDim()<<1.0,0.0,0.0).finished());
                    temp.push_back((VectorDim()<<0.0,1.0,0.0).finished());
                    temp.push_back((VectorDim()<<0.0,0.0,1.0).finished());
                }
                
                // Constrain simple nodes to move normal to tangent
                if(this->is_simple() && temp.size()==1) // simple node on one glide plane
                {
                    Eigen::Matrix<double,dim,1> T(this->get_T());
                    double normT(T.norm());
                    if (normT>FLT_EPSILON)
                    {
                        temp.push_back(T/normT);
                    }
                }
                
            }
            else if (meshLocation()==onMeshBoundary)
            { // DislocationNode is on mesh boundary, constrain by boundaryNormal
                temp.push_back(boundaryNormal);
            }
            else
            {
                //std::cout<<"DislocationNode "<<this->sID<< " at "<<this->get_P().transpose()<<" is outside mesh."<<std::endl;
                assert(0 && "DISLOCATION NODE FOUND OUTSIDE MESH."); //RE-ENABLE THIS
            }
            //showTemp(temp);
            GramSchmidt::orthoNormalize(temp);
            assert(temp.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
            return temp;
        }
        
        
        
        /**********************************************************************/
        void make_projectionMatrix()
        {
            
            
                // Add normal to glide and sessile planes
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
        
        /**********************************************************************/
        void set_V(const VectorDofType& vNew)
        {
            vOld=velocity; // store current value of velocity before updating
            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
           

            if(use_velocityFilter)
            {
                const double filterThreshold=0.05*velocity.norm()*vOld.norm();
 
                if(velocity.dot(vOld)<-filterThreshold)
                {
                    velocityReductionCoeff*=velocityReductionFactor;
                }
                else if(velocity.dot(vOld)>filterThreshold)
                {
                    velocityReductionCoeff*=velocityIncreaseFactor;
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
        bool invertedMotion() const
        {/*! The nodal velocity vector
          */
            return velocity.template segment<dim>(0).dot( vOld.template segment<dim>(0) ) < 0.0;
        }
        
        /**********************************************************************/
        void move(const double & dt)
        {
            const VectorDim P_old(this->get_P());
            
            VectorDim dX=velocity.template segment<dim>(0)*dt;
            
            //Limit dX for boundaryNodes bec
            const double dXnorm(dX.norm());
            if((isBoundaryNode() || isConnectedToBoundaryNodes()) && dXnorm>10.0)
            {
                dX*=10.0/dXnorm;
            }
            
            switch (_confiningPlanes.size())
            {
                case 1:
                    dX=_confiningPlanes[0]->n.snapToLattice(dX);
                    break;
                    
                case 2:
                {
                    PlanePlaneIntersection ppi(*_confiningPlanes[0],*_confiningPlanes[1]);
                    LatticeLine line(ppi.P,ppi.d);
                    dX=line.d.snapToDirection(dX);
                    break;
                }
                    
                default:
                    dX.setZero();
                    break;
            }
            
            //            if(isPureBoundaryNode())
            //            {
            //                dX.setZero();
            //            }
            
            //			if (dX.squaredNorm()>0.0 && (meshLocation()!=onMeshBoundary || shared.use_bvp==0)) // move a node only if |v|!=0 and if not on mesh boundary
            if (dX.squaredNorm()>0.0) // move a node only if |v|!=0
            {
                if(shared.use_boundary) // using confining mesh
                {
                    // See if the new position is inside mesh
                    std::set<const Simplex<dim,dim>*> path;
                    const std::pair<bool,const Simplex<dim,dim>*> temp(DislocationSharedObjects<dim>::mesh.searchWithGuess(true,this->get_P()+dX,p_Simplex,path));
                    //p_Simplex=temp.second;
                    
                    
                    if(isBoundaryNode()) // node already a boundary node, it must remain on boundary
                    {
                        const VectorDim bndNrml=SimplexBndNormal::get_boundaryNormal(this->get_P()+dX,*temp.second,bndDistance); // boundary normal at new position
                        if(temp.first && bndNrml.squaredNorm()>0.0) // new position is on boundary
                        {
                            p_Simplex=temp.second;
                            set(L+LatticeVectorType(dX));
                            //                            boundaryNormal=bndNrml;
                            make_bndNormal();
                            assert(meshLocation()==onMeshBoundary);
                        }
                        else // new position is not on boundary
                        {
                            moveToBoundary(boundaryNormal,temp,dX);
                        }
                    }
                    else // not a boundary node
                    {
                        if(!temp.first) // node moved outside or already on boundary
                        {
                            
                            
                            VectorDim outDir=boundaryNormal;
                            if(outDir.squaredNorm()==0.0)
                            {// node is exiting for the first time, we need a tentative boundary normal to identify the "outside direction"
                                for(const auto& simplex : path) // loop over the pathof simplices  connecting P to P+dX
                                {
                                    outDir=SimplexBndNormal::get_boundaryNormal(this->get_P()+dX,*simplex,dX.norm()+FLT_EPSILON);
                                    
                                    if(outDir.squaredNorm()>0.0)
                                    {
                                        break;
                                    }
                                    
                                }
                                
                            }
                            assert(outDir.squaredNorm()>0 && "COULD NOT DETERMINE OUTDIR");
                            moveToBoundary(outDir,temp,dX);
                        }
                        else // node is internal and remains internal
                        {
                            //                            if(temp.second->region->regionID!=p_Simplex->region->regionID)
                            //                            {// node is crossing regions
                            //                                assert(0 && "RE-ENABLE THIS");
                            //                                //                            //// ////std::cout<<"DislocationNode "<<this->sID<<" crossing region. Path size="<<path.size()<<std::endl;
                            //                                //
                            //                                //                            int faceID=-1;
                            //                                ////                            const Simplex<dim,dim>* regBndSimplex=(const Simplex<dim,dim>*) NULL;
                            //                                //                            const Simplex<dim,dim>* regBndSimplex(NULL);
                            //                                //
                            //                                //                            Eigen::Matrix<double,dim+1,1> faceInt(Eigen::Matrix<double,dim+1,1>::Zero());
                            //                                //
                            //                                //                            for(const auto& simplex : path) // loop over the pathof simplices  connecting P to P+dX
                            //                                //                            {
                            //                                //                                const Eigen::Matrix<double,dim+1,1> baryOld(simplex->pos2bary(this->get_P()));
                            //                                //                                const Eigen::Matrix<double,dim+1,1> baryNew(simplex->pos2bary(this->get_P()+dX));
                            //                                //                                for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                            //                                //                                {
                            //                                //                                    if(simplex->child(f).isRegionBoundarySimplex())
                            //                                //                                    {
                            //                                //                                        faceInt=simplex->faceLineIntersection(baryOld,baryNew,f);
                            //                                //                                        //                                    //////std::cout<<"DislocationNode "<<this->sID<<", baryMin="<<faceInt.minCoeff()<<std::endl;
                            //                                //
                            //                                //                                        if(faceInt.minCoeff()>=-FLT_EPSILON) // faceInt belongs to triangle
                            //                                //                                        {
                            //                                //                                            regBndSimplex=simplex; // current simplex is the region boundary simplex wanted
                            //                                //                                            faceID=f; // intersection face is f
                            //                                //                                            break;
                            //                                //                                        }
                            //                                //                                    }
                            //                                //                                }
                            //                                //
                            //                                //                                if(faceID>=0)
                            //                                //                                {
                            //                                //                                    break;
                            //                                //                                }
                            //                                //                            }
                            //                                //
                            //                                //                            assert(faceID>=0 && "FACE INTERSECTION NOT FOUND");
                            //                                //                            this->set(regBndSimplex->bary2pos(faceInt)); // move node to intersction with region boundary
                            //                                //                            regionBndNormal=regBndSimplex->nda.col(faceID).normalized();
                            //                            }
                            //                            else // node not crossing regions
                            //                            {
                            p_Simplex=temp.second;
                            set(L+LatticeVectorType(dX));
                            //boundaryNormal=SimplexBndNormal::get_boundaryNormal(this->get_P(),*p_Simplex,bndDistance); // check if node is now on a boundary
                            make_bndNormal();
                            //                                //                                L+=LatticeVectorType(dX);
                            //                                //                                this->set(this->get_P()+dX); // move node
                            //                                //
                            //                                //                            PROBLEM HERE, BOUNDSARY NODES ARE ASSIGNED ZERO NORMAL
                            //
                            //                            }
                            
                        }
                    }
                    
                    make_projectionMatrix();
                    
                }
                else // move node freely
                {
                    set(L+LatticeVectorType(dX));
                    //                    L+=LatticeVectorType(dX);
                    //                    this->set(this->get_nodeDof()+dX);
                }
            }
            
            // Store actual velocity
            if(dt>0.0)
            {
                velocity=(this->get_P()-P_old)/dt;
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
        bool isPureBoundaryNode() const
        {
            return isBoundaryNode() && isConnectedToBoundaryNodes();
        }
        
        /**********************************************************************/
        bool isConnectedToBoundaryNodes() const
        {
            bool temp(!this->is_isolated());
            for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
            {
                if (std::get<2>(neighborIter->second)) // not self
                {
                    temp*=std::get<0>(neighborIter->second)->isBoundaryNode();
                }
            }
            
            return temp;
        }
        
        /**********************************************************************/
        std::pair<std::deque<std::pair<size_t,size_t> >,std::deque<std::pair<size_t,size_t> >> edgeDecomposition() const
        {
            
            std::deque<std::pair<size_t,size_t> > firstLinkDeq;
            std::deque<std::pair<size_t,size_t> > secondLinkDeq;
            
            if(this->Neighborhood.size()>2 && meshLocation()==insideMesh)
            {
                std::deque<std::pair<int,int> > linkDeq;
                
                
                const size_t neighSize=this->Neighborhood.size()-1;
                
                Eigen::MatrixXd temp(neighSize,dim);
                
                int k=0;
                for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
                {
                    if (std::get<2>(neighborIter->second)) // not self
                    {
                        LinkType* pL(std::get<1>(neighborIter->second));
                        temp.row(k)=pL->Burgers;
                        linkDeq.emplace_back(pL->source->sID,pL->sink->sID);
                        k++;
                    }
                }
                
                //                if(temp.rows()>7)
                //                {
                //                    model::cout<<"Dislocation Node "<<this->sID<<std::endl;
                //                }
                auto vecVec=EdgePermutations::edgeStats(temp);
                if(vecVec.size()==2)
                {
                    if((vecVec[0].array()*vecVec[1].array()).matrix().squaredNorm()<FLT_EPSILON) // decompisition is unique
                    {
                        //                       //std::cout<<"DislocationNode "<<this->sID<<std::endl;
                        //                       //std::cout<<"Matrix of Burgers= "<<std::endl<<temp<<std::endl;
                        //                       //std::cout<<"null configurations= "<<std::endl;
                        //                        for(auto v : vecVec)
                        //                        {
                        //                           //std::cout<<v<<std::endl;
                        //                        }
                        for(size_t n=0;n<neighSize;++n)
                        {
                            if(vecVec[0](n)!=0)
                            {
                                firstLinkDeq.push_back(linkDeq[n]);
                            }
                            else
                            {
                                secondLinkDeq.push_back(linkDeq[n]);
                            }
                        }
                        
                    }
                }
            }
            
            return make_pair(firstLinkDeq,secondLinkDeq);
        }
        
        /******************************************************************************/
        void neighborsAt(const LatticeVectorType& L0, std::set<size_t>& temp) const
        {/*!\param[in] P0 position to be serached
          * \param[out]temp set of IDs of neighbors of this which are located at P0 (possibly including *this)
          * \param[in] tol tolerance used to detect position overlap
          */
            for (typename Derived::NeighborContainerType::const_iterator nIiter =this->neighborhood().begin();
                 /*                                                   */ nIiter!=this->neighborhood().end();
                 /*                                                   */ ++nIiter)
            { // loop over neighborhood
                if((std::get<0>(nIiter->second)->get_L()-L0).squaredNorm()==0)
                { // a neighbor of I exists at P0
                    temp.insert(std::get<0>(nIiter->second)->sID);
                }
            }
        }
        
        /**********************************************************************/
        const double& velocityReduction() const
        {
            return velocityReductionCoeff;
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const NodeType& ds)
        {
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.get_P().transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.get_T().transpose()<<"\t"
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
    double DislocationNode<_dim,corder,InterpolationType,QuadratureRule>::velocityIncreaseFactor=1.1;
    
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    double DislocationNode<_dim,corder,InterpolationType,QuadratureRule>::bndDistance=2.0;
    
    
} // close namespace
#endif


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
//                                const PlanePlaneIntersection ppi(*_confiningPlanes[0],*_confiningPlanes[1]);
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


//        /**********************************************************************/
//        void make_projectionMatrix()
//        {
//
//            // Add normal to glide and sessile planes
//            Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
//            VectorOfNormalsType  CN;
//            for(const auto& plane : _confiningPlanes)
//            {
//                CN.push_back(plane->n.cartesian().normalized());
//            }
//
//            // Add normal to mesh boundary
//            if(meshLocation()==onMeshBoundary)
//            {
//                const Eigen::Matrix<double,dim+1,1> bary(p_Simplex->pos2bary(this->get_P()));
//                for (int i=0;i<dim+1;++i)
//                {
//                    if(std::fabs(bary(i))<FLT_EPSILON && p_Simplex->child(i).isBoundarySimplex())
//                    {
//                        CN.push_back(p_Simplex->nda.col(i));
//                    }
//                }
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
//        }
