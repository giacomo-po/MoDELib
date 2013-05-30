/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNODE_H_
#define model_DISLOCATIONNODE_H_

#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <model/Network/Operations/EdgeExpansion.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationConsts.h>
#include <model/Geometry/Splines/SplineNodeBase.h>
#include <model/BVP/SearchData.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/Math/GramSchmidt.h>
#include <model/DislocationDynamics/DislocationEnergyRules.h>
#include <model/BVP/Tetrahedron.h>

namespace model {
	
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNode : public SplineNodeBase<DislocationNode<_dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>,
	/*                                         */ _dim,corder,InterpolationType>{
		
	public:
		
        // make dim available outside class
        enum{dim=_dim};
        
        // define Derived to use NetworkTypedefs.h
		typedef DislocationNode       <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule> Derived;
#include <model/Network/NetworkTypedefs.h>
		
		
		typedef SplineNodeBase<NodeType,dim,corder,InterpolationType> NodeBaseType;
		
		using NodeBaseType::NdofXnode;
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
		typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
		
		
		
	private:
		
		DislocationSharedObjects<LinkType> shared;
		
        //! The mesh ID containing this
		int  currentMeshID;
		
		//! The std::vector containing the glidePlaneNormal(s) of the connected DislocationSegment(s)
		VectorOfNormalsType planenormals;
		
		
		/* initMeshLocation ***************************************************/
		void initMeshLocation()
        {
			if (shared.boundary_type){
				model::SearchData<dim> SD(this->get_P());
				shared.domain.findIncludingTet(SD);
				
				nodeMeshLocation = SD.nodeMeshLocation;
				currentMeshID = SD.newMeshID;
				
				if (SD.nodeMeshLocation == onMeshBoundary){
					boundaryNormal=SD.outwardFaceNormal;
					triIndex=SD.triIndex;
				}
				if (SD.nodeMeshLocation == outsideMesh) {
                    std::cout<< "NODE "<<this->sID<<" IS OUTSIDE DOMAIN AT " << this->get_P().transpose() << std::endl;
                    assert(0 && "DISLOCATION NODE CREATED OUTSIDE DOMAIN.");
                }
			}
			
		}
		
		
		VectorDofType velocity;
		VectorDofType vOld;
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		
		
		int nodeMeshLocation; // 1=inside
		VectorDim   boundaryNormal;
        
        
		MatrixDim bvpStress;
		
		unsigned int triIndex;
        
		/* Constructor ********************************************************/
		DislocationNode(const VectorDofType& Qin) :
        /* base constructor */ NodeBaseType::SplineNodeBase(Qin),
        /* init list        */ velocity(VectorDofType::Zero()),
		/* init list        */ vOld(VectorDofType::Zero()),
		/* init list        */ nodeMeshLocation(insideMesh),
		/* init list        */ boundaryNormal(VectorDim::Zero()),
		/* init list        */ bvpStress(MatrixDim::Zero())
        {/*! Constructor from DOF
          */
			initMeshLocation();
		}
		
		/* Constructor ********************************************************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const double& u) :
        /* base constructor */ NodeBaseType::SplineNodeBase(pL,u),
        /* init list        */ velocity((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
		/* init list        */ vOld((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
		/* init list        */ nodeMeshLocation(insideMesh),
		/* init list        */ boundaryNormal(VectorDim::Zero()),
		/* init list        */ bvpStress(MatrixDim::Zero()) // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        {/*! Constructor from ExpandingEdge and parameter along link
          */
			initMeshLocation();
		}
		
		/* Constructor ********************************************************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const VectorDofType& Qin) :
        /* base constructor */ NodeBaseType::SplineNodeBase(pL,Qin),
        /* init list        */ velocity((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
		/* init list        */ vOld((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
		/* init list        */ nodeMeshLocation(insideMesh),
		/* init list        */ boundaryNormal(VectorDim::Zero()),
		/* init list        */ bvpStress(MatrixDim::Zero()) // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        {/*! Constructor from ExpandingEdge and DOF
          */
			initMeshLocation();
		}
		
		/* Constructor from Link and position along link **********************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const VectorDofType& Qin, const VectorDofType& Vin)
        /* base constructor */ : NodeBaseType::SplineNodeBase(pL,Qin),
        /* init list        */ velocity(Vin),
		/* init list        */ vOld(velocity), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
		/* init list        */ nodeMeshLocation(insideMesh),
		/* init list        */ boundaryNormal(VectorDim::Zero()),
		/* init list        */ bvpStress(MatrixDim::Zero()) // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        {    
            initMeshLocation();
		}
        

        /* make_planeNormals **************************************************/
		void make_planeNormals()
        {
            //! 1- Clear and re-builds the std::vector planenormals
			planenormals.clear();
			for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter){
				if (std::get<2>(neighborIter->second)){
					LinkType* pL(std::get<1>(neighborIter->second));
					if (std::find(planenormals.begin(),planenormals.end(), pL->glidePlaneNormal )==planenormals.end() &&
						std::find(planenormals.begin(),planenormals.end(),-pL->glidePlaneNormal )==planenormals.end()   ){
						planenormals.push_back(pL->glidePlaneNormal );
					}
					if(pL->sessilePlaneNormal.norm()>FLT_EPSILON){
						planenormals.push_back(pL->sessilePlaneNormal);
					}
					
				}
			}
			
			//! 2- Compute projectionMatrix
			make_projectionMatrix();
        }

        /* make_planeNormals **************************************************/
        void removeFromNeighborhood(LinkType* const pL)
        {/*! @param[in] pL A pointer to the DislocationSegment being disconnected
          *  Overwrites NetworkNode::removeFromNeighborhood in order to modify
          *  planeNormals and tangent after the DislocationSegment is disconnected
          */
            NodeBaseType::removeFromNeighborhood(pL);
            make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
            NodeBaseType::make_T();
        }
		
		/* meshID *************************************************************/
		const int& meshID() const
        {
            return currentMeshID;
        }
        
		/* planeNormals *******************************************************/
		const VectorOfNormalsType& planeNormals() const
        {
            return planenormals;
        }
		
		/* constraintNormals **************************************************/
		VectorOfNormalsType constraintNormals() const
        {
			GramSchmidt<dim> GS;
			if (nodeMeshLocation==insideMesh){ // DislocationNode is inside mesh
				if (this->is_balanced()){ // DislocationNode is balanced
					GS.orthoNormalize(planenormals); //  constrain by planenormals
				}
				else{ // DislocationNode is unbalanced and is inside mesh: fix
					GS.push_back((VectorDim()<<1.0,0.0,0.0).finished());
					GS.push_back((VectorDim()<<0.0,1.0,0.0).finished());
					GS.push_back((VectorDim()<<0.0,0.0,1.0).finished());
				}
			}
			else if (nodeMeshLocation==onMeshBoundary){ // DislocationNode is on mesh boundary:
				VectorOfNormalsType PN(planenormals); // constrain by planenormals
				PN.push_back(boundaryNormal); // constrain by boundaryNormal
				GS.orthoNormalize(PN);
			}
			else{
                std::cout<<"DislocationNode "<<this->sID<< " at "<<this->get_P().transpose()<<" is outside mesh."<<std::endl;
                assert(0 && "DISLOCATION NODE FOUND OUTSIDE MESH."); //RE-ENABLE THIS
			}
			assert(GS.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
			return GS;
		}
		
		///////////////////////////////
		bool is_removable() const
        {
			return (this->is_simple() && constraintNormals().size()==1);
		}
		
		
        /* make_projectionMatrix **********************************************/
		void make_projectionMatrix()
        {
			Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
			VectorOfNormalsType  CN = planenormals;
			CN.push_back(boundaryNormal);
			GramSchmidt<dim> GS(CN);
			this->prjM.setIdentity();
			for (size_t k=0;k<GS.size();++k){
				this->prjM*=( I-GS[k]*GS[k].transpose() );
			}
		}
		
		/* updateBvpStress ****************************************************/
		void updateBvpStress()
        {
			bvpStress=shared.domain.tetContainer[currentMeshID].getStress(); // stress is constant in the element because of Linear Shape Functions
		}
		
		
		/***************************************/
		VectorDim deformedPosition() const {
			VectorDim temp(this->get_P());
			Eigen::Matrix<double,4,1> barycoord(shared.domain.tetContainer[currentMeshID].getBarycentric(this->get_P()));
			for(int k=0;k<4;++k){
				temp+=shared.domain.tetContainer[currentMeshID].eleNodes[0]->u*barycoord(k);
			}
			return temp;
		}
		
		
		/***************************************/
		void set_V(const VectorDofType& vNew)
        {
            vOld=velocity; // store current value of velocity before updating
            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
		}
        
        /***************************************/
        void implicitStep()
        {
            velocity= (velocity+vOld)*0.5; // this brings back
        }
		
		/***************************************/
		const VectorDofType& get_V() const
        {/*! The nodal velocity vector
          */
			return velocity;
		}
        
		
		
		/***************************************/
		void move(const double & dt, const double & dt_old)
        {
			
			VectorDim dX=velocity.template segment<dim>(0)*dt - vOld.template segment<dim>(0)*dt_old;
			if (dX.squaredNorm()>0.0 && (nodeMeshLocation!=onMeshBoundary || shared.use_bvp==0)) // move a node only if |v|>0 and if not on mesh boundary
            { 
				if(shared.boundary_type)
                {
					
					VectorDim newP, dir;
					if (nodeMeshLocation){ // inside=1 or boundary=2
						newP=this->get_P()+dX;
						dir=dX;
					}
					else{ // outside=0
						newP=this->get_P();
						dir=dX;
					}
					
    				
					model::SearchData<dim> SD(newP,dir,currentMeshID,nodeMeshLocation,triIndex,boundaryNormal);
                    //					model::SearchData<dim> SD(newP,dir,currentMeshID,nodeMeshLocation,triIndex); // OLD
					
					shared.domain.SearchMovingNode(SD);
					
					
					currentMeshID=SD.newMeshID;
					nodeMeshLocation=SD.nodeMeshLocation;
                    
                    
					switch (SD.nodeMeshLocation)
                    {
						case onMeshBoundary:
							boundaryNormal=SD.outwardFaceNormal;
							assert(dir.cross(SD.projectedP-this->get_P()).norm()<FLT_EPSILON && "CORRECTION NOT ALIGNED WITH DIR");
							velocity=(SD.projectedP-this->get_P())/dt;
							this->set(SD.projectedP);
							triIndex=SD.triIndex;
							break;
						case insideMesh:
							boundaryNormal=VectorDim::Zero();
							this->set(this->get_nodeDof()+velocity*dt);
							break;
						default:
							assert(0);
							break;
					}
					make_projectionMatrix();
				}
				else
                {
					this->set(this->get_nodeDof()+velocity*dt);
				}
			}
			
		}
		
        /* operator<< *********************************************************/
        template <class T>
		friend T& operator << (T& os, const NodeType& ds)
        {
			os  << ds.sID<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.get_P().transpose()<<" "
			/**/<< std::setprecision(15)<<std::scientific<<ds.get_T().transpose()<<" "
            /**/<< ds.pSN()->sID<<" ";
//			if (ds.shared.use_bvp)
//            { //output in deformed configuration
//					os << std::setprecision(15)<<std::scientific<<ds.deformedPosition().transpose()<<" ";
//				}
//				else{
//					os<< VectorDim::Zero().transpose();
//				}
				//os << "\n";
			return os;
        }
		
	}; // close DislocationNode
	
} // close namespace model
#endif
   