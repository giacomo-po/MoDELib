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

namespace model
{
	
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//	class DislocationNode : public SplineNodeBase<DislocationNode<_dim,corder,InterpolationType,qOrder,QuadratureRule>,
//	/*                                         */ _dim,corder,InterpolationType>
    class DislocationNode : public SplineNode<DislocationNode<_dim,corder,InterpolationType,qOrder,QuadratureRule>,
    /*                                         */ _dim,corder,InterpolationType>

    {
		
	public:
		
        // make dim available outside class
        constexpr static int dim=_dim;
        
        
        // define Derived to use NetworkTypedefs.h
		typedef DislocationNode       <dim,corder,InterpolationType,qOrder,QuadratureRule> Derived;
#include <model/Network/NetworkTypedefs.h>
		
		
		typedef SplineNode<NodeType,dim,corder,InterpolationType> NodeBaseType;
		
		using NodeBaseType::NdofXnode;
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;
		typedef Eigen::Matrix<double,dim,dim> MatrixDim;
		typedef Eigen::Matrix<double,NdofXnode,1> VectorDofType;
		typedef std::vector<VectorDim,Eigen::aligned_allocator<VectorDim> > VectorOfNormalsType;
		
        static bool use_velocityFilter;
        static double velocityReductionFactor;

        
	private:
		
		DislocationSharedObjects<LinkType> shared;
		
        //! A pointer to the Simplex containing *this
        const Simplex<dim,dim>* p_Simplex;
        
		//! The std::vector containing the glidePlaneNormal(s) of the connected DislocationSegment(s)
		VectorOfNormalsType planenormals;
		
        //! The current velocity vector of *this DislocationNode
        VectorDofType velocity;
        
        //! The previous velocity vector of *this DislocationNode
		VectorDofType vOld;
        
        double velocityReductionCoeff;
        
        //! The normal unit vector of the boundary on which *this DislocationNode is moving on
        VectorDim boundaryNormal;
        
        
        /**********************************************************************/
        const Simplex<dim,dim>* get_includingSimplex(const Simplex<dim,dim>* const guess) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if (DislocationSharedObjects<LinkType>::use_boundary)
            {
                if (guess==NULL)
                {
                    temp=DislocationSharedObjects<LinkType>::mesh.search(this->get_P());
                }
                else
                {
                    temp=DislocationSharedObjects<LinkType>::mesh.searchWithGuess(this->get_P(),guess);
                }
                
                if(!temp.first) // DislocationNode not found inside mesh
                {
                    // Detect if the DislocationNode is sligtly outside the boundary
                    int faceID;
                    const double baryMin(temp.second->pos2bary(this->get_P()).minCoeff(&faceID));
                    const bool isApproxOnBoundary(std::fabs(baryMin)<FLT_EPSILON && temp.second->child(faceID).isBoundarySimplex());
                    assert(isApproxOnBoundary && "DISLOCATION NODE CREATED OUTSIDE MESH.");
                }
            }
            
            return temp.second;
        }
        

		
        
		
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
		/* Constructor ********************************************************/
		DislocationNode(const VectorDofType& Qin,const Simplex<dim,dim>* guess=(const Simplex<dim,dim>*) NULL) :
        /* base constructor */ NodeBaseType(Qin),
        /* init list        */ p_Simplex(get_includingSimplex(guess)),
        /* init list        */ velocity(VectorDofType::Zero()),
		/* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ boundaryNormal(get_boundaryNormal())
        {/*! Constructor from DOF
          */
		}
		
		/* Constructor ********************************************************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const double& u) :
        /* base constructor */ NodeBaseType(pL,u),
        /* init list        */ p_Simplex(get_includingSimplex(pL.E.source->includingSimplex())),
        /* init list        */ velocity((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ boundaryNormal(get_boundaryNormal())
        {/*! Constructor from ExpandingEdge and parameter along link
          */
		}
		
		/* Constructor ********************************************************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const VectorDofType& Qin) :
        /* base constructor */ NodeBaseType(pL,Qin),
        /* init list        */ p_Simplex(get_includingSimplex(pL.E.source->includingSimplex())),
        /* init list        */ velocity((pL.E.source->velocity+pL.E.sink->velocity)*0.5), // TO DO: this should be calculated using shape functions from source and sink nodes of the link
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ boundaryNormal(get_boundaryNormal())
        {/*! Constructor from ExpandingEdge and DOF
          */
		}
		
		/* Constructor from Link and position along link **********************/
		DislocationNode(const ExpandingEdge<LinkType>& pL, const VectorDofType& Qin, const VectorDofType& Vin)
        /* base constructor */ : NodeBaseType(pL,Qin),
        /* init list        */ p_Simplex(get_includingSimplex(pL.E.source->includingSimplex())),
        /* init list        */ velocity(Vin),
        /* init list        */ vOld(velocity),
        /* init list        */ velocityReductionCoeff(1.0),
        /* init list        */ boundaryNormal(get_boundaryNormal())
        {
		}
        
        
        /**********************************************************************/
		void make_planeNormals()
        {
            //! 1- Clear and re-builds the std::vector planenormals
			planenormals.clear();
			for (typename NeighborContainerType::const_iterator neighborIter=this->Neighborhood.begin();neighborIter!=this->Neighborhood.end();++neighborIter)
            {
				if (std::get<2>(neighborIter->second))
                {
					LinkType* pL(std::get<1>(neighborIter->second));
					if (std::find(planenormals.begin(),planenormals.end(), pL->glidePlaneNormal )==planenormals.end() &&
						std::find(planenormals.begin(),planenormals.end(),-pL->glidePlaneNormal )==planenormals.end()   )
                    {
						planenormals.push_back(pL->glidePlaneNormal );
					}
					if(pL->sessilePlaneNormal.norm()>FLT_EPSILON)
                    {
						planenormals.push_back(pL->sessilePlaneNormal);
					}
					
				}
			}
			
			//! 2- Compute projectionMatrix
			make_projectionMatrix();
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
            make_planeNormals();
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
            make_planeNormals();
            DislocationEnergyRules<dim>::template findEdgeConfiguration<NodeType>(*this);
            NodeBaseType::make_T();
        }
        
        /**********************************************************************/
        VectorDim get_boundaryNormal() const
        {
            VectorDim temp(VectorDim::Zero());
            if (shared.use_boundary)
            {
                const Eigen::Matrix<double,dim+1,1> bary(p_Simplex->pos2bary(this->get_P()));
                int faceID;
                const double baryMin(bary.minCoeff(&faceID)); // this also writes faceID
                if (std::fabs(baryMin)<FLT_EPSILON && p_Simplex->child(faceID).isBoundarySimplex())
                {
                    temp=p_Simplex->nda.col(faceID).normalized();
                }
            }
            return temp;
        }
        
        /**********************************************************************/
		const VectorOfNormalsType& planeNormals() const
        {
            return planenormals;
        }
		
        /**********************************************************************/
		VectorOfNormalsType constraintNormals() const
        {
			GramSchmidt<dim> GS;
			if (meshLocation()==insideMesh)
            { // DislocationNode is inside mesh
				if (this->is_balanced())
                { // DislocationNode is balanced
					GS.orthoNormalize(planenormals); //  constrain by planenormals
				}
				else
                { // DislocationNode is unbalanced and is inside mesh: fix
					GS.push_back((VectorDim()<<1.0,0.0,0.0).finished());
					GS.push_back((VectorDim()<<0.0,1.0,0.0).finished());
					GS.push_back((VectorDim()<<0.0,0.0,1.0).finished());
				}
			}
			else if (meshLocation()==onMeshBoundary)
            { // DislocationNode is on mesh boundary:
				VectorOfNormalsType PN(planenormals); // constrain by planenormals
				PN.push_back(boundaryNormal); // constrain by boundaryNormal
				GS.orthoNormalize(PN);
			}
			else
            {
                std::cout<<"DislocationNode "<<this->sID<< " at "<<this->get_P().transpose()<<" is outside mesh."<<std::endl;
                assert(0 && "DISLOCATION NODE FOUND OUTSIDE MESH."); //RE-ENABLE THIS
			}
			assert(GS.size()>=1 && "GLIDING NODE MUST HAVE AT LEAST ONE CONSTRAINT.");
			return GS;
		}
		
        /**********************************************************************/
		bool is_removable() const
        {
			return (this->is_simple() && constraintNormals().size()==1);
		}
		
        /**********************************************************************/
		void make_projectionMatrix()
        {
			Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
			VectorOfNormalsType  CN = planenormals;
            
            if(meshLocation()==onMeshBoundary)
            {
                const Eigen::Matrix<double,dim+1,1> bary(p_Simplex->pos2bary(this->get_P()));
                for (int i=0;i<dim+1;++i)
                {
                    if(std::fabs(bary(i))<FLT_EPSILON && p_Simplex->child(i).isBoundarySimplex())
                    {
                        CN.push_back(p_Simplex->nda.col(i));
                    }
                }
            }
            
			//CN.push_back(boundaryNormal);
			GramSchmidt<dim> GS(CN);
			this->prjM.setIdentity();
			for (size_t k=0;k<GS.size();++k)
            {
				this->prjM*=( I-GS[k]*GS[k].transpose() );
			}
		}
		
        /**********************************************************************/
		void set_V(const VectorDofType& vNew)
        {
            vOld=velocity; // store current value of velocity before updating
            velocity=this->prjM*vNew; // kill numerical errors from the iterative solver
            
            if(use_velocityFilter)
            {
                if(velocity.dot(vOld)<0.0)
                {
                    velocityReductionCoeff*=velocityReductionFactor;
                }
                else
                {
                    velocityReductionCoeff/=velocityReductionFactor;
                }
                if(velocityReductionCoeff>1.0)
                {
                    velocityReductionCoeff=1.0;
                }
                if(velocityReductionCoeff<0.01)
                {
                    velocityReductionCoeff=0.01;
                }
                velocity*=velocityReductionCoeff;
            }
            
		}
        
//        /**********************************************************************/
//        void implicitStep()
//        {
//            velocity= (velocity+vOld)*0.5; // this brings back
//        }
		
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
//		void move(const double & dt, const double & dt_old)
        void move(const double & dt)
        {
			
			VectorDim dX=velocity.template segment<dim>(0)*dt;
            
            //			if (dX.squaredNorm()>0.0 && (meshLocation()!=onMeshBoundary || shared.use_bvp==0)) // move a node only if |v|!=0 and if not on mesh boundary
            if (dX.squaredNorm()>0.0) // move a node only if |v|!=0
            {
				if(shared.use_boundary) // using confining mesh
                {
                    // See if the new position is inside mesh
                    const std::pair<bool,const Simplex<dim,dim>*> temp(DislocationSharedObjects<LinkType>::mesh.searchWithGuess(this->get_P()+dX,p_Simplex));
                    p_Simplex=temp.second;
                    
                    if(temp.first) // new position is inside mesh
                    {
                        this->set(this->get_P()+dX); // move node
                        boundaryNormal=get_boundaryNormal(); // check if node is now on a boundary
                    }
                    else // new position is outside mesh
                    {
                        
                        int faceID;
                        const Eigen::Matrix<double,dim+1,1> baryNew(p_Simplex->pos2bary(this->get_P()+dX));
                        const double baryMin(baryNew.minCoeff(&faceID)); // this also finds faceID
                        assert(p_Simplex->child(faceID).isBoundarySimplex() && "FACE MUST BE A BOUNDARY FACE");
                        
                        if(baryMin>-FLT_EPSILON) // DislocationNode is sligtly outside the boundary
                        {
                            this->set(this->get_P()+dX); // move node
                            boundaryNormal=get_boundaryNormal(); // check if node is now on a boundary
                        }
                        else // Node is completely outside the boundary. We bring it back to the boundary.
                        {
                            const Eigen::Matrix<double,dim+1,1> baryOld(p_Simplex->pos2bary(this->get_P()));
                            const Eigen::Matrix<double,dim+1,1> faceInt(p_Simplex->faceLineIntersection(baryOld,baryNew,faceID));
                            dX=p_Simplex->bary2pos(faceInt)-this->get_P();
                            this->set(this->get_P()+dX);
                            boundaryNormal=p_Simplex->nda.col(faceID).normalized();
                            velocity=dX/dt; // correct stored velocity
                        }
                        
                        assert(meshLocation()==onMeshBoundary && "NODE MUST NOW BE ON MESH BOUNDARY.");
                    }
                    
                    make_projectionMatrix();
                    
				}
				else // move node freely
                {
					this->set(this->get_nodeDof()+dX);
				}
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
		
        /**********************************************************************/
        std::deque<std::pair<size_t,size_t> > edgeDecomposition() const
        {
            
            std::deque<std::pair<size_t,size_t> > firstLinkDeq;
            
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
                
                if(temp.rows()>7)
                {
                    std::cout<<"Dislocation Node "<<this->sID<<std::endl;
                }
                auto vecVec=EdgePermutations::edgeStats(temp);
                if(vecVec.size()==2)
                {
                    if((vecVec[0].array()*vecVec[1].array()).matrix().squaredNorm()<FLT_EPSILON) // decompisition is unique
                    {
//                        std::cout<<"DislocationNode "<<this->sID<<std::endl;
//                        std::cout<<"Matrix of Burgers= "<<std::endl<<temp<<std::endl;
//                        std::cout<<"null configurations= "<<std::endl;
//                        for(auto v : vecVec)
//                        {
//                            std::cout<<v<<std::endl;
//                        }
                        for(int n=0;n<neighSize;++n)
                        {
                            if(vecVec[0](n)!=0)
                            {
                                firstLinkDeq.push_back(linkDeq[n]);
                            }
                        }
                        
                    }
                }
            }
            
            return firstLinkDeq;
        }
        
    };
    
    
    // static data
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    bool DislocationNode<_dim,corder,InterpolationType,qOrder,QuadratureRule>::use_velocityFilter=true;

    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
    double DislocationNode<_dim,corder,InterpolationType,qOrder,QuadratureRule>::velocityReductionFactor=0.75;

    
} // close namespace
#endif
