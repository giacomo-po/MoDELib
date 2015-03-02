/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

namespace model {
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // SPECIALIZATION CORDR=0, INTRPOLATION=Hermite
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,0,Hermite> : public NetworkNode<Derived>{
        
        
    public:
        
        
        enum {corder=0, NdofXnode=dim*(corder+1)};
        
#include "model/Geometry/Splines/SplineNodeBase_corder0.h"
        
        
        
        typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
        typedef Eigen::Matrix<int, NdofXnode, 1>    VectorDofIDType;
        
        
        //////////////////////////////////////////////////
        // Constructor
        SplineNode(VectorDofType Q_in) :  P(Q_in.segment<dim>(0)){
        }
        
        //////////////////////////////////////////////////
        // set
        void set(VectorDofType Q_in){
            
            //! 1- Sets P=Q_in
            this->P=Q_in;
            
            //! 2- Updates all SplineSegments connected to this
            typedef void (LinkType::*link_member_function_pointer_type)(void);
            link_member_function_pointer_type Lmfp;
            Lmfp=&LinkType::update;
            this->linkTransmit(Lmfp,1);
        }
        
        
        //////////////////////////////////////////////////
        // get_nodeDof
        VectorDofType get_nodeDof() const {
            return this->P;
        }
        
        //////////////////////////////////////////////////
        // dofID
        Eigen::VectorXi dofID() const {
            VectorDofIDType VectorDofID;
            
            for (int k=0;k<NdofXnode;++k){
                VectorDofID(k)=NdofXnode*this->snID()+k;
            }
            
            return VectorDofID;
        }
        
        
        
    };
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // SPECIALIZATION CORDER=1, INTRPOLATION=Hermite
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,1,Hermite> : public NetworkNode<Derived>{
        
        /*!	\brief  Implementation of Hermite SplineNode for the case of C1 continuity. Degrees of freedom are
         *	the nodal position vector and nodal tangent vector.
         *  (tangent continuity at each node).
         *
         *	Spline Nodes of type SplineNode<Derived,dim,1,Hermite> guarantee that spline
         *	segments in the spline network meet at each node with continuity of both position
         *	and parametric tangent. Each node is defined by a position and a tangent vector which
         *  are passed together in VectorDofType using the set(VectorDofType Q_in) member function.
         
         *
         */
        
    public:
        
        enum {corder=1, NdofXnode=dim*(corder+1)};
        
#include "model/Geometry/Splines/SplineNodeBase_corder1.h"
        
        typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
        typedef Eigen::Matrix<double, NdofXnode, NdofXnode> MatrixDofType;
        typedef Eigen::Matrix<int, NdofXnode, 1> VectorDofIDType;
        
        
        
        //////////////////////////////////////////////////
        // Constructor
        SplineNode(VectorDofType Q_in) :  P(Q_in.template segment<dim>(0)),	 T(Q_in.template segment<dim>(dim)){
            prjM.setIdentity();
        }
        
        
        //////////////////////////////////////////////////
        // set
        void set(VectorDofType Q_in){
            /*!	1- Update of the nodal degrees of freedom (position and tangent).
             *	\param[Q_in] is a vector of length dim*(corder+1) = dim*2.
             *  The first dim values are the position vector while the second block of dim values
             *	are the tangent vector. Tangent vector must have non-vanishing norm.
             */
            this->P=Q_in.template segment<dim>(0);
            this->T=Q_in.template segment<dim>(dim);
            assert(T.squaredNorm()>0.0);
            
            //! 2- Updates all SplineSegments connected to this
            typedef void (LinkType::*link_member_function_pointer_type)(void);
            link_member_function_pointer_type Lmfp;
            Lmfp=&LinkType::update;
            this->linkTransmit(Lmfp,1);
            
            
        }
        
        //////////////////////////////////////////////////
        // set
        VectorDofType get_nodeDof() const {
            return (VectorDofType() << this->P, this->T).finished();
        }
        
        //////////////////////////////////////////////////
        // dofID
        VectorDofIDType dofID() const {
            VectorDofIDType VectorDofID;
            
            for (int k=0;k<NdofXnode;++k){
                VectorDofID(k)=NdofXnode*this->snID()+k;
            }
            
            return VectorDofID;
        }
        
        MatrixDofType W2H()const {
            return MatrixDofType::Identity();
        }
        
        
        
        //		//////////////////////////////////////////////////
        //		// GnDof
        //		size_t GnDof(){
        //			return NdofXnode*NetworkNode<SubNetworkType, Derived, SegmentType>::Naddresses();
        //		}
        
        //		//////////////////////////////////////////////////
        //		// get_NdofXnode
        //		size_t get_NdofXnode(){
        //			return NdofXnode;
        //		}
        
    };
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // SPECIALIZATION CORDER=2, INTRPOLATION=Hermite
    template <typename Derived, short unsigned int dim>
    class SplineNode<Derived, dim,2,Hermite> : public NetworkNode<Derived>{
        
        
        
    public:
        
        enum {corder=2, NdofXnode=dim*(corder+1)};
        
#include "model/Geometry/Splines/SplineNodeBase_corder2.h"
        
        
        
        
        
        typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
        typedef Eigen::Matrix<int,    NdofXnode, 1> VectorDofIDType;
        
        
        
        
        //////////////////////////////////////////////////
        // Constructor 
        SplineNode(VectorDofType Q_in) : P(Q_in.segment<dim>(0)),
        /*														 */	 T(Q_in.segment<dim>(dim)),
        /*														 */	 K(Q_in.segment<dim>(2*dim)){
            prjM.setIdentity();
        }
        
        
        //////////////////////////////////////////////////
        // set
        void set(VectorDofType Q_in){
            this->P=Q_in.segment<dim>(0);
            this->T=Q_in.segment<dim>(dim);
            assert(T.squaredNorm()>0.0);
            this->K=Q_in.segment<dim>(2*dim);
            
            //! 2- Updates all SplineSegments connected to this
            typedef void (LinkType::*link_member_function_pointer_type)(void); 
            link_member_function_pointer_type Lmfp;
            Lmfp=&LinkType::update;
            this->linkTransmit(Lmfp,1);
            
        }
        
        //////////////////////////////////////////////////
        // set
        VectorDofType get_nodeDof() const {
            return (VectorDofType() << this->P, this->T, this->K).finished();
        }
        
        //////////////////////////////////////////////////
        // dofID
        Eigen::VectorXi dofID() const {	
            VectorDofIDType VectorDofID;
            
            for (int k=0;k<NdofXnode;++k){
                VectorDofID(k)=NdofXnode*this->snID()+k;
            }
            
            return VectorDofID;
        }
        
        
        
        
    };
    
    
}


