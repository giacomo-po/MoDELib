/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


template <typename Derived, short unsigned int dim>
class SplineNodeBase<Derived, dim,1,CatmullRom> :	public NetworkNode<Derived>{ 
	
	
	
public:
    
	Eigen::VectorXi edgeConfiguration;
    
	
	enum  {corder=1};
	enum  {NdofXnode=dim};
	
	typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
	
#include "model/Geometry/Splines/SplineNodeBase_corder1.h"
	typedef SplineNodeBase<Derived, dim,corder,CatmullRom> NodeType;
	
	
public:
	
	/* Constructor from position **************************************************************/
	SplineNodeBase(const VectorDim& P_in) : 	
	/* init list */ P(P_in),
	/* init list */ T(VectorDim::Zero()),	
	/* init list */ prjM(MatrixDim::Identity()){		
		set(get_nodeDof()); // trigger calculation of tangents		
	}
	
	/* Constructor from EdgeExpansion and parameter along edge ********************************/
	//SplineNodeBase(const EdgeExpansion<LinkType>& pL, const double& u) : 
	SplineNodeBase(const ExpandingEdge<LinkType>& pL, const double& u) : 
	/* init list */ NetworkNode<Derived>::NetworkNode(pL),	
	/* init list */ P(pL.E.get_r(u)),
	/* init list */ prjM(MatrixDim::Identity()){
		set(get_nodeDof()); // trigger calculation of tangents
	}
	
	/* Constructor from EdgeExpansion and position along edge *********************************/
    //	SplineNodeBase(const EdgeExpansion<LinkType>& pL, const VectorDim& P_in) : 
	SplineNodeBase(const ExpandingEdge<LinkType>& pL, const VectorDim& P_in) : 
	/* init list */ NetworkNode<Derived>::NetworkNode(pL),
	/* init list */ P(P_in),
	/* init list */ prjM(MatrixDim::Identity()){
		set(get_nodeDof()); // trigger calculation of tangents
	}
	
	
	
	
	//////////////////////////////////////////////////
	// set
	void set(const VectorDim& P_in){
		/*! Because of the CatmullRom rule for parametric tangents, changing the position of this 
		 *  CatmullRom node must also change its parametric tangent and the
		 *  parametric tangents of its first-neighboring nodes. This in turn affects the shape
		 *  of all SplineSegments attached to the first-neighboring nodes. The update strategy is as follows:
		 */
		
		//! 1- Sets P=P_in
		this->P=P_in;
		
		//! 2- Transmits 'make_T' on the neighbor nodes of level (corder+1=2), that is this node and its first-neighbors
		typedef void (Derived::*node_member_function_pointer_type)(void); 
		node_member_function_pointer_type Nmfp;
		//	Nmfp=&Derived::update;
		Nmfp=&Derived::make_T;
		//		this->nodeTransmit(Nmfp,corder+2);
        //		this->nodeTransmit(Nmfp,corder+1); // transmit to 1st neighbors (corder+1=2)
		this->depthFirstNodeExecute(Nmfp,corder+1); // transmit to 1st neighbors (corder+1=2)
	}
	
	
	
	
	//////////////////////////////////////////////////
	// topologyChangeActions
	void topologyChangeActions(){
		/*! Because of the CatmullRom rule for parametric tangents, changing the connectivity of this 
		 *  CatmullRom node must also change its parametric tangent. This in turn affects the shape
		 *  of all SplineSegments attached to it. The update strategy is as follows:
		 */
		
		
        //		energyRule();
		this->p_derived()->findEdgeConfiguration();
		make_T();
		
		
		//		//! 1- Recalculates tangent  
		//		typedef void (Derived::*node_member_function_pointer_type)(void); 
		//		node_member_function_pointer_type Nmfp;
		//		//	Nmfp=&Derived::update;
		//		Nmfp=&Derived::make_T;
		//		this->nodeTransmit(Nmfp,corder+2);
		
		// MAKE_CRNEIGHNORS AND MAKE_T SHOULD BE SEPARATED SO THAT HERE ONE TRANSMITS MAKE_CRNEIGHNORS AND MAKE_T,
		// WHILE IN SET ONLY MAKE_T IS TRANSMITTED
		
		//		//! 2- Transmits 'update' on the neighbor links of level (corder=1), that is the first neighbor links
		//		typedef void (LinkType::*link_member_function_pointer_type)(void); 
		//		link_member_function_pointer_type Lmfp;
		//		Lmfp=&LinkType::update;
		//		this->linkTransmit(Lmfp,corder+2);
	}
	
	//////////////////////////////////////////////////
	// make_T
	void make_T(){
		
		//! Calls make_CR2H()
		make_CR2H();
		
		//! computes T=(CR2H*VectorDof).segment(dim,dim)
        this->T=(CR2H*VectorDof).template segment<dim>(dim);
        
 //       print182();
        
	}
	
	
	//////////////////////////////////////////////////
	// get_nodeDof
	VectorDim get_nodeDof() const {
		return this->P;
	}	
	
	//////////////////////////////////////////////////
	// dofID
	Eigen::Matrix<int,dim,1> node_dofID() const {	
		return ( (Eigen::Array<int,dim,1>() << 0, 1, 2).finished()+this->snID()*dim).matrix();
	}
	
	//////////////////////////////////////////////////
	// get_dof
	Eigen::VectorXd get_dof() const {
		return VectorDof;
	}	
	
	//////////////////////////////////////////////////
	// dofID
	Eigen::VectorXi dofID() const {	
		
		Eigen::VectorXi VectorDofID;
		VectorDofID.resize(Ndof);
		
		size_t k=0;
		for (constNeighborIteratorType neighborIter=this->neighborhood().begin(); neighborIter != this->neighborhood().end(); ++neighborIter) {
			VectorDofID.segment<dim>(k*dim)	=boost::tuples::get<0>(neighborIter->second)->node_dofID();
			++k;
		}
		
		return VectorDofID;
	}
	
	//////////////////////////////////////////////////
	// W2H
	const Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> & W2H() const {	
		
		//			make_CR2H();
		return CR2H;
	}
	
	//////////////////////////////////////////////////
	// get_CRneighbors
    //	NeighborContainerType get_CRneighbors() const {
    //		return CRneighbors;
    //	}
	
	
	
//    void print182(){
//        if (this->sID==182){
//            std::cout<<"Node "<<this->sID<<" is_balanced="<<this->is_balanced()<<", is_isolated="<<this->is_isolated()<<std::endl; 
//            std::cout<<"\n"<<CR2H<<std::endl;
//            //         std::cout<<this->prjM<<std::endl;
//            //          std::cout<<"Node "<<this->sID<< " making T="<<this->T.transpose()<<std::endl;
//        }
//    }
    
	
	
private:
	
	
	//	size_t Nneighbors;
	size_t Ndof;
	
	Eigen::VectorXd VectorDof;
	//	Eigen::VectorXi VectorDofID;
	
	Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> CR2H;
	
	
    
    /* findEdgeConfiguration (possibly overwritten by Derived) ************************/
	void findEdgeConfiguration(){
        // 
		assert(0 && "NEED TO RE-IMPLEMENT ORIGINAL CR-RULE");
	}
    
    
    
	
	/* make_CRneighbors ************************************/
	void make_CRneighbors(){
		//size_t Nneighbors=this->neighborhood().size();
		Ndof=dim*this->neighborhood().size();
		VectorDof.resize(Ndof);			
		size_t k=0;
		for (constNeighborIteratorType neighborIter=this->neighborhood().begin(); neighborIter != this->neighborhood().end(); ++neighborIter) {
			VectorDof.segment<dim>(k*dim)	=boost::tuples::get<0>(neighborIter->second)->get_P();
			++k;
		}
	}
	
	
	
	
	
	//////////////////////////////////////////////////
	// make_CR2H
	void make_CR2H(){
		
		//! 1- Call make_CRneighbors() 
		make_CRneighbors();
		
		//! 2- Resize CR2H to dim*(corder+1) x Ndof
		CR2H.setZero(dim*(corder+1),Ndof);
		
		
		
		if (this->is_isolated()){
            //            std::cout<<"Node "<<this->sID<<" is isolated"<<std::endl;
			CR2H.template block<dim,dim>(0,0).setIdentity();
			CR2H.template block<dim,dim>(dim,0).setZero();		
		}
		else{
			if (this->is_balanced()){
                //                        std::cout<<"Node "<<this->sID<<" is not isolated and balanced "<<std::endl;
				make_CR2H_central();
			}		
			else{
                //             std::cout<<"Node "<<this->sID<<" is not isolated and not balanced "<<std::endl;
                
				int k=0;				
				for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
					const short int dir(boost::tuples::get<2>(neighborIter->second));
					switch ( dir ) {
						case  0:	// self
							CR2H.template block<dim,dim>(0,k*dim).setIdentity();
							//							CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( (int(this->outOrder())-int(this->inOrder() ))/CPLT - CPLDPinv + CPLARinv)/(CRneighbors.size()-2);						
							break;
					}
					++k;
				}
			}
		}
        
	}
	
    
	
    
	
	
	//////////////////////////////////////////////////
	// make_CR2H_central
	void make_CR2H_central(){
		
		//	double CPL;
		//	double CPLDPinv=0.0;
		//	double CPLARinv=0.0;
		double CPLT=0.0;
		double sjT(0.0);
		double sjOverGjT(0.0);
		
		int sgnID(0);
		for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
			switch ( boost::tuples::get<2>(neighborIter->second) ) {
				case  0:	// self
					break;
				default:	// neighbor
					double CPL=boost::tuples::get<1>(neighborIter->second)->chordParametricLength();	// chord parametric length
					CPLT+=CPL;
					sjT+=edgeConfiguration(sgnID);
					sjOverGjT+=edgeConfiguration(sgnID)/CPL;
					sgnID++;
					break;
			}
		}
		
		double CPLTinv=1.0/CPLT;
		//Ndof=dim*this->neighborhood().size();
		//		short int dir;			
		int k=0;
		sgnID=0;
		
		for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
			int dir=boost::tuples::get<2>(neighborIter->second);
			switch ( dir ) {
				case  0:	// self
					CR2H.template block<dim,dim>(0,k*dim).setIdentity();
					CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( sjT/CPLT -  sjOverGjT)/(this->neighborhood().size()-2);						
					break;
				default:	// departing or arriving
					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/boost::tuples::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(this->neighborhood().size()-2);
					sgnID++;
					break;
			}
			++k;
		}
		
	}
	
    
    
    
}; 

