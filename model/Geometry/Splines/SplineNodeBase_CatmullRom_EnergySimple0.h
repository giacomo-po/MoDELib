/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 Benjamin Ramirez <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


//#include <mmdl/CompileTimeMath/Pow.h>

//#include <algorithm> // or find_if


//bool isPositive (int i) {
//	return ((i%2)==1);
//}

/*******************************************************************************/
/* class template specialization: SplineNodeBase<Derived,dim,1,CatmullRom>     */
/*******************************************************************************/



template <unsigned int N>
struct EdgeConfigs{
	
	static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci(){
		Eigen::Matrix<int,Pow<2,N-1>::value,N> temp;
		temp.template block<Pow<2,N-2>::value,N-1>(0,1)=EdgeConfigs<N-1>::get_Ci();
		temp.template block<Pow<2,N-2>::value,N-1>(Pow<2,N-2>::value,1)=-EdgeConfigs<N-1>::get_Ci();
		temp.col(0).setOnes();
		return temp;
	}
	
//	static const Eigen::Matrix<int,Pow<2,N-1>::value,N> Ci;
};

//! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
//template<unsigned int N>
//const Eigen::Matrix<int,Pow<2,N-1>::value,N> EdgeConfigs<N>::Ci=EdgeConfigs<N>::get_Ci();


template <>
struct EdgeConfigs<1>{
	enum {N=1};
	
	static Eigen::Matrix<int,Pow<2,N-1>::value,N> get_Ci(){
//		Eigen::Matrix<int,Pow<2,N-1>::value,N> temp;
//		temp<<1;
//		return temp;
		return Eigen::Matrix<int,Pow<2,N-1>::value,N>::Ones();
	}
	
//	static const Eigen::Matrix<int,Pow<2,N-1>::value,N> Ci;
};

//! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
//template<>
//const Eigen::Matrix<int,1,1> EdgeConfigs<1>::Ci=EdgeConfigs<1>::get_Ci();




template <typename Derived, short unsigned int dim>
class SplineNodeBase<Derived, dim,1,CatmullRom> :	public NetworkNode<Derived>{ 
	
	Eigen::VectorXi edgeConfiguration;
	
	
public:
	
	enum  {corder=1};
	enum  {NdofXnode=dim};
	
	typedef Eigen::Matrix<double, NdofXnode, 1> VectorDofType;
	
#include "mmdl/Geometry/Splines/SplineNodeBase_corder1.h"
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
//		depthFirstExecute
		
		//		//! 3- Transmits 'update' on the neighbor links of level (corder+1=2), that is the first neighbor links of first neighbor nodes
		//		typedef void (LinkType::*link_member_function_pointer_type)(void); 
		//		link_member_function_pointer_type Lmfp;
		//		Lmfp=&LinkType::update; //! CatmullRom: WRONG!!! UPDATE IS ONLY IN THE DISLOCATION LAYER!!!! 
		//		this->linkTransmit(Lmfp,corder+2);
	}
	
	
	
	
	//////////////////////////////////////////////////
	// topologyChangeActions
	void topologyChangeActions(){
		/*! Because of the CatmullRom rule for parametric tangents, changing the connectivity of this 
		 *  CatmullRom node must also change its parametric tangent. This in turn affects the shape
		 *  of all SplineSegments attached to it. The update strategy is as follows:
		 */
		
		
		energyRule();
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
		for (constNeighborIteratorType neighborIter=CRneighbors.begin(); neighborIter != CRneighbors.end(); ++neighborIter) {
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
	NeighborContainerType get_CRneighbors() const {
		return CRneighbors;
	}
	
	
	
	
	
private:
	
	NeighborContainerType CRneighbors;
	
	//	size_t Nneighbors;
	size_t Ndof;
	
	Eigen::VectorXd VectorDof;
	//	Eigen::VectorXi VectorDofID;
	
	Eigen::Matrix<double, dim*(corder+1), Eigen::Dynamic> CR2H;
	
	//	Eigen::MatrixXd BB;
	//	std::vector<typename NetworkNode<Derived>::FlowType>  BsV;
	
	
	//	Eigen::MatrixXd CR2K0;
	//	Eigen::MatrixXd K0source;
	//	Eigen::MatrixXd K0sink;
	
	
	/* make_CRneighbors ************************************/
	void make_CRneighbors(){
		//size_t tempKey;
		// Derived* tempNode;
		//		typedef typename NetworkNode<Derived>::LinkType LinkType;
		//		LinkType* tempLink=NULL;
		//		NeighborContainerType NeighborOFNeighbor;
		//! 1- The Catmul-Rom Neighbors starts as the this->neighborhood()
		CRneighbors = this->neighborhood();
		//		if (this->is_source() || this->is_sink()) {
		//			// 2- Modify the Catmul-Rom Neighbors in case of source or sink node
		//			for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
		//				if (boost::tuples::get<0>(NeighborIter->second)->is_central()){
		//					NeighborOFNeighbor=boost::tuples::get<0>(NeighborIter->second)->neighborhood();
		//					for (NeighborIteratorType NeighborIter2=NeighborOFNeighbor.begin();NeighborIter2!=NeighborOFNeighbor.end();++NeighborIter2){
		//						size_t tempKey=boost::tuples::get<0>(NeighborIter2->second)->sID;
		//						Derived* tempNode=boost::tuples::get<0>(NeighborIter2->second);
		//						CRneighbors.insert(std::make_pair(tempKey,boost::tuples::make_tuple(tempNode,tempLink,0)));
		//					}
		//				}
		//				
		//			}
		//		}
		size_t Nneighbors=CRneighbors.size();
		Ndof=dim*Nneighbors;
		VectorDof.resize(Ndof);			
		size_t k=0;
		for (NeighborIteratorType neighborIter=CRneighbors.begin(); neighborIter != CRneighbors.end(); ++neighborIter) {
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
			CR2H.template block<dim,dim>(0,0).setIdentity();
			CR2H.template block<dim,dim>(dim,0).setZero();		
		}
		else{
			if (this->is_balanced()){
				make_CR2H_central();
			}		
			else{
				//				if (this->is_source() || this->is_sink()){
				int k=0;
				
				for (NeighborIteratorType neighborIter=CRneighbors.begin();neighborIter!=CRneighbors.end();++neighborIter){
					short int dir=boost::tuples::get<2>(neighborIter->second);
					switch ( dir ) {
						case  0:	// self
							CR2H.template block<dim,dim>(0,k*dim).setIdentity();
							//							CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( (int(this->outOrder())-int(this->inOrder() ))/CPLT - CPLDPinv + CPLARinv)/(CRneighbors.size()-2);						
							break;
							//						default:	// departing or arriving
							//							CR2H.template block<dim,dim>(dim,k*dim)= dir*this->prjM*( 1.0/boost::tuples::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(CRneighbors.size()-2);
							//							break;
					}
					++k;
				}
			}
		}
		
		
		//		if (this->is_source()) {
		//			make_CR2H_source();
		//		}
		//		else if (this->is_sink()) {
		//			make_CR2H_sink();
		//		}
		//		else if (this->is_central()) {
		//			make_CR2H_central();
		//		}
		//		else {
		//			//assert(0);
		//			CR2H.template block<dim,dim>(0,0).setIdentity();
		//			CR2H.template block<dim,dim>(dim,0).setZero();
		//			
		//		}
		
	}
	
	
	
	//////////////////////////////////////////////////////////////
	void energyRule(){
		
		if (this->is_isolated()){
			edgeConfiguration.setZero(0);
		}
		else{
			std::vector<typename NetworkNode<Derived>::FlowType>  BsV;
			//		BsV.clear();
			
			
			
			// Collect Burgers as if all edges were exiting the node		
			//		for (constNeighborIteratorType neighborIter=this->outNeighborhood().begin();neighborIter!=this->outNeighborhood().end();++neighborIter){
			//			BsV.push_back(+1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
			//		}	
			//		for (constNeighborIteratorType neighborIter=this-> inNeighborhood().begin();neighborIter!=this-> inNeighborhood().end();++neighborIter){
			//			BsV.push_back(-1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
			//		}	
			
			
			for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
				int dir=boost::tuples::get<2>(neighborIter->second);
				switch ( dir ) {
					case   1:
						BsV.push_back(+1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
						break;
					case  -1:
						BsV.push_back(-1.0*(boost::tuples::get<1>(neighborIter->second)->flow));
						break;
					default:	// self
						break;
				}
			}
			
			
			
			
			
			Eigen::MatrixXd BB(Eigen::MatrixXd::Zero(BsV.size(),BsV.size()));
			
			
			//BB.setZero(BsV.size(),BsV.size());		
			for (unsigned int i=0;i<BsV.size();++i){
				for (unsigned int j=0;j<BsV.size();++j){
					BB(i,j)=BsV[i].dot(BsV[j]);
				}
			}
			
//			enum {N1=1, N2=2, N3=3, N4=4, N5=5, N6=6, N7=7, N8=8};
			
			Eigen::MatrixXi Ci;
			switch (BsV.size()) {
				case 1:
					Ci=EdgeConfigs<1>::get_Ci();
					break;
				case 2:
//					Ci.resize(Pow<2,N2-1>::value,N2);
					Ci=EdgeConfigs<2>::get_Ci();
					break;
				case 3:
					Ci=EdgeConfigs<3>::get_Ci();
					break;
				case 4:
					Ci=EdgeConfigs<4>::get_Ci();
					break;
				case 5:
					Ci=EdgeConfigs<5>::get_Ci();
					break;
				case 6:
					Ci=EdgeConfigs<6>::get_Ci();
					break;
				case 7:
					Ci=EdgeConfigs<7>::get_Ci();
					break;
				case 8:
					Ci=EdgeConfigs<8>::get_Ci();
					break;
				case 9:
					Ci=EdgeConfigs<9>::get_Ci();
					break;
				case 10:
					Ci=EdgeConfigs<10>::get_Ci();
					break;
				case 11:
					Ci=EdgeConfigs<11>::get_Ci();
					break;
				case 12:
					Ci=EdgeConfigs<12>::get_Ci();
					break;
					
					
				default:
					std::cout<<"Node "<<this->sID<<": #edges="<<BsV.size()<<std::endl;
					for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
						int dir=boost::tuples::get<2>(neighborIter->second);
						switch ( dir ) {
							case   1:
								std::cout<<this->sID<<"->"<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"\n";
								//<<std::end;
								break;
							case  -1:
								std::cout<<(boost::tuples::get<0>(neighborIter->second)->sID)<<"->"<<this->sID<<"\n";			
								break;
							default:	// self
								break;
						}
					}
					
					assert(0 && "NEED TO IMPLEMENT MULTIPLE JUNCTIONS");
					break;
			}
			
			int sigCi=1;
			for (constNeighborIteratorType neighborIter=this->neighborhood().begin();neighborIter!=this->neighborhood().end();++neighborIter){
				int dir=boost::tuples::get<2>(neighborIter->second);
				if (dir!=0){
					sigCi=dir;
					break;
				}
			}
			Ci*=sigCi;
			
			// Compute and sort energy levels
			std::multimap<double,Eigen::VectorXi> ELev;	// MAP DOES NOT ACCEPT EQUAL VALUES SO multimap IS USED TO STORE  DEGENERATE STATES
			for (int k=0;k<Ci.rows();++k){
				ELev.insert(std::make_pair(Ci.row(k).cast<double>()*BB*Ci.row(k).cast<double>().transpose(),Ci.row(k)));
			}
			
			
//			bool energyLevelsAreValid(false);
			//		Eigen::VectorXi edgeConfiguration(Ci.row(0)*0);
			edgeConfiguration.setZero(Ci.cols());
			
			if (this->is_balanced()){
				assert(std::fabs(BB.determinant())<FLT_EPSILON && "BB MATRIX MUST BE SINGULAR FOR BALANCED NODES");
				assert(ELev.size()>1 && "MORE THAN ONE ENERGY LEVELS MUST BE FOUND FOR A BALANCED NODE");	
				
				//std::map<double,Eigen::VectorXi>::const_iterator firstPositiveE=std::find_if(ELev.begin(),ELev.size.end(),isPositive)
//				assert(ELev.begin()->first>=0.0 &&		  "FIRST ENERGY LEVEL MUST BE 0.");
				assert(std::fabs(ELev.begin()->first)<FLT_EPSILON && "FIRST ENERGY LEVEL MUST BE 0.");

				std::multimap<double,Eigen::VectorXi>::const_iterator firstNonZero(ELev.lower_bound(FLT_EPSILON)); // the first element that compares >=FLT_EPSILON
				assert(firstNonZero!=ELev.end() && "AT LEAST ONE POSITIVE ENERGY LEVEL MUST EXIST");
				std::multimap<double,Eigen::VectorXi>::const_iterator nextNonZero(firstNonZero);
				++nextNonZero;
				if (nextNonZero==ELev.end()){	// MINIMUM ENERGY LEVEL IS LAST ONE SO IT IS UNIQUE
					edgeConfiguration=firstNonZero->second;
				}
				else{
					if(std::fabs(nextNonZero->first-firstNonZero->first)>=FLT_EPSILON){ // MINIMUM ENERGY LEVEL IS UNIQUE
						edgeConfiguration=firstNonZero->second;
					}
				}
				

				
				
			}
			
			
			
			
			//		if (energyLevelsAreValid){
			// set sourceTfactor/sinkTfactor to 0
			unsigned int kk(0);
			NeighborContainerType tempNeigh=this->neighborhood();
			for (NeighborIteratorType neighborIter=tempNeigh.begin();neighborIter!=tempNeigh.end();++neighborIter){
				int dir=boost::tuples::get<2>(neighborIter->second);
				switch ( dir ) {
					case   1:	// out
						//std::cout<<"Writing sourceTfactor="<<edgeConfiguration(kk)<<" (kk="<<kk<<")"<<
						//					" in edge ("<<boost::tuples::get<1>(neighborIter->second)->source->sID<<"->"<<
						//					boost::tuples::get<1>(neighborIter->second)->sink->sID<<")"<<std::endl;
						boost::tuples::get<1>(neighborIter->second)->sourceTfactor=edgeConfiguration(kk);
						//					std::cout<<"sourceTfactor is now: "<<boost::tuples::get<1>(neighborIter->second)->sourceTfactor<<std::endl;
						kk++;
						break;
					case  -1:	// in
						//					std::cout<<"Writing sinkTfactor="<<edgeConfiguration(kk)<<" (kk="<<kk<<")"<<
						//					" in edge ("<<boost::tuples::get<1>(neighborIter->second)->source->sID<<"->"<<
						//					boost::tuples::get<1>(neighborIter->second)->sink->sID<<")"<<std::endl;
						boost::tuples::get<1>(neighborIter->second)->  sinkTfactor=edgeConfiguration(kk);
						//					std::cout<<"sinkTfactor is now: "<<boost::tuples::get<1>(neighborIter->second)->sinkTfactor<<std::endl;
						kk++;
						break;
					default:	// self
						break;
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
		for (NeighborIteratorType neighborIter=CRneighbors.begin();neighborIter!=CRneighbors.end();++neighborIter){
			switch ( boost::tuples::get<2>(neighborIter->second) ) {
				case  0:	// departing
					break;
				default:	// self
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
		
		for (NeighborIteratorType neighborIter=CRneighbors.begin();neighborIter!=CRneighbors.end();++neighborIter){
			int dir=boost::tuples::get<2>(neighborIter->second);
			switch ( dir ) {
				case  0:	// self
					CR2H.template block<dim,dim>(0,k*dim).setIdentity();
					CR2H.template block<dim,dim>(dim,k*dim)=this->prjM*( sjT/CPLT -  sjOverGjT)/(CRneighbors.size()-2);						
					break;
				default:	// departing 
					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/boost::tuples::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(CRneighbors.size()-2);
					sgnID++;
					break;
					//				case -1:	// arriving 
					//					CR2H.template block<dim,dim>(dim,k*dim)= edgeConfiguration(sgnID)*this->prjM*( 1.0/boost::tuples::get<1>(neighborIter->second)->chordParametricLength()-CPLTinv)/(CRneighbors.size()-2);
					//					sgnID++;
					//					break;
			}
			++k;
		}
		
	}
	
}; 



//		for (NeighborIteratorType neighborIter=tempNeigh.begin();neighborIter!=tempNeigh.end();++neighborIter){
//			int dir=boost::tuples::get<2>(neighborIter->second);
//			switch ( dir ) {
//				case   0:	
//					break;
//				default:	// self
//					std::cout<<" Edge ("<<boost::tuples::get<1>(neighborIter->second)->source->sID<<"->"<<
//					boost::tuples::get<1>(neighborIter->second)->sink->sID<<") "
//					<<" "<<boost::tuples::get<1>(neighborIter->second)->sourceTfactor
//					<<" "<<boost::tuples::get<1>(neighborIter->second)->sinkTfactor<<std::endl;
//					break;
//			}
//		}


//		for (constNeighborIteratorType neighborIter=this->outNeighborhood().begin();neighborIter!=this->outNeighborhood().end();++neighborIter){
//			boost::tuples::get<1>(neighborIter->second)->sourceTfactor=edgeConfiguration(kk);
//			kk++;
//		}	
//		for (constNeighborIteratorType neighborIter=this-> inNeighborhood().begin();neighborIter!=this-> inNeighborhood().end();++neighborIter){
//			boost::tuples::get<1>(neighborIter->second)->  sinkTfactor=edgeConfiguration(kk);
//			kk++;
//		}


//		}
//		else {
//			// set sourceTfactor/sinkTfactor to 0
//			for (constNeighborIteratorType neighborIter=this->outNeighborhood().begin();neighborIter!=this->outNeighborhood().end();++neighborIter){
//				boost::tuples::get<1>(neighborIter->second)->sourceTfactor=0.0;
//			}	
//			for (constNeighborIteratorType neighborIter=this-> inNeighborhood().begin();neighborIter!=this-> inNeighborhood().end();++neighborIter){
//				boost::tuples::get<1>(neighborIter->second)->  sinkTfactor=0.0;
//			}
//		}


//			if(this->neighborhood().size()==4){
//				std::cout<<"Node "<<this->sID<<": edgeConfiguration = "<<edgeConfiguration<<std::endl;
//			}


//				 eIter(ELev.begin());
//				for (std::map<double,Eigen::VectorXi>::const_iterator eIter=ELev.begin();eIter!=ELev.end()){
//				
//					if(eIter==ELev.begin()){
//						assert(std::fabs(eIter->first)<FLT_EPSILON && "FIRST ENERGY LEVEL MUST BE 0.");
//					}
//					else{
//						std::map<double,Eigen::VectorXi>::const_iterator nextEiter(eIter);
//						++nextEiter;
//						
//						if ()
//					
//					}
//				
//				}


// Select the energy level		

//				switch (ELev.size()) {
//					case 0:
//						assert(0 && "MAYBE THIS HAPPENS WHEN A NODE IS LEFT ISOLATED??");
//						break;
//					case 1:
//						assert(0 && "A SINGLE ENERGY LEVEL CANNOT BE FOUND FOR A BALANCED NODE");
//						break;
//					default:{
//						std::map<double,Eigen::VectorXi>::const_iterator firstEnergyIter(ELev.begin());
//						assert(std::fabs(firstEnergyIter->first)<FLT_EPSILON && "FIRST ENERGY LEVEL MUST BE 0.");
//						
//						std::map<double,Eigen::VectorXi>::const_iterator secondEnergyIter(firstEnergyIter);
//						secondEnergyIter++;
//						std::cout<<"Second Energy Level is"<<secondEnergyIter->first<<std::endl;
//						assert(std::fabs(secondEnergyIter->first)>FLT_EPSILON && "SECOND ENERGY LEVEL MUST BE > 0.");
//						
//						if (ELev.size()==2){
//							edgeConfiguration=secondEnergyIter->second;
////							energyLevelsAreValid=true;
//						}
//						else{
//							std::map<double,Eigen::VectorXi>::const_iterator thirdEnergyIter(secondEnergyIter);
//							thirdEnergyIter++;
//							if(std::fabs(thirdEnergyIter->first-secondEnergyIter->first)>FLT_EPSILON){ // second and third energy levels must be distinct
//								edgeConfiguration=secondEnergyIter->second;
////								energyLevelsAreValid=true;
//							}
//						}
//					}
//						break;
//				}



//	//////////////////////////////////////////////////
//	// make_CRN2N
//	Eigen::MatrixXd CRNs2Ns(Eigen::MatrixXd ENDRULE){
//		
//		Eigen::MatrixXd CRN2N;	// CatmulRom Neighbors to Neigbors matrix
//		CRN2N.setZero(ENDRULE.cols(),dim*CRneighbors.size());
//		
//		
//		Eigen::MatrixXd W2HofCentral;	
//		NeighborContainerType NeighborOFCentral;
//		
//		int curRow=0;
//		int curCol;
//		
//		for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
//			//std::cout<<"node "<<boost::tuples::get<0>(NeighborIter->second)->sID<<std::endl;
//			switch (boost::tuples::get<0>(NeighborIter->second)->is_central()) {
//				case 0:
//					curCol=dim*std::distance(CRneighbors.begin(), CRneighbors.find(boost::tuples::get<0>(NeighborIter->second)->sID) );
//					CRN2N.template block<dim,dim>(curRow,curCol).setIdentity();
//					curRow+=dim;
//					break;
//				case 1:
//					//! since the creation of source/sink nodes requires W2H from neighbor central nodes (which are not updated yet) update them
//					//boost::tuples::get<0>(NeighborIter->second)->update(); 
//					boost::tuples::get<0>(NeighborIter->second)->make_T(); 
//					W2HofCentral=boost::tuples::get<0>(NeighborIter->second)->W2H();
//					NeighborOFCentral=boost::tuples::get<0>(NeighborIter->second)->neighborhood();
//					int k=0;
//					for (NeighborIteratorType NeighborIter2=NeighborOFCentral.begin();NeighborIter2!=NeighborOFCentral.end();++NeighborIter2){
//						curCol=dim*std::distance(CRneighbors.begin(), CRneighbors.find(boost::tuples::get<0>(NeighborIter2->second)->sID) );
//						CRN2N.template block<2*dim,dim>(curRow,curCol)=W2HofCentral.template block<2*dim,dim>(0,dim*k);
//						++k;
//					}
//					curRow+=2*dim;
//					break;
//					
//			}
//			
//			
//		}
//		
//		return CRN2N;
//	}




//	//////////////////////////////////////////////////
//	// make_CR2H_source
//	void make_CR2H_source(){
//		
//		int N= this->Neighborhood.size()-1;
//		int C=0;
//		int S=0;
//		
//		int chi=1+!(N>1);
//		int SOURCERULEcols=0;
//		double gT=0.0;
//		double gCinv=0.0;
//		double gSinv=0.0;
//		
//		
//		
//		for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
//			
//			if (!boost::tuples::get<2>(NeighborIter->second)==0){
//				
//				gT+=boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//				
//				if (boost::tuples::get<0>(NeighborIter->second)->is_central()){
//					SOURCERULEcols+=2*dim;
//					gCinv+=1.0/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//					C+=1;
//				}
//				else if (boost::tuples::get<0>(NeighborIter->second)->is_sink()){
//					SOURCERULEcols+=dim;
//					gSinv+=1.0/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//					S+=1;
//				}
//				else{
//					assert(0);
//				}
//				
//			}
//			else{ //yourself
//				SOURCERULEcols+=dim;
//			}
//		}
//		
//		assert(N==C+S);
//		
//		Eigen::Matrix<double,2*dim,Eigen::Dynamic> SOURCERULE;
//		SOURCERULE.setZero(2*dim, SOURCERULEcols);
//		
//		
//		int currentCol=0;
//		for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
//			switch (boost::tuples::get<2>(NeighborIter->second)) {
//				case 0: // yourself
//					SOURCERULE.template block<dim,dim>(0,currentCol).setIdentity()*= (chi*N-1);
//					SOURCERULE.template block<dim,dim>(dim,currentCol)=this->prjM * ( (1.5*C+S)/gT - chi*(1.5*gCinv+gSinv) );
//					currentCol+=dim;
//					break;
//				default: // other
//					
//					if (boost::tuples::get<0>(NeighborIter->second)->is_central()){	// central
//						SOURCERULE.template block<dim,dim>(dim,currentCol)	  =this->prjM * 1.5*( chi/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength() - 1.0/gT );
//						SOURCERULE.template block<dim,dim>(dim,currentCol+dim)=this->prjM * 0.5*( boost::tuples::get<1>(NeighborIter->second)->chordParametricLength()/gT - chi );
//						currentCol+=2*dim;
//					}
//					else if (boost::tuples::get<0>(NeighborIter->second)->is_sink()){ // sink
//						SOURCERULE.template block<dim,dim>(dim,currentCol)=this->prjM * ( chi/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength() - 1.0/gT );
//						currentCol+=dim;
//					}
//					else{
//						assert(0);
//					}
//					break;
//			}
//		}
//		SOURCERULE/=(chi*N-1);
//		
//		
//		
//		
//		CR2H=SOURCERULE*CRNs2Ns(SOURCERULE);
//		
//		
//		
//		
//	}





//	//////////////////////////////////////////////////
//	// make_CR2H_sink
//	void make_CR2H_sink(){
//		
//		int N= this->Neighborhood.size()-1;
//		int C=0;
//		int S=0;
//		
//		int chi=1+!(N>1);
//		int SINKRULEcols=0;
//		double gT=0.0;
//		double gCinv=0.0;
//		double gSinv=0.0;
//		
//		
//		
//		for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
//			
//			if (!boost::tuples::get<2>(NeighborIter->second)==0){
//				
//				gT+=boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//				
//				if (boost::tuples::get<0>(NeighborIter->second)->is_central()){
//					SINKRULEcols+=2*dim;
//					gCinv+=1.0/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//					C+=1;
//				}
//				else if (boost::tuples::get<0>(NeighborIter->second)->is_source()){
//					SINKRULEcols+=dim;
//					gSinv+=1.0/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength();
//					S+=1;
//				}
//				else{
//					assert(0);
//				}
//				
//			}
//			else{ //yourself
//				SINKRULEcols+=dim;
//			}
//		}
//		
//		assert(N==C+S);
//		
//		Eigen::Matrix<double,2*dim,Eigen::Dynamic> SINKRULE;
//		SINKRULE.setZero(2*dim, SINKRULEcols);
//		
//		
//		int currentCol=0;
//		for (NeighborIteratorType NeighborIter=this->Neighborhood.begin();NeighborIter!=this->Neighborhood.end();++NeighborIter){
//			switch (boost::tuples::get<2>(NeighborIter->second)) {
//				case 0: // yourself
//					SINKRULE.template block<dim,dim>(0,currentCol).setIdentity()*= (chi*N-1);
//					SINKRULE.template block<dim,dim>(dim,currentCol)=this->prjM * ( chi*(1.5*gCinv+gSinv) - (1.5*C+S)/gT );
//					currentCol+=dim;
//					break;
//				default: // other
//					
//					if (boost::tuples::get<0>(NeighborIter->second)->is_central()){	// central
//						SINKRULE.template block<dim,dim>(dim,currentCol)	  =this->prjM * 1.5*(  1.0/gT - chi/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength() );
//						SINKRULE.template block<dim,dim>(dim,currentCol+dim)  =this->prjM * 0.5*( boost::tuples::get<1>(NeighborIter->second)->chordParametricLength()/gT - chi );
//						currentCol+=2*dim;
//					}
//					else if (boost::tuples::get<0>(NeighborIter->second)->is_source()){ // source
//						SINKRULE.template block<dim,dim>(dim,currentCol)=this->prjM * ( 1.0/gT - chi/boost::tuples::get<1>(NeighborIter->second)->chordParametricLength() );
//						currentCol+=dim;
//					}
//					else{
//						assert(0);
//					}
//					break;
//			}
//		}
//		SINKRULE/=(chi*N-1);
//		
//		
//		CR2H=SINKRULE*CRNs2Ns(SINKRULE);
//		
//	}





//	/* set ************************************/
//	void update(){		
//		
//		// MAKE THIS SET AND CHANGE TANGENTS ON NEIGHBORS
//		
//		make_T();
//	}

//	//////////////////////////////////////////////////
//	// Constructor with position
//	SplineNodeBase(const InitializationData<Derived>& idV) : NetworkNode<Derived>::NetworkNode(idV), P(idV.Q){
//		
//		set(get_nodeDof());
//		prjM.setIdentity();
//	}

//		//////////////////////////////////////////////////
//		// make_CR2H_source
//		void make_CR2H_free(){
//			CR2H.setZero(2*dim,dim);
//			CR2H.template block<dim,dim>(0,0).setIdentity();
//		}


//		//////////////////////////////////////////////////
//		// get_T
//		const VectorDim & get_T() const {
////			make_T();
//			return this->T;
//		}


//	std::cout<<"I'm here 2"<<std::endl;
//	std::cout<<SOURCERULE<<std::endl;
//	std::cout<<std::endl;
//	std::cout<<CRNs2Ns(SOURCERULE)<<std::endl;
//CR2H=SOURCERULE*CRN2N;
