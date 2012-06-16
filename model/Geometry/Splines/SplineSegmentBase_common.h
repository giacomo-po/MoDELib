/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */





/*!
 A SplineSegmentBase is a special kind of parametric curve where the position vector can be written
 as:
 \verbatim
 r(u)=UPOW(u)*H*q
 \endverbatim
 
 Each shape function is a pOrder-polynomial:
 \verbatim
 [N1 ... Nn]= [1 u ... n^porder] * SFC
 \endverbatim
 where SFC is a n x n matrix of shape function coefficients that in general dependes on the knot vector. It's posible
 to take advantage of this particular form of shape functions to calculate ruu and ru therefore these functions are
 redefined here (virtual in ParametricCurve) and the flow chart of their implementation is shown below:
 \verbatim
 make_UPOWuu()	  make_UPOWu()	make_UPOW()
 ^				  ^				^
 |				  |				| 
 make_SFuu()	  make_SFu()	make_SF()
 ^				  ^				^
 |				  |				|
 make_ruu()		  make_ru()		make_r()	(this level is CRTP in ParametricCurve)
 \endverbatim
 */





public:
EIGEN_MAKE_ALIGNED_OPERATOR_NEW


typedef typename NetworkLink<Derived>::NodeType NodeType;
typedef typename NetworkLink<Derived>::LinkType LinkType;
typedef typename NetworkLink<Derived>::FlowType FlowType;

//typedef model::SplineSegmentBase<Derived,dim,corder,alpha,qOrder,QuadratureRule> SplineSegmentBaseType;
typedef model::SplineSegmentBase<Derived,dim,corder,alpha> SplineSegmentBaseType;


typedef model::ParametricCurve<SplineSegmentBaseType,dim> ParametricCurveType;

//#include<model/Geometry/Splines/SplineEnums.h>



//protected:
//
//Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;

int sourceTfactor; // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of source node, which happens before constructor of this
int sinkTfactor;   // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of   sink node, which happens before constructor of this

VectorDim sourceT() const {
	return this->source->get_T()*sourceTfactor;
}

VectorDim sinkT() const {
	//AFTER INTRODUCING THE ENERGY CRITERION CHANGE THIS IN A -1
	return -this->sink->get_T()*sinkTfactor;
	//return this->sink->get_T()*sinkTfactor;
}

//////////////////////////////////////////////////////////////
protected:
//////////////////////////////////////////////////////////////		
//model::SplineIntersection<dim,corder> SI;

// Eigen::Matrix<double, Ndof, Eigen::Dynamic> G2H;	
















public:
//typedef typename TypeTraits<Derived>::FlowType FlowType;

/////////////////////////////////////////////////////////////////////////////////
// Constructor with Subnetwork* pair of Nodes*
//template <typename FlowType>
SplineSegmentBase(const std::pair<NodeType*,NodeType*> & nodePair_in, 
/*             */ const FlowType& Fin) : NetworkLink<Derived>::NetworkLink(nodePair_in, Fin){}
///*                                    */ sourceTfactor(1),
///*                                    */ sinkTfactor(1){}


SplineSegmentBase(const std::pair<NodeType*,NodeType*> & nodePair_in, 
/*             */ const ExpandingEdge<LinkType>& ee) : NetworkLink<Derived>::NetworkLink(nodePair_in, ee){}
///*             */ const LinkType* const pL) : NetworkLink<Derived>::NetworkLink(nodePair_in, pL){}
///*                                        */ sourceTfactor(1),
///*                                        */ sinkTfactor(1){}


//////////////////////////////////////////////////////////////
//! Returns the chord vector (source -> sink)
VectorDim chord() const {
	return this->sink->get_P()-this->source->get_P();
}

//////////////////////////////////////////////////////////////
//! Returns the length of the chord vector
double chordLength() const {
	return chord().norm();
}

//////////////////////////////////////////////////////////////
//! Returns the length of the chord vector to the power alpha
double chordParametricLength() const {
	return std::pow(chordLength(),alpha);;
}


/////////////////////////////////////////////////////////////////////////////////
// SF, SFu, SFuu
RowNcoeff get_SF(const double & uin) const {
	//make_UPOW(uin);
	//	std::cout<<"SplineSegmentBase::make_SF"<<std::endl;
	return get_UPOW(uin)*get_SFCH();
}

RowNcoeff get_SFu(const double & uin) const {
	
	//make_UPOWu(uin);
	
	//	std::cout<<"SplineSegmentBase::make_SF"<<std::endl;
	return get_UPOWu(uin)*get_SFCH().template block<Ncoeff-1,Ncoeff>(1,0);
	
}

RowNcoeff get_SFuu(const double & uin) const {
	//make_UPOWuu(uin);
	//	std::cout<<"SplineSegmentBase::make_SF"<<std::endl;
	return  get_UPOWuu(uin)*get_SFCH().template block<Ncoeff-2,Ncoeff>(2,0);
}

//////////////////////////////////////////////////////////////
//! SF*qH
VectorDim get_r(const double & u) const {
	/*! The position vector at parameter u
	 *  @param[in] u the parametrization variable in [0:1]
	 *	\f[
	 *		\mathbf{r} = \mathbf{q}\mathbf{H}\mathbf{u} \\
	 *		r_i = q_{ik}H_{km}u^{m} = q_{ik} N_k
	 *	\f]
	 * with i=0..dim-1, k,m = 0... Ncoeff
	 * ACTUALLY IN THE CODE WE HAVE THE TRANSPOSE OF THIS !!!!
	 */
	return get_SF(u)*get_qH();
}

//void make_r(const int & k){
//	make_r(this->abscissa(k));
//	SFgauss.row(k)=SF;
//}

//////////////////////////////////////////////////////////////
//! ru=SFu*qH
VectorDim get_ru(const double & uin) const {
	//make_SFu(uin);
	return get_SFu(uin)*get_qH();
}

//////////////////////////////////////////////////////////////
//! get_rmu
VectorDim get_rmu(const double & uin) const {
	return this->get_ru(uin)/chordParametricLength();
}

//////////////////////////////////////////////////////////////
//! ruu=SFuu*qH
VectorDim get_ruu(const double & uin) const {
	//make_SFuu(uin);
	return get_SFuu(uin)*get_qH();
}

//////////////////////////////////////////////////////////////
//! get_rmumu
VectorDim get_rmumu(const double & uin) const {
	return this->get_ruu(uin)/std::pow(chordParametricLength(),2);
}

//////////////////////////////////////////////////////////////
//! get_G2H
Eigen::Matrix<double, Ndof, Eigen::Dynamic>  get_G2H() const {
	//make_G2H();
	
	size_t gDof(this->pSN()->nodeOrder()*this->source->NdofXnode); // CHANGE HERE, NdofXnode should be available directly
	
	Eigen::Matrix<double, Ndof, Eigen::Dynamic> G2H(Eigen::Matrix<double, Ndof, Eigen::Dynamic>::Zero(Ndof,gDof));	
	
	//G2H.setZero(Ndof,gDof);
	
	Eigen::VectorXi dofid(this->source->dofID());
	Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> M(this->source->W2H());
	M.block(dim,0,dim,M.cols())*=sourceTfactor;
//	std::cout<<"M source=\n"<<M<<std::endl;
	
	for (int k=0;k<dofid.size();++k){
		G2H.template block<Ndof/2,1>(0,dofid(k))=M.col(k);
	}
	
	dofid=this->sink->dofID();
	M=this->sink->W2H();
	M.block(dim,0,dim,M.cols())*=(-sinkTfactor);
//	std::cout<<"M sink=\n"<<M<<std::endl;
	
	for (int k=0;k<dofid.size();++k){
		G2H.template block<Ndof/2,1>(Ndof/2,dofid(k))=M.col(k);
	}
	
	return G2H;
}

//////////////////////////////////////////////////////////////
//! intersects
//std::set<std::pair<double,double> > intersectWith( const Derived* const p_other  , const double& tol=FLT_EPSILON) const {
//	
//	
//	
//	std::set<std::pair<double,double> > intersectionParameters;
//	
//    VectorDim center1(0.5*(this->source->get_P()+this->sink->get_P()));
//    VectorDim center2(0.5*(p_other->source->get_P()+p_other->sink->get_P()));
//    double R1((this->source->get_P()-center1).norm());
//    double R2((p_other->source->get_P()-center2).norm());
//    double dist((center1-center2).norm());
//    
//    if(dist < 1.5*(R1+R2)){
//		VectorDim unUsed;		
//		assert((SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol)));
//		assert((SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff(),unUsed,tol)));
//		
//		//    std::cout<<"intersecting spheres"<<std::endl;
//  //       int degThis = SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol) + SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol);
//   //     int degOther= SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff(),unUsed,tol) + SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff(),tol);
//        
//		//		std::cout<<"This is "<< this->source->sID << "->" <<this->sink->sID<<std::endl;
//		//		std::cout<<"This is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol)<<std::endl;
//		//		std::cout<<"This is line? "<<SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol)<<std::endl;
//		////		
//		////		//			std::cout<<p_other->hermiteCoefficients()<<std::endl;
//		//		std::cout<<"Other is "<< p_other->source->sID << "->" <<p_other->sink->sID<<std::endl;
//		//		std::cout<<"Other is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff()   ,unUsed,tol)<<std::endl;
//		//		std::cout<<"Other is line? "<<SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff()   ,tol)<<std::endl;
//		//		
//		//		
//		//		std::cout<<"degThis= "<<degThis<<std::endl;
//		//		std::cout<<"degOther= "<<degOther<<std::endl;        
//		//        if(degThis==0 && degOther==0){
//		//			
//		//
//		//        }
//	//	if(degThis==1 && degOther==1){
//            SplineIntersection<dim,pOrder,1,1> SI(this->polynomialCoeff(),p_other->polynomialCoeff(),tol);
//            intersectionParameters=SI.intersectionParameters;
//			//			std::cout<<"Case 1 1"<<std::endl;
//      //  }
////        else if(degThis==1 && degOther==2){
////            SplineIntersection<dim,pOrder,1,2> SI(this->polynomialCoeff(),p_other->source->get_P(),p_other->sink->get_P(),tol);
////            intersectionParameters=SI.intersectionParameters;
////			//						std::cout<<"Case 1 2"<<std::endl;
////        }
////        else if(degThis==2 && degOther==1){
////            SplineIntersection<dim,pOrder,1,2> SI(p_other->polynomialCoeff(),this->source->get_P(),this->sink->get_P(),tol);
////            std::set<std::pair<double,double> > temp=SI.intersectionParameters;
////            for(std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter){ //switch first-second
////                intersectionParameters.insert(std::make_pair(paramIter->second,paramIter->first));
////            }
////			//						std::cout<<"Case 2 1"<<std::endl;
////        }
////        else if(degThis==2 && degOther==2){
////            SplineIntersection<dim,pOrder,2,2> SI(this->source->get_P(),this->sink->get_P(),p_other->source->get_P(),p_other->sink->get_P(),tol);
////            intersectionParameters=SI.intersectionParameters;
////			//						std::cout<<"Case 2 2"<<std::endl;
////        }
////        else{
////			
////			//			std::cout<<this->hermiteCoefficients()<<std::endl;
////			//			std::cout<<"This is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////			//			std::cout<<"This is line? "<<SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol)<<std::endl;
////			//			
////			//			//			std::cout<<p_other->hermiteCoefficients()<<std::endl;
////			//			std::cout<<"Other is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////			//			std::cout<<"Other is line? "<<SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff()   ,tol)<<std::endl;
////			
////			
////			
////            assert(0 && "NOT READY YET");
////			
////            
////        }
//        
//        
//		
//        
//	}
//	return intersectionParameters;
//	//return SI.intersect(this->polynomialCoeff(),p_other->polynomialCoeff(),Alpha);
//}




//std::set<std::pair<double,double> > intersectWith( const Derived* const p_other  , const double& tol=FLT_EPSILON) const {
//	
//	
//	
//	std::set<std::pair<double,double> > intersectionParameters;
//	
//    VectorDim center1(0.5*(this->source->get_P()+this->sink->get_P()));
//    VectorDim center2(0.5*(p_other->source->get_P()+p_other->sink->get_P()));
//    double R1((this->source->get_P()-center1).norm());
//    double R2((p_other->source->get_P()-center2).norm());
//    double dist((center1-center2).norm());
//    
//    if(dist < 1.5*(R1+R2)){
//    //    std::cout<<"intersecting spheres"<<std::endl;
//        VectorDim unUsed;
//        int degThis = SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol) + SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol);
//        int degOther= SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff(),unUsed,tol) + SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff(),tol);
//        
////		std::cout<<"This is "<< this->source->sID << "->" <<this->sink->sID<<std::endl;
////		std::cout<<"This is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////		std::cout<<"This is line? "<<SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol)<<std::endl;
//////		
//////		//			std::cout<<p_other->hermiteCoefficients()<<std::endl;
////		std::cout<<"Other is "<< p_other->source->sID << "->" <<p_other->sink->sID<<std::endl;
////		std::cout<<"Other is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////		std::cout<<"Other is line? "<<SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff()   ,tol)<<std::endl;
////		
////		
////		std::cout<<"degThis= "<<degThis<<std::endl;
////		std::cout<<"degOther= "<<degOther<<std::endl;        
////        if(degThis==0 && degOther==0){
////			
////
////        }
//		if(degThis==1 && degOther==1){
//            SplineIntersection<dim,pOrder,1,1> SI(this->polynomialCoeff(),p_other->polynomialCoeff(),tol);
//            intersectionParameters=SI.intersectionParameters;
////			std::cout<<"Case 1 1"<<std::endl;
//        }
//        else if(degThis==1 && degOther==2){
//            SplineIntersection<dim,pOrder,1,2> SI(this->polynomialCoeff(),p_other->source->get_P(),p_other->sink->get_P(),tol);
//            intersectionParameters=SI.intersectionParameters;
////						std::cout<<"Case 1 2"<<std::endl;
//        }
//        else if(degThis==2 && degOther==1){
//            SplineIntersection<dim,pOrder,1,2> SI(p_other->polynomialCoeff(),this->source->get_P(),this->sink->get_P(),tol);
//            std::set<std::pair<double,double> > temp=SI.intersectionParameters;
//            for(std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paramIter!=temp.end();++paramIter){ //switch first-second
//                intersectionParameters.insert(std::make_pair(paramIter->second,paramIter->first));
//            }
////						std::cout<<"Case 2 1"<<std::endl;
//        }
//        else if(degThis==2 && degOther==2){
//            SplineIntersection<dim,pOrder,2,2> SI(this->source->get_P(),this->sink->get_P(),p_other->source->get_P(),p_other->sink->get_P(),tol);
//            intersectionParameters=SI.intersectionParameters;
////						std::cout<<"Case 2 2"<<std::endl;
//        }
//        else{
//			
//			//			std::cout<<this->hermiteCoefficients()<<std::endl;
////			std::cout<<"This is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(this->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////			std::cout<<"This is line? "<<SplineDegeneracy<dim,pOrder>::isLine(this->polynomialCoeff()   ,tol)<<std::endl;
////			
////			//			std::cout<<p_other->hermiteCoefficients()<<std::endl;
////			std::cout<<"Other is planar? "<<SplineDegeneracy<dim,pOrder>::isPlanar(p_other->polynomialCoeff()   ,unUsed,tol)<<std::endl;
////			std::cout<<"Other is line? "<<SplineDegeneracy<dim,pOrder>::isLine(p_other->polynomialCoeff()   ,tol)<<std::endl;
//			
//			
//			
//            assert(0 && "NOT READY YET");
//			
//            
//        }
//        
//        
//
//        
//	}
//	return intersectionParameters;
//	//return SI.intersect(this->polynomialCoeff(),p_other->polynomialCoeff(),Alpha);
//}





/////////////////////////////////////////////////////////////////////////////////
//! make_G2H
//void make_G2H(){
//	
//	
//	size_t gDof=this->pSN()->nodeOrder()*this->source->NdofXnode;
//	
//	Eigen::Matrix<double, Ndof, Eigen::Dynamic> G2H;	
//
//	G2H.setZero(Ndof,gDof);
//	
//	Eigen::VectorXi dofid=this->source->dofID();
//	Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> M=this->source->W2H();
//	
//	
//	
//	for (int k=0;k<dofid.size();++k){
//		G2H.template block<Ndof/2,1>(0,dofid(k))=M.col(k);
//	}
//	
//	dofid=this->sink->dofID();
//	M=this->sink->W2H();
//	
//	for (int k=0;k<dofid.size();++k){
//		G2H.template block<Ndof/2,1>(Ndof/2,dofid(k))=M.col(k);
//	}
//	
//	
//}





//	switch (std::make_pair(degThis,degOther)) {
//		case std::make_pair(2,2):{
//			SplneIntersection<dim,pOrder,2,2> SI(this->source->get_P(),this->sink->get_P(),p_other->source->get_P(),p_other->sink->get_P(),tol);
//			intersectionParameters=SI.intersectionParameters;
//			break;}
//			
//		case std::make_pair(1,1):{
////			SplneIntersection<dim,pOrder,1,1> SI(this->polynomialCoeff(),p_other->polynomialCoeff(),tol);
////			intersectionParameters=SI.intersectionParameters;
//			break;}
//			
//		case std::make_pair(1,2):{
////			SplneIntersection<dim,pOrder,1,2> SI(this->polynomialCoeff(),p_other->source->get_P(),p_other->sink->get_P(),tol);
////			intersectionParameters=SI.intersectionParameters;
//			break;}
//			
//			
//		case std::make_pair(2,1):{
//			SplneIntersection<dim,pOrder,1,2> SI(p_other->polynomialCoeff(),this->source->get_P(),this->sink->get_P(),tol);
//			std::set<std::pair<double,double> > temp=SI.intersectionParameters;
//			for(std::set<std::pair<double,double> >::const_iterator paramIter=temp.begin();paraItere!=temp.end();++paramIter){ //switch first-second
//				intersectionParameters.insert(std::make_pair(paramIter.second,paramIter.first));
//			}
//			break;}
//
//		default:
//			
//			break;
//	}


//RowNcoeff	UPOW, SF, SFu, SFuu;
//RowNcoeffu	UPOWu;
//RowNcoeffuu UPOWuu;

// Vector of Hermite DOF for the segment

//MatrixNcoeffDim qH;



//! Matrix of Hermite shape functions coefficients 
//MatrixNcoeff SFCH;	//! THIS SHOULD NOT BE STORED but calculated in get_SFCH(), maybe make template<alpha> ? 




/////////////////////////////////////////////////////////////////////////////////
// make_SFCH_qH
//	void make_SFCH_qH();

/////////////////////////////////////////////////////////////////////////////////
// UPOW, WPOWu, UPOWuu
//	void make_UPOW(const double & uin);
//	void make_UPOWu(const double & uin);
//	void make_UPOWuu(const double & uin);


//! 3- transmit the SubNetwork 
//	typedef void (NodeType::*node_member_function_pointer_type)(void); 
//	node_member_function_pointer_type Nmfp;
//	Nmfp=&NodeType::update;
//	
//	typedef void (Derived::*link_member_function_pointer_type)(void); 
//	link_member_function_pointer_type Lmfp;
//	Lmfp=&Derived::update;

//transmit(Nmfp,Lmfp,void);

// the problem here is that if you call source->transmit(Nmfp,Lmfp,void) then the source will update and then this link is
// updated and then the sink is updated. Therefore the link is updated before the sink. Now if the tangent of the sink is calculated 
// by update this means that the link is updated with the wrong sink tangent
// if get_T instead calculates the tabgent this is ok but there are a lot of redundant calculations


//! Trigger the transmit of update()
//this->source->set(this->source->get_nodeDof());
//this->sink  ->set(this->sink  ->get_nodeDof());

//update();



/////////////////////////////////////////////////////////////////////////////////
// Destructor
//~SplineSegmentBase(){
//	
//	//! Trigger the transmit of update()
//	
//	// THIS IS A PROBLEM HERE!!!!! because the connectivity is changed by the NetworkLink destructor (executed after this)
//	//	this->source->set(this->source->get_nodeDof());
//	//	this->sink  ->set(this->sink  ->get_nodeDof());
//	
//}


/////////////////////////////////////////////////////////////////////////////////
//! update		
//void update(){
//
////	std::cout<<"SplineSegment ("<<this->source->sID<<"->"<<this->sink->sID<< ") updating"<<std::endl;
//	
////	make_SFCH_qH();
//	
////	std::cout<<"Out of SplineSegmentBase make_SFCH_qH"<<std::endl;
//	
//	// resize SFgauss
//	//SFgauss.resize(this->get_QuadratureOrder(),Ncoeff);
//	
//	//! update the parametric curve
//	//! SplineSegmentBase_common !!!!! THIS CAN BE REMOVED !!!!!
//	//ParametricCurveType::update();	
//	
//	
//}


