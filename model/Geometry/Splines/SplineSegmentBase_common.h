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
/*             */ const FlowType& Fin) : NetworkLink<Derived>::NetworkLink(nodePair_in, Fin)
{/*! Constructor with Nodes and flow
  */

}

///*                                    */ sourceTfactor(1),
///*                                    */ sinkTfactor(1){}


SplineSegmentBase(const std::pair<NodeType*,NodeType*> & nodePair_in, 
/*             */ const ExpandingEdge<LinkType>& ee) :
/* init list */ NetworkLink<Derived>::NetworkLink(nodePair_in, ee)
{/*! Constructor from ExpandingEdge
  */
    
    
//    assert(this->source!=ee.E.sink);
//    assert(this->sink!=ee.E.source);
//    
//    if (this->source==ee.E.source) // first new link in expansion
//    {
//        sourceTfactor=ee.E.sourceTfactor;
//        sinkTfactor= -ee.E.sinkTfactor;
//    }
//    else // second new link in expansion
//    {
//        sourceTfactor=-ee.E.sourceTfactor;
//        sinkTfactor=ee.E.sinkTfactor;
//    }
    
    
//    if (this->source->sID==88 && this->sink->sID==1)
//    {
//        std::cout<<"sourceTfactor="<<sourceTfactor<<std::endl;
//        std::cout<<"  sinkTfactor="<<sinkTfactor<<std::endl;
//    }
    
}



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


/* isCommonNeighborAt *************************************************/
std::pair<bool,size_t> isCommonNeighborAt(const VectorDim& P0) const {
    std::pair<bool,size_t> temp(false,0);
    for (typename Derived::NeighborContainerType::const_iterator nIiter=this->source->neighborhood().begin();nIiter!=this->source->neighborhood().end();++nIiter){ // loop over neighborhood of source
        if (std::get<0>(nIiter->second)->sID!=this->source->sID && std::get<0>(nIiter->second)->sID!=this->sink->sID){ // neighbor is neither source nor sink
            if((std::get<0>(nIiter->second)->get_P()-P0).norm()<FLT_EPSILON){ // a neighbor of I exists at P0
                const size_t p0ID(std::get<0>(nIiter->second)->sID); // the sID of the neighbor at P0
                for (typename Derived::NeighborContainerType::const_iterator nJiter=this->sink->neighborhood().begin();nJiter!=this->sink->neighborhood().end();++nJiter){ // loop over neighborhood of sink
                    if(std::get<0>(nJiter->second)->sID==p0ID){
                        temp=std::pair<bool,size_t>(true,p0ID);
                    }
                }
            }
        }
    }
    return temp;
}