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

//typedef SplineSegmentBase<Derived,dim,corder,alpha,qOrder,QuadratureRule> SplineSegmentBaseType;
//typedef SplineSegmentBase<Derived,dim,corder,alpha> SplineSegmentBaseType;
typedef SplineSegmentBase<Derived,dim,corder> SplineSegmentBaseType;


typedef ParametricCurve<SplineSegmentBaseType,dim> ParametricCurveType;

//#include<model/Geometry/Splines/SplineEnums.h>



//protected:
//
//Eigen::Matrix<double,qOrder,Ncoeff> SFgauss;

int sourceTfactor; // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of source node, which happens before constructor of this
int sinkTfactor;   // LEAVE THIS UNINITIALIZED: this is calculated in TopologyChangeActions of   sink node, which happens before constructor of this

VectorDim sourceT() const
{
    return this->source->get_T()*sourceTfactor;
}

VectorDim sinkT() const
{
    //AFTER INTRODUCING THE ENERGY CRITERION CHANGE THIS IN A -1
    return -this->sink->get_T()*sinkTfactor;
    //return this->sink->get_T()*sinkTfactor;
}

public:

/******************************************************************************/
SplineSegmentBase(const std::pair<NodeType*,NodeType*> & nodePair_in,
/*             */ const FlowType& Fin) : NetworkLink<Derived>::NetworkLink(nodePair_in, Fin)
{/*! Constructor with Nodes and flow
  */
    
}

/******************************************************************************/
SplineSegmentBase(const std::pair<NodeType*,NodeType*> & nodePair_in,
/*             */ const EdgeRef<LinkType>& ee) :
/* init list */ NetworkLink<Derived>::NetworkLink(nodePair_in, ee)
{/*! Constructor from EdgeRef
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

/******************************************************************************/
VectorDim chord() const
{/*!\returns the chord vector (source -> sink)
  */
    return this->sink->get_P()-this->source->get_P();
}

/******************************************************************************/
double chordLength() const
{/*!\returns the length of the chord vector
  */
    return chord().norm();
}

/******************************************************************************/
RowNcoeff get_SF(const double & uin) const
{
    return get_UPOW(uin)*get_SFCH();
}

/******************************************************************************/
RowNcoeff get_SFu(const double & uin) const
{
    return get_UPOWu(uin)*get_SFCH().template block<Ncoeff-1,Ncoeff>(1,0);
}

/******************************************************************************/
RowNcoeff get_SFuu(const double & uin) const
{
    return  get_UPOWuu(uin)*get_SFCH().template block<Ncoeff-2,Ncoeff>(2,0);
}

/******************************************************************************/
VectorDim get_r(const double & u) const
{/*!\returns The position vector at parameter u
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

/******************************************************************************/
VectorDim get_ru(const double & uin) const
{
    return get_SFu(uin)*get_qH();
}

/******************************************************************************/
VectorDim get_rmu(const double & uin) const
{
    return this->get_ru(uin)/chordParametricLength();
}

/******************************************************************************/
VectorDim get_ruu(const double & uin) const
{
    return get_SFuu(uin)*get_qH();
}

/******************************************************************************/
VectorDim get_rmumu(const double & uin) const
{
    return this->get_ruu(uin)/std::pow(chordParametricLength(),2);
}

/******************************************************************************/
Eigen::Matrix<double, Ndof, Eigen::Dynamic>  get_G2H() const
{
    //make_G2H();
    
    size_t gDof(this->pSN()->nodeOrder()*this->source->NdofXnode); // CHANGE HERE, NdofXnode should be available directly
    
    Eigen::Matrix<double, Ndof, Eigen::Dynamic> G2H(Eigen::Matrix<double, Ndof, Eigen::Dynamic>::Zero(Ndof,gDof));
    
    //G2H.setZero(Ndof,gDof);
    
    Eigen::VectorXi dofid(this->source->dofID());
    Eigen::Matrix<double, Ndof/2, Eigen::Dynamic> M(this->source->W2H());
    M.block(dim,0,dim,M.cols())*=sourceTfactor;
    //	std::cout<<"M source=\n"<<M<<std::endl;
    
    for (int k=0;k<dofid.size();++k)
    {
        G2H.template block<Ndof/2,1>(0,dofid(k))=M.col(k);
    }
    
    dofid=this->sink->dofID();
    M=this->sink->W2H();
    M.block(dim,0,dim,M.cols())*=(-sinkTfactor);
    //	std::cout<<"M sink=\n"<<M<<std::endl;
    
    for (int k=0;k<dofid.size();++k)
    {
        G2H.template block<Ndof/2,1>(Ndof/2,dofid(k))=M.col(k);
    }
    
    return G2H;
}

/******************************************************************************/
std::pair<double,std::pair<double,VectorDim> > closestPoint(const VectorDim& P0) const
{/*!@param[in] P0 reference point
  * \returns The closesest point to P0 along this segment. The return value is a
  * pair, where pair.first is the parameter value u, and pair.second is the position P
  * of the closest point, so that get_r(u)=P.
  *
  * 
  */
    
    // solve (P-P0)*dP/du=0
    
    // The polynomial coefficients of this spine segment
    Eigen::Matrix<double,dim,Ncoeff> coeffs(this->polynomialCoeff());
    coeffs.col(0)-=P0;
    
    
    // The derivative of the polynomial coefficients
    Eigen::Matrix<double,dim,Ncoeff-1> dcoeffs(Eigen::Matrix<double,dim,Ncoeff-1>::Zero());
    for (int i=0;i<Ncoeff-1;++i)
    {
        dcoeffs.col(i)=(i+1)*coeffs.col(i+1);
    }
    
    Eigen::Matrix<double,1,2*Ncoeff-2> pcoeffs(Eigen::Matrix<double,1,2*Ncoeff-2>::Zero()); // degree of product = pOrder+(pOrder-1)=2*pOrder-1. nCoeffs of product = 2*pOrder-1+1= 2*pOrder = 2*Ncoeff-2
    
    // The polynomial coefficients of (P-P0)*dP/du, in reverse order
    for (int i=0;i<Ncoeff;++i)
    {
        for (int j=0;j<Ncoeff-1;++j)
        {
            pcoeffs(2*Ncoeff-3-i-j) += coeffs.col(i).dot(dcoeffs.col(j));
        }
    }
    
    // Compute roots using the eigenvalue method
    MatrixCompanion<2*Ncoeff-3> mc(pcoeffs);
    
    // sort roots according to distance to P0
    std::map<double,std::pair<double,VectorDim> > rootMap;
    
    //    for (int k=0;k<2*Ncoeff-3;++k)
    for (size_t k=0;k<mc.rootSize;++k)
    {
        if (std::fabs(mc.root(k).imag())<FLT_EPSILON && mc.root(k).real()>0.0 && mc.root(k).real()<1.0 )
        {
            
            VectorDim P(this->get_r(mc.root(k).real()));
            rootMap.insert(std::make_pair((P-P0).norm(), std::make_pair(mc.root(k).real(),P) ));
            
        }
        
    }
    
    // check distance to limits of interval
    rootMap.insert(std::make_pair((this->source->get_P()-P0).norm(), std::make_pair(0.0,this->source->get_P()) ));
    rootMap.insert(std::make_pair((this->  sink->get_P()-P0).norm(), std::make_pair(1.0,this->  sink->get_P()) ));
    
    return *rootMap.begin();
    
}
