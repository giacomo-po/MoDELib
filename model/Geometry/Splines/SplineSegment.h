/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SplineSegment_H_
#define model_SplineSegment_H_

#include <math.h>
#include <Eigen/Dense>




//#include <NetworkLink.h>
#include <NetworkLink.h>
#include <ParametricCurve.h>
//#include <SplineIntersection.h>
#include <Coeff2Hermite.h>
#include <SplineSegmentBase.h>
//#include <EdgeFinder.h>
//#include <Eigen/Sparse>
#include <MatrixCompanion.h>

namespace model
{
    
    /*!
     A SplineSegment is a special kind of parametric curve where the position vector can be written
     as:
     \verbatim
     r(u)=N(u)*q
     \endverbatim
     where N(u) is a (dim x NdofXsegment) matrix of shape functions and q is a (NdofXsegment x 1) column
     vector of degrees of freedom, where NdofXsegment= dim x n and n=porder+1 is the number of shape functions:
     \verbatim
     | N1  0  0		Nn  0  0 |
     N= |  0 N1  0	...  0 Nn  0 |
     |  0  0 N1		 0  0 Nn |
     \endverbatim
     
     Each shape function is a porder:
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
     make_Nuu()		  make_Nu()		make_N()
     ^				  ^				^
     |				  |				|
     make_ruu()		  make_ru()		make_r()	(this level is virtual in ParametricCurve)
     \endverbatim
     */
    
    
    
    
    /************************************************************************/
    /* SplineSegment, general case **************************************/
    /************************************************************************/
    template <short unsigned int _dim,short unsigned int corder>
    class SplineSegment : public ParametricCurve<SplineSegment<_dim,corder>,_dim>
    {
        
        static_assert(_dim>=1 && _dim <=3,"DIMENSION MUST BE 1, 2 or 3."); // requires c++11
        static_assert(corder>=0 && corder <=2,"CONTINUITY ORDER MUST BE 0, 1 or 2."); // requires c++11
        
        
    public:
        //enum  {dim=_dim};
        static constexpr int dim= _dim;
        static constexpr int Ncoeff= 2*(corder+1);      // number of Hermite coefficients
        static constexpr int pOrder= 2*corder+1;
        static constexpr int Ndof  = dim*Ncoeff;        // number of Hermite dofs
        static constexpr int eigenSize=pOrder*pOrder;
        
        typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
        typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
        typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
        typedef Eigen::Matrix<double, Ndof,1> VectorNdof;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double, dim, Ndof> MatrixDimNdof;
        typedef Eigen::Matrix<double,Ncoeff,1>     VectorNcoeff;
        typedef Eigen::Matrix<double,Ndof,Ndof>	MatrixNdof;
        typedef SplineSegmentBase<dim,corder> SplineSegmentBaseType;
        typedef ParametricCurve<SplineSegment<_dim,corder>,dim> ParametricCurveType;
        
        
        
        const VectorDim& sourceP;
        const VectorDim& sinkP;
        
    private:

        VectorDim _chord;
        double _chordLength2;
        double _chordLength;
        VectorDim _unitDirection;
        
    public:
        

        
        static double alpha;
        
        
        /**********************************************************************/
        SplineSegment(const VectorDim& sourceP_in,
                      const VectorDim& sinkP_in) :
        /* init */ sourceP(sourceP_in)
        /* init */,sinkP(sinkP_in)
        /* init */,_chord(sinkP-sourceP)
        /* init */,_chordLength2(_chord.squaredNorm())
        /* init */,_chordLength(sqrt(_chordLength2))
        /* init */,_unitDirection(_chordLength>FLT_EPSILON? (_chord/_chordLength).eval() : VectorDim::Zero())
        {/*! Constructor with Nodes and flow
          */
            
        }
        
        /**********************************************************************/
        void updateGeometry()
        {
            _chord=sinkP-sourceP;
            _chordLength2=_chord.squaredNorm();
            _chordLength=sqrt(_chordLength2);
            _unitDirection=_chordLength>FLT_EPSILON? (_chord/_chordLength).eval() : VectorDim::Zero();
        }
        
        /**********************************************************************/
        MatrixNcoeffDim hermiteDofs() const
        {
            return SplineSegmentBaseType::hermiteDofs(*this);
        }
        
        /**********************************************************************/
        std::pair<Eigen::Matrix<double,Ncoeff,Eigen::Dynamic>,Eigen::Matrix<double,Eigen::Dynamic,dim>> hermite2posMatrix() const
        {
            return SplineSegmentBaseType::hermite2posMatrix(*this);
        }
        
        /**********************************************************************/
        typename SplineSegmentBaseType::H2PmapType hermite2posMap() const
        {
            return SplineSegmentBaseType::hermite2posMap(*this);
        }

        /**********************************************************************/
        Eigen::Matrix<double,dim,Ncoeff> hermiteCoefficients() const
        {
            return SplineSegmentBaseType::hermiteCoefficients(*this);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,Ncoeff> polynomialCoeff() const
        {/*! The matrix of coefficients of the polynomial associated to this
          *  SplineSegmentBase. If C=polynomialCoeff() then the polynomial is:
          *  P(u)=C.col(0)+u*C.col(1)+u^2*C.col(2)+...
          */
            return Coeff2Hermite<pOrder>::template h2c<dim>(hermiteCoefficients());
        }
        
        /**********************************************************************/
        MatrixNcoeff sfCoeffs() const
        {
            return SplineSegmentBaseType::sfCoeffs(parametricChordLength());
        }
        
//        /**********************************************************************/
//        VectorDim chord() const
//        {/*!\returns the chord vector (source -> sink)
//          */
//            return this->sink->get_P()-this->source->get_P();
//        }
//        
//        /**********************************************************************/
//        double chordLength() const
//        {/*!\returns the length of the chord vector
//          */
//            return chord().norm();
//        }

        /**********************************************************************/
        const VectorDim& chord() const
        {/*!\returns the chord vector (source -> sink)
          */
            return _chord;
        }
        
        /**********************************************************************/
        const double& chordLength() const
        {/*!\returns the length of the chord vector
          */
            return _chordLength;
        }
        
        /**********************************************************************/
        const double& chordLengthSquared() const
        {/*!\returns the length of the chord vector
          */
            return _chordLength2;
        }
        
        /**********************************************************************/
        const VectorDim& unitDirection() const
        {/*!\returns the chord vector (source -> sink)
          */
            return _unitDirection;
        }
        
        /**********************************************************************/
        double parametricChordLength() const
        {//!\returns  the length of the chord vector to the power alpha
            return std::pow(chordLength(),alpha);
        }
        
        /**********************************************************************/
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
            return sf(u)*hermiteDofs();
        }
        
        /**********************************************************************/
        RowNcoeff sf(const double& u) const
        {/*\param[in] uin parameter value in [0,1]
          *\returns The row vector of shape function at uin
          */
            return SplineSegmentBaseType::sf(u,parametricChordLength());
        }
        
        /**********************************************************************/
        VectorDim get_ru(const double & uin) const
        {
            return SplineSegmentBaseType::sfDiff1(uin,parametricChordLength())*hermiteDofs();
        }
        
        /**********************************************************************/
        VectorDim get_rmu(const double & uin) const
        {
            return this->get_ru(uin)/parametricChordLength();
        }
        
        /**********************************************************************/
        VectorDim get_ruu(const double & uin) const
        {
            return SplineSegmentBaseType::sfDiff2(uin,parametricChordLength())*hermiteDofs();
        }
        
        /**********************************************************************/
        VectorDim get_rmumu(const double & uin) const
        {
            return this->get_ruu(uin)/std::pow(parametricChordLength(),2);
        }
        
        /**********************************************************************/
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
            std::map<double,std::pair<double,VectorDim>> rootMap;
            
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
            rootMap.insert(std::make_pair((sourceP-P0).norm(), std::make_pair(0.0,sourceP)));
            rootMap.insert(std::make_pair((sinkP-P0).norm(), std::make_pair(1.0,sinkP)));
            
            return *rootMap.begin();
            
        }
        
    };
    
    
    //static data
    template <short unsigned int dim,short unsigned int corder>
    double SplineSegment<dim,corder>::alpha=0.5;
    
}
#endif

