/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SPLINESEGMENTBASE_H_
#define model_SPLINESEGMENTBASE_H_

#include <math.h>
#include <Eigen/Dense>



#include <model/Geometry/Splines/SplineConsts.h>

#include <model/Network/NetworkLink.h>
#include <model/Geometry/ParametricCurve.h>
//#include <model/Geometry/Splines/Intersection/SplineIntersection.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Network/Operations/EdgeFinder.h>
#include <Eigen/Sparse>
#include <model/Math/MatrixCompanion.h>

namespace model {
	
	/*!
	 A SplineSegmentBase is a special kind of parametric curve where the position vector can be written
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
	/* SplineSegmentBase, general case **************************************/
	/************************************************************************/
	template <typename Derived, short unsigned int dim, short unsigned int corder>
	class SplineSegmentBase {};
	
	
	/************************************************************************/
	/* SplineSegmentBase, template specialization corder=0 ******************/
	/************************************************************************/
	template <typename Derived, short unsigned int dim> 
	class SplineSegmentBase<Derived,dim,0> : public NetworkLink<Derived>,
	/*	                                        */ public ParametricCurve<Derived, dim>
    {
		
		
		enum {corder=0};	
#include<model/Geometry/Splines/SplineEnums.h>
		
		
        ////////////////////////////////////////////////////////////////
        //! Returns the length of the chord vector to the power alpha
        double chordParametricLength() const
        {
        	return chordLength();
        }
        
	public:
		RowNcoeff get_UPOW(const double & uin) const {
			return (RowNcoeff()<<1.0, uin).finished();
		}
		
		//////////////////////////////////////////////////////////////
		
		RowNcoeffu get_UPOWu(const double & uin) const {
			return (RowNcoeffu()<<1.0).finished();
		}
		

		
		RowNcoeffuu get_UPOWuu(const double & uin) const {
			// ????????????????????????
		}
		
		static MatrixNcoeff get_SFCH()
        {
			/*! The matrix of shape function coefficients in Hermite form of this 
			 *  spline segment.
			 */
			/*                         P0   P1 */
			return (MatrixNcoeff()<<  1.0, 0.0,             // u^0
					/*            */ -1.0, 1.0).finished(); // u^1;
		}
		
		MatrixNcoeffDim get_qH() const {
			return (MatrixNcoeffDim()<< this->source->get_P().transpose(), 
					/*               */ this->  sink->get_P().transpose()).finished();
		}
		

		
		//////////////////////////////////////////////////////////////
		//BezierCoefficients: Bezier coefficients (uniform parametrization)
		Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
			Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
			double w =1.0;
			BzCf.col(0)<< this->source->get_P(), w;
			BzCf.col(2)<< this->sink  ->get_P(), w;
			return BzCf;
		}
		
#include "SplineSegmentBase_common.h"
		
	};
	
	
	
	
	
	
	
	/************************************************************************/
	/* SplineSegmentBase, template specialization corder=1 ******************/
	/************************************************************************/
	template <typename Derived, short unsigned int dim>
	class SplineSegmentBase<Derived,dim,1> : public NetworkLink<Derived>,
	/*	                                        */ public ParametricCurve<Derived,dim> {
		
				
		enum {corder=1};
#include<model/Geometry/Splines/SplineEnums.h>
		

		
	public:
        
        
        static double alpha;

		
        ////////////////////////////////////////////////////////////////
        //! Returns the length of the chord vector to the power alpha
        double chordParametricLength() const
        {
        	return std::pow(chordLength(),alpha);;
        }
		
		
		RowNcoeff get_UPOW(const double & uin) const {
			return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3)).finished();
		}
		
		
		
		RowNcoeffu get_UPOWu(const double & uin) const {
			return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2)).finished();
		}
		
		
		RowNcoeffuu get_UPOWuu(const double & uin) const {
			return (RowNcoeffuu()<<2.0, 6.0*uin).finished();
		}
		
		
		//////////////////////////////////////////////////////////////
		//get_SFCH
		 MatrixNcoeff get_SFCH() const
        {
			/*! The matrix of shape function coefficients in Hermite form of this 
			 *  spline segment.
			 */
			const double g(chordParametricLength());
			/*                         P0      T0    P1   T1  */
			return (MatrixNcoeff()<<  1.0,    0.0,  0.0, 0.0,             // u^0
					/*            */  0.0,      g,  0.0, 0.0,             // u^1
					/*            */ -3.0, -2.0*g,  3.0,  -g,             // u^2
					/*            */  2.0,      g, -2.0,   g).finished(); // u^3;
		}
		
		//////////////////////////////////////////////////////////////
		//get_qH
		MatrixNcoeffDim get_qH() const {
			/*! The matrix of Hermite dof of this spline segment.
			 *  [P0x P0y P0z;T0x T0y T0z;P1x P1y P1z;T1x T1y T1z]
			 */
			return (MatrixNcoeffDim()<< this->source->get_P().transpose(),
					/*            */	sourceT().transpose(),
					/*            */	this->  sink->get_P().transpose(),
					/*            */	sinkT().transpose()).finished();
		}
		
		//////////////////////////////////////////////////////////////
		//BezierCoefficients: Bezier coefficients (uniform parametrization)
		Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
			/*! The [dim x 4] matrix of Bezier coefficients of this spline segment
			 * [P0 P0+T0/3 P1-T1/3 P1]
			 */
			Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
			double w =1.0;
			BzCf.col(0)<< this->source->get_P(), w;
			BzCf.col(1)<< this->source->get_P()+sourceT()/3.0*chordParametricLength(), w;
			BzCf.col(2)<< this->sink  ->get_P()-  sinkT()/3.0*chordParametricLength(), w;
			BzCf.col(3)<< this->sink  ->get_P(), w;
			return BzCf;
		}
		
		//change		
		//////////////////////////////////////////////////////////////
		//hermiteCoefficients: Hermite coefficients (uniform parametrization)
		Eigen::Matrix<double,dim,Ncoeff> hermiteCoefficients() const
        {
			Eigen::Matrix<double,dim,Ncoeff> HrCf;
			HrCf.col(0)= this->source->get_P();
			HrCf.col(1)= sourceT()*chordParametricLength();
			HrCf.col(2)= this->sink->get_P();
			HrCf.col(3)= sinkT()*chordParametricLength();
			return HrCf;
		}
		

		

		
        /************************************************************************/
		Eigen::Matrix<double,dim,Ncoeff> polynomialCoeff() const
        {/*! The matrix of coefficients of the polynomial associated to this 
			 *  SplineSegmentBase. If C=polynomialCoeff() then the polynomial is:
			 *  P(u)=C.col(0)+u*C.col(1)+u^2*C.col(2)+...
			 */
			return Coeff2Hermite<pOrder>::template h2c<dim>(hermiteCoefficients());
		}
		
		
		

		
#include "SplineSegmentBase_common.h"
		
	};
    
    //static data 
    template <typename Derived, short unsigned int dim>
	double SplineSegmentBase<Derived,dim,1>::alpha=0.5;

	
	/************************************************************************/
	/* SplineSegmentBase, template specialization corder=2 ******************/
	/************************************************************************/
	template <typename Derived, short unsigned int dim>
	class SplineSegmentBase<Derived,dim,2> : public NetworkLink<Derived>,
	/*	                                        */ public ParametricCurve<Derived,dim> {
		
		
		enum {corder=2};
#include<model/Geometry/Splines/SplineEnums.h>
		

        
        ////////////////////////////////////////////////////////////////
        //! Returns the length of the chord vector to the power alpha
        double chordParametricLength() const
        {
        	return std::pow(chordLength(),alpha);;
        }
		
	public:
		
        static double alpha;

		
		RowNcoeff get_UPOW(const double & uin){
			return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3), std::pow(uin,4), std::pow(uin,5)).finished();
		}
		
		
		RowNcoeffu get_UPOWu(const double & uin){
			return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2), 4.0*std::pow(uin,3), 5.0*std::pow(uin,4)).finished();
		}
		
		
		RowNcoeffuu get_UPOWuu(const double & uin){
			return (RowNcoeffuu()<<2.0, 6.0*uin, 12.0*std::pow(uin,2), 20.0*std::pow(uin,3)).finished();
		}
		

		//////////////////////////////////////////////////////////////
		//get_SFCH
		MatrixNcoeff get_SFCH() const {
			//! !!!!!!!! FINISH HERE !!!!!!!! //
			assert(0);
			return MatrixNcoeff::Zero();
		}
		
		MatrixNcoeffDim get_qH() const {
			return (MatrixNcoeffDim()<< this->source->get_P().transpose(),
					/*               */ sourceT().transpose(), 
					/*               */ this->source->get_K().transpose(),
					/*               */ this->  sink->get_P().transpose(),
					/*               */ sinkT().transpose(),
					/*               */ this->  sink->get_K().transpose()).finished();
		}
		
#include "SplineSegmentBase_common.h"
		
	};
    
    //static data
    template <typename Derived, short unsigned int dim>
	double SplineSegmentBase<Derived,dim,2>::alpha=0.5;
//	double SplineSegmentBase<Derived,dim,2>::alpha=1.0;

	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif

