/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SplineShapeFunction_H_
#define model_SplineShapeFunction_H_

#include <math.h>
#include <Eigen/Dense>


namespace model
{
    
    /************************************************************************/
    /* SplineShapeFunction, general case **************************************/
    /************************************************************************/
    template <short unsigned int dim, short unsigned int corder>
    class SplineShapeFunction {};
    
    
    /************************************************************************/
    /* SplineShapeFunction, template specialization corder=1 ******************/
    /************************************************************************/
    template <short unsigned int dim>
    struct SplineShapeFunction<dim,1>
    {
        
        static constexpr int corder = 1;
        static constexpr int Ncoeff= 2*(corder+1);
        typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
        typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
        typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
        
        /**********************************************************************/
        static RowNcoeff powers(const double & uin)
        {
            return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3)).finished();
        }
        
        static RowNcoeffu powersDiff1(const double & uin)
        {
            return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2)).finished();
        }
        
        static RowNcoeffuu powersDiff2(const double & uin)
        {
            return (RowNcoeffuu()<<2.0, 6.0*uin).finished();
        }
        
        
        /**********************************************************************/
        static MatrixNcoeff sfCoeffs(const double& g)
        {/*!\returns The matrix of shape function coefficients in Hermite form of this
          *  spline segment.
          */
            /*                         P0      T0    P1   T1  */
            return (MatrixNcoeff()<<  1.0,    0.0,  0.0, 0.0,             // u^0
                    /*            */  0.0,      g,  0.0, 0.0,             // u^1
                    /*            */ -3.0, -2.0*g,  3.0,  -g,             // u^2
                    /*            */  2.0,      g, -2.0,   g).finished(); // u^3;
        }
        
        /******************************************************************************/
        static RowNcoeff sf(const double& uin,const double& g)
        {
            return powers(uin)*sfCoeffs(g);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff1(const double& uin,const double& g)
        {
            return powersDiff1(uin)*sfCoeffs(g).template block<Ncoeff-1,Ncoeff>(1,0);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff2(const double& uin,const double& g)
        {
            return  powersDiff2(uin)*sfCoeffs(g).template block<Ncoeff-2,Ncoeff>(2,0);
        }
        
        //
        //        /**********************************************************************/
        //        Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
        //            /*!\returns The [dim x 4] matrix of Bezier coefficients of this spline segment
        //             * [P0 P0+T0/3 P1-T1/3 P1]
        //             */
        //            Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
        //            double w =1.0;
        //            BzCf.col(0)<< this->source->get_P(), w;
        //            BzCf.col(1)<< this->source->get_P()+sourceT()/3.0*chordParametricLength(), w;
        //            BzCf.col(2)<< this->sink  ->get_P()-  sinkT()/3.0*chordParametricLength(), w;
        //            BzCf.col(3)<< this->sink  ->get_P(), w;
        //            return BzCf;
        //        }
        //
        //        /**********************************************************************/
        //        Eigen::Matrix<double,dim,Ncoeff> hermiteCoefficients() const
        //        {
        //            Eigen::Matrix<double,dim,Ncoeff> HrCf;
        //            HrCf.col(0)= this->source->get_P();
        //            HrCf.col(1)= sourceT()*chordParametricLength();
        //            HrCf.col(2)= this->sink->get_P();
        //            HrCf.col(3)= sinkT()*chordParametricLength();
        //            return HrCf;
        //        }
        //
        //        /************************************************************************/
        //        Eigen::Matrix<double,dim,Ncoeff> polynomialCoeff() const
        //        {/*!\returns The matrix of coefficients of the polynomial associated to this
        //          *  SplineShapeFunction. If C=polynomialCoeff() then the polynomial is:
        //          *  P(u)=C.col(0)+u*C.col(1)+u^2*C.col(2)+...
        //          */
        //            return Coeff2Hermite<pOrder>::template h2c<dim>(hermiteCoefficients());
        //        }
        
    };
    
    //    //static data
    //    template <typename Derived, short unsigned int dim>
    //    double SplineShapeFunction<Derived,dim,1>::alpha=0.5;
    
    
    //	/************************************************************************/
    //	/* SplineShapeFunction, template specialization corder=0 ******************/
    //	/************************************************************************/
    //	template <typename Derived, short unsigned int dim>
    //	class SplineShapeFunction<Derived,dim,0> : public NetworkLink<Derived>,
    //    /*	                              */ public ParametricCurve<Derived, dim>
    //    /*	                              */ public ParametricCurve<Derived, dim>
    //    {
    //
    //
    //		enum {corder=0};
    //#include<model/Geometry/Splines/SplineEnums.h>
    //
    //
    //        ////////////////////////////////////////////////////////////////
    //        //! Returns the length of the chord vector to the power alpha
    //        double chordParametricLength() const
    //        {
    //        	return chordLength();
    //        }
    //
    //	public:
    //		RowNcoeff powers(const double & uin) const {
    //			return (RowNcoeff()<<1.0, uin).finished();
    //		}
    //
    //		//////////////////////////////////////////////////////////////
    //
    //		RowNcoeffu powersu(const double & uin) const {
    //			return (RowNcoeffu()<<1.0).finished();
    //		}
    //
    //
    //
    //		RowNcoeffuu powersuu(const double & uin) const {
    //			// ????????????????????????
    //		}
    //
    //		static MatrixNcoeff get_SFCH()
    //        {
    //			/*! The matrix of shape function coefficients in Hermite form of this
    //			 *  spline segment.
    //			 */
    //			/*                         P0   P1 */
    //			return (MatrixNcoeff()<<  1.0, 0.0,             // u^0
    //					/*            */ -1.0, 1.0).finished(); // u^1;
    //		}
    //
    //		MatrixNcoeffDim get_qH() const {
    //			return (MatrixNcoeffDim()<< this->source->get_P().transpose(),
    //					/*               */ this->  sink->get_P().transpose()).finished();
    //		}
    //
    //
    //
    //		//////////////////////////////////////////////////////////////
    //		//BezierCoefficients: Bezier coefficients (uniform parametrization)
    //		Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
    //			Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
    //			double w =1.0;
    //			BzCf.col(0)<< this->source->get_P(), w;
    //			BzCf.col(2)<< this->sink  ->get_P(), w;
    //			return BzCf;
    //		}
    //
    //#include "SplineShapeFunction_common.h"
    //
    //	};
    
    
    
    
    
    
    
    
    
    
    //	/************************************************************************/
    //	/* SplineShapeFunction, template specialization corder=2 ******************/
    //	/************************************************************************/
    //	template <typename Derived, short unsigned int dim>
    //	class SplineShapeFunction<Derived,dim,2> : public NetworkLink<Derived>,
    //	/*	                                        */ public ParametricCurve<Derived,dim> {
    //
    //
    //		enum {corder=2};
    //#include<model/Geometry/Splines/SplineEnums.h>
    //
    //
    //
    //        ////////////////////////////////////////////////////////////////
    //        //! Returns the length of the chord vector to the power alpha
    //        double chordParametricLength() const
    //        {
    //        	return std::pow(chordLength(),alpha);;
    //        }
    //
    //	public:
    //
    //        static double alpha;
    //
    //
    //		RowNcoeff powers(const double & uin){
    //			return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3), std::pow(uin,4), std::pow(uin,5)).finished();
    //		}
    //
    //
    //		RowNcoeffu powersu(const double & uin){
    //			return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2), 4.0*std::pow(uin,3), 5.0*std::pow(uin,4)).finished();
    //		}
    //
    //
    //		RowNcoeffuu powersuu(const double & uin){
    //			return (RowNcoeffuu()<<2.0, 6.0*uin, 12.0*std::pow(uin,2), 20.0*std::pow(uin,3)).finished();
    //		}
    //
    //
    //		//////////////////////////////////////////////////////////////
    //		//get_SFCH
    //		MatrixNcoeff get_SFCH() const {
    //			//! !!!!!!!! FINISH HERE !!!!!!!! //
    //			assert(0);
    //			return MatrixNcoeff::Zero();
    //		}
    //		
    //		MatrixNcoeffDim get_qH() const {
    //			return (MatrixNcoeffDim()<< this->source->get_P().transpose(),
    //					/*               */ sourceT().transpose(), 
    //					/*               */ this->source->get_K().transpose(),
    //					/*               */ this->  sink->get_P().transpose(),
    //					/*               */ sinkT().transpose(),
    //					/*               */ this->  sink->get_K().transpose()).finished();
    //		}
    //		
    //#include "SplineShapeFunction_common.h"
    //		
    //	};
    //    
    //    //static data
    //    template <typename Derived, short unsigned int dim>
    //	double SplineShapeFunction<Derived,dim,2>::alpha=0.5;
    ////	double SplineShapeFunction<Derived,dim,2>::alpha=1.0;
    //
    //	
    //	//////////////////////////////////////////////////////////////s
    //} // namespace model
}
#endif

