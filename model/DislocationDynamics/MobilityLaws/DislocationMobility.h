/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobility_h_
#define _model_DislocationMobility_h_

#include <iostream>
#include <cmath>
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/HexLattice.h>

#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/TerminalColors.h> // defines mode::cout

namespace model
{
    

    
    template <typename CrystalStructure>
    struct DislocationMobility
    {
        //        static_assert(false, "CrystalStructure must be FCC or BCC");
    };



    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct DislocationMobility<Hexagonal>
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        const double B0;
        const double B1;
        
        /**********************************************************************/
        constexpr DislocationMobility(const double& b_real,
                                      const double& mu_real,
                                      const double& cs_real,
                                      const double& B1e_real,
                                      const double& B1s_real) :
        
        /* init */ B0(0.0),
        /* init */ B1(0.5*(B1e_real+B1s_real)*cs_real/(mu_real*b_real))
        {/*! Empty constructor is required by constexpr
          */
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& , // xi
                        const VectorDim& n,
                        const double& T) const
        {
            return std::fabs(b.transpose()*S*n)/(B0+B1*T);
        }
        
        
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct DislocationMobility<FCC>
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        const double B0;
        const double B1;
        
        /**********************************************************************/
        constexpr DislocationMobility(const double& b_real,
                                      const double& mu_real,
                                      const double& cs_real,
                                      const double& B1e_real,
                                      const double& B1s_real) :
        
        /* init */ B0(0.0),
        /* init */ B1(0.5*(B1e_real+B1s_real)*cs_real/(mu_real*b_real))
        {/*! Empty constructor is required by constexpr
          */
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& , // xi
                        const VectorDim& n,
                        const double& T) const
        {
            return std::fabs(b.transpose()*S*n)/(B0+B1*T);
        }
        
        
        
    };
    
    //    /**************************************************************************/
    //    /**************************************************************************/
    //    template <>
    //    struct DislocationMobility<FCC>
    //    {
    //
    //        typedef Eigen::Matrix<double,3,3> MatrixDim;
    //        typedef Eigen::Matrix<double,3,1> VectorDim;
    //
    //        const double B0;
    //        const double B1e;
    //        const double B1s;
    //
    //        /**********************************************************************/
    //        constexpr DislocationMobility(const double& b_real,
    //                                      const double& mu_real,
    //                                      const double& cs_real,
    //                                      const double& B1e_real,
    //                                      const double& B1s_real) :
    //
    //        /* init */ B0(0.0),
    //        /* init */ B1e(B1e_real*cs_real/(mu_real*b_real)),
    //        /* init */ B1s(B1s_real*cs_real/(mu_real*b_real))
    //        {/*! Empty constructor is required by constexpr
    //          */
    //        }
    //
    //        /**********************************************************************/
    //        double velocity(const MatrixDim& S,
    //                        const VectorDim& b,
    //                        const VectorDim& xi, // tangent vector
    //                        const VectorDim& n,
    //                        const double& T) const
    //        {
    //            const double ve=std::fabs(b.transpose()*S*n)/(B1e*T);
    //            const double vs=std::fabs(b.transpose()*S*n)/(B1s*T);
    //            const double cos2=std::pow(b.normalized().dot(xi),2);
    //            const double sin2=1.0-cos2;
    //            return vs*cos2+ve*sin2;
    //        }
    //        
    //    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct DislocationMobility<BCC>
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        //! Boltzmann constant in [eV]
        static constexpr double kB=8.617e-5;
        
        const double h;
        const double w;
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double Bk;
        const double dH0;
        const double p;
        const double q;
        const double T0;
        const double tauC;
        //        const std::array<double,5> A;
        const double a0;
        const double a1;
        const double a2;
        const double a3;
        const double a4;
        //        const double a5;
        
        
        
        /**********************************************************************/
        constexpr DislocationMobility(const double& b_real,
                                      const double& mu_real,
                                      const double& cs_real,
                                      const double& B0e_real, const double& B1e_real,
                                      const double& B0s_real, const double& B1s_real,
                                      const double& Bk_real,
                                      const double& dH0_real,
                                      const double& p_in,
                                      const double& q_in,
                                      const double& T0_in,
                                      const double& tauC_in,
                                      const double& a0_in,
                                      const double& a1_in,
                                      const double& a2_in,
                                      const double& a3_in,
                                      const double& a4_in) :
        //                                      const std::array<double,5>& A_in) :
        /* init */ h(2.0*sqrt(2.0)/3.0), // units of b
        /* init */ w(25.0), // units of b
        /* init */ B0e(B0e_real*cs_real/(mu_real*b_real)),
        /* init */ B1e(B1e_real*cs_real/(mu_real*b_real)),
        /* init */ B0s(B0s_real*cs_real/(mu_real*b_real)),
        /* init */ B1s(B1s_real*cs_real/(mu_real*b_real)),
        /* init */ Bk(  Bk_real*cs_real/(mu_real*b_real)),
        /* init */ dH0(dH0_real),
        /* init */ p(p_in),
        /* init */ q(q_in),
        /* init */ T0(T0_in),
        /* init */ tauC(tauC_in/mu_real),
        //        /* init */ A(A_in)
        /* init */ a0(a0_in),
        /* init */ a1(a1_in),
        /* init */ a2(a2_in),
        /* init */ a3(a3_in),
        /* init */ a4(a4_in)
        
        {/*! Empty constructor is required by constexpr
          */
        }
        
        static double sigmoid(const double & x)
        {
            return 1.0/(1.0+exp(-x));
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T) const
        {
            
            
            const double bNorm=b.norm();
            const VectorDim s = b/bNorm;
            const VectorDim n1 = Eigen::AngleAxisd(M_PI/3.0,s)*n;
            
            //            // Compute components of non-Schmid model
            const double tau=std::fabs(s.transpose()*S*n); // magnitude of resolved shear stress
            const double tauOrt=n.cross(s).transpose()*S*n;
            const double tau1=std::fabs(s.transpose()*S*n1); // resolved schear stress on
            const double tauOrt1=n1.cross(s).transpose()*S*n1;
            
            const double num=tau+a1*tau1;
            //
                        assert(num>=0.0 && "num must be >= 0.");
//            const double den=(1.0+0.5*A[1])*tauC*(A[4]+A[0]*sigmoid(-(A[2]*tauOrt+A[3]*tauOrt1)/tauC));
            const double den=(1.0+0.5*a1)*tauC*(a4+a0*sigmoid(-(a2*tauOrt+a3*tauOrt1)/tauC));
            assert(den>0.0 && "den must be > 0.");

            const double Theta=num/den;
            const double dg = (Theta<1.0)? (std::pow(1.0-std::pow(Theta,p),q)-T/T0) : 0.0;
            const double dg1 = (dg>0.0)? dg : 0.0;
            const double expCoeff = exp(-dH0*dg1/(2.0*kB*T));

            // Compute screw drag coeff
            const double sgm=sigmoid((0.05-dg1)/0.05);
            const double Bs=Bk*w/(2.0*h)*(1.0-sgm)+(B0s+B1s*T)*sgm;
            
            // Compute screw velocity
            const double vs=tau*bNorm/Bs*expCoeff;
            
            // Compute edge velocity
            const double ve=tau*bNorm/(B0e+B1e*T);
            
            // Interpolate ve and vs
            const double cos2=std::pow(b.normalized().dot(xi),2);
            const double sin2=1.0-cos2;
            return vs*cos2+ve*sin2;
        }
        
    };
    
}

#endif


