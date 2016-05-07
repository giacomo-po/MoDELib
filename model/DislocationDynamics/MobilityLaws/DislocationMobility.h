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
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
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
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct DislocationMobility<BCC>
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double Bk;
        const double DH0;
        const double p;
        const double q;
        const double T0;
        
        /**********************************************************************/
        constexpr DislocationMobility(const double& b_real,
                                      const double& mu_real,
                                      const double& cs_real,
                                      const double& B0e_real, const double& B1e_real,
                                      const double& B0s_real, const double& B1s_real,
                                      const double& Bk_real,
                                      const double& DH0_real,
                                      const double& p_in,
                                      const double& q_in,
                                      const double& T0_in) :
        /* init */ B0e(B0e_real*cs_real/(mu_real*b_real)),
        /* init */ B1e(B1e_real*cs_real/(mu_real*b_real)),
        /* init */ B0s(B0s_real*cs_real/(mu_real*b_real)),
        /* init */ B1s(B1s_real*cs_real/(mu_real*b_real)),
        /* init */ Bk(Bk_real*cs_real/(mu_real*b_real)),
        /* init */ DH0(DH0_real*0),
        /* init */ p(p_in),
        /* init */ q(q_in),
        /* init */ T0(T0_in)
        {/*! Empty constructor is required by constexpr
          */
        }
        
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                                  const VectorDim& b,
                                  const VectorDim& xi,
                                  const VectorDim& n,
                                  const double& T) const
        {
            
            // magnitude of resolved shear stress
            const double taub=std::fabs(b.transpose()*S*n);
            
            // Compute edge velocity
            const double ve=taub/(B0e*B1e*T);

            // Compute screw velocity
            const double vs=ve/50.0;

            
            // Interpolate
            const double cos2=std::pow(b.normalized().dot(xi),2);
            const double sin2=1.0-cos2;

            return vs*cos2+ve*sin2;
        }
        
    };
    
}

#endif



//    template <short unsigned int dim>  // non-type template
//    class DislocationMobility
//    {
//
////        const double k_b ;  // boltzman's constant [ev/K]
//
////        enum{Z=74};
//
//
//    public:
//
//        typedef Eigen::Matrix<double, dim, 1> VectorDim;
//        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
//
////        const double& dH0;
//        const VectorDim& glidePlaneNormal;
//        const VectorDim& Burgers;
//        const Eigen::Matrix<double,1,2> dH0;
//        const double& tauP; // Peierls stress [Pa]
////        const double& Ta;
//        const double& p;
//        const double& q;
//        const double Be;
//        const double Bs;
//        const double& Ta; // Peierls stress [Pa]
//
//
//        /**********************************************************************/
//        DislocationMobility(const VectorDim& pN,const VectorDim& b) :
////        /* init list */ k_b(8.617e-5),
//        /* init list */ glidePlaneNormal(pN),
//        /* init list */ Burgers(b),
////        /* init list */ dH0(PeriodicElement<Z,Isotropic>::dH0),
//        /* init list       */ dH0(Material<Isotropic>::dH0.row(CrystalOrientation<dim>::planeID(glidePlaneNormal))),
//        /* init list */ tauP(Material<Isotropic>::tauP),
////        /* init list */ Ta(Material<Isotropic>::Ta),
//        /* init list */ p(Material<Isotropic>::p),
//        /* init list */ q(Material<Isotropic>::q),
//        /* init list */ Be(Material<Isotropic>::Be),
//        /* init list */ Bs(Material<Isotropic>::Bs),
//        /* init list */ Ta(Material<Isotropic>::Ta)
//        {/*! Contructor initializes data members.
//          */
//
////            std::cout<<"planeNormal="<<glidePlaneNormal.transpose()<<"\n"
////            <<", dH0="<<dH0<<"\n"
////            <<", k_b*T="<<Material<Isotropic>::kT*Material<Isotropic>::T<<"\n"
////            <<std::endl;
//
//        }
//
//        /**********************************************************************/
//        VectorDim getVelocity(
//                              //const VectorDim& pK,
//                           const MatrixDim& sigma,
////                           const VectorDim& b,
//                           const VectorDim& t
////                           const VectorDim& n,
////                           const double& T
//                           ) const
//        {/*!@param[in] sigma the stress tensor
//          * @param[in] b Burgers vector
//          * @param[in] t tangent vector
//          * @param[in] n normal vector
//          * @param[in] T temperature
//          \returns The dislocation velocity vector.
//          */
//
//            VectorDim v(VectorDim::Zero()); // initialize velocity to 0-vector
//
////            const double pKn(pK.norm()); // norm of pK
//
//            const VectorDim pK=(sigma*Burgers).cross(t);
//            const double pKnorm(pK.norm());
//
//            if(pKnorm>0.0 && Material<Isotropic>::T>0.0)
//            {
//                v = pK/pKnorm; // unit vector in the direction of pK
//
//                // 1- multiply by pre-exponential factor
//                //v *= 1.0-exp(-sqrt(pKn/(A*Material<Isotropic>::T)));
//
//                const double cos2Theta=std::pow(t.dot(Burgers.normalized()),2);
//
//                const double vE=pKnorm/Be;
//                      double vS=pKnorm/Bs;
//
//
//
//                // v *= pKn/(A*Material<Isotropic>::T);
//                // 2- multiply by exponential factor
//                // 2.1- compute activation energy
//                const double dG(deltaG(sigma,t));
//                if (dG>0.0)
//                {
//                        vS *= exp(-dG/(Material<Isotropic>::kb*Material<Isotropic>::T));
//                }
//
//                v *= vS*cos2Theta+vE*(1.0-cos2Theta);
//
//
//            }
//
//
//            return v;
//        }
//
//
//        /**********************************************************************/
//        double deltaG(const MatrixDim& sigma, const VectorDim& t) const
//        {/*!
//          */
//
////            const double trss = (sigma*Burgers.normalized()).dot(glidePlaneNormal);
//            const double trss = std::fabs((sigma*Burgers.normalized()).dot(glidePlaneNormal));
//
//            const VectorDim m = glidePlaneNormal.cross(Burgers.normalized());
//            const double s=(sigma*Burgers.normalized()).dot(m);
//
//
//            double thetaTwin=M_PI/3.0;
//
//            const double BdotT = Burgers.dot(t);
//
//            if(BdotT<0.0)
//            {// A "negative" dislocation
//                thetaTwin*=-1.0;
//            }
//
//            //const MatrixDim R = ;
//            const VectorDim n1 = Eigen::AngleAxisd(thetaTwin, Burgers.normalized())*glidePlaneNormal;
//            const double trss1 = (sigma*Burgers.normalized()).dot(n1);
//
//
//            const VectorDim m1 = n1.cross(Burgers.normalized());
//            const double s1=(sigma*Burgers.normalized()).dot(m1);
//
//
//            const double a1=0.0;
//            const double a2=0.56;
//            const double a3=0.75;
//
//            const double tauS = std::fabs(trss+a1*trss1+a2*s+a3*s1)/tauP;
//
//            return dH0(0)*(std::pow(1.0-std::pow(tauS,p),q)-Material<Isotropic>::T/Ta);
//
//        }
//
//
//
//
//    };

