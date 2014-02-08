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
#include <vector>
#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/Material.h>

namespace model {
    
    
    template <short unsigned int dim>  // non-type template
    class DislocationMobility
    {
        
//        const double k_b ;  // boltzman's constant [ev/K]
        
//        enum{Z=74};
        
        
    public:
        
        typedef Eigen::Matrix<double, dim, 1> VectorDim;
//        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        
//        const double& dH0;
        const VectorDim& glidePlaneNormal;
        const VectorDim& Burgers;
        const Eigen::Matrix<double,1,2> dH0;
        const double& tauP; // Peierls stress [Pa]
//        const double& Ta;
        const double& p;
        const double& q;
        const double& A;
        
        
        /**********************************************************************/
        DislocationMobility(const VectorDim& pN,const VectorDim& b) :
//        /* init list */ k_b(8.617e-5),
        /* init list */ glidePlaneNormal(pN),
        /* init list */ Burgers(b),
//        /* init list */ dH0(PeriodicElement<Z,Isotropic>::dH0),
        /* init list       */ dH0(Material<Isotropic>::dH0.row(CrystalOrientation<dim>::planeID(glidePlaneNormal))),
        /* init list */ tauP(Material<Isotropic>::tauP),
//        /* init list */ Ta(Material<Isotropic>::Ta),
        /* init list */ p(Material<Isotropic>::p),
        /* init list */ q(Material<Isotropic>::q),
        /* init list */ A(Material<Isotropic>::A)
        {/*! Contructor initializes data members.
          */
            
//            std::cout<<"planeNormal="<<glidePlaneNormal.transpose()<<"\n"
//            <<", dH0="<<dH0<<"\n"
//            <<", k_b*T="<<Material<Isotropic>::kT*Material<Isotropic>::T<<"\n"
//            <<std::endl;
            
        }
        
        /**********************************************************************/
        VectorDim getVelocity(const VectorDim& pK,
//                           const MatrixDim& sigma,
//                           const VectorDim& b,
                           const VectorDim& t
//                           const VectorDim& n,
//                           const double& T
                           ) const
        {/*!@param[in] sigma the stress tensor
          * @param[in] b Burgers vector
          * @param[in] t tangent vector
          * @param[in] n normal vector
          * @param[in] T temperature
          \returns The dislocation velocity vector.
          */
            
            VectorDim v(VectorDim::Zero()); // initialize velocity to 0-vector
            
            const double pKn(pK.norm()); // norm of pK
            
            if(pKn>0.0 && Material<Isotropic>::T>0.0)
            {
                v = pK/pKn; // unit vector in the direction of pK
                
                // 1- multiply by pre-exponential factor
                v *= 1.0-exp(-sqrt(pKn/(A*Material<Isotropic>::T)));
                // v *= pKn/(A*Material<Isotropic>::T);
                // 2- multiply by exponential factor
                // 2.1- compute activation energy
                const double aE(activationEnergy(pKn));
                if (aE>0.0)
                {
//                    if (T>0.0)
//                    {
                        v *= exp(-aE/(Material<Isotropic>::kT*Material<Isotropic>::T));
//                    }
//                    else // T==0
//                    {
//                        v *= 0.0;
//                    }
                }
//                else
//                {
//                    v *= 0.0;
//                }

                
//                const double v0(1.0-exp(pKn/(A*Material<Isotropic>::T)));
                
                
//                const double trss(std::fabs((sigma*b.normalized()).dot(n.normalized()))); //resolved shear stress
//                const double tau_star = transitionStress(T);
                //            assert(tau_star >= tau_o && "tau_star must be >= tau_o"); // will make sure transition stress is always greater than unlocking stress, just an internal check
                
                
//                const double v0(trss*PeriodicElement<Z,Isotropic>::b/(A*Material<Isotropic>::T));
                
                
//                const double cs(sqrt(PeriodicElement<Z,Isotropic>::mu/PeriodicElement<Z,Isotropic>::rho));
                
//                double v(cs*(1.0-exp(-v0/cs)));
                
                //            if(trss < tau_star)
                //            {
                
            }

            return v;
        }
        
        
        /**********************************************************************/
        double activationEnergy(const double& trss) const
        {/*!
          */
//            const double tau0(transitionStress(Material<Isotropic>::T));
            
            const double T0(1000.0);
            const double tauS(tauP*std::pow(1.0-Material<Isotropic>::T/T0,2));
            
            const double DH(dH0(0));
            
//            return (trss>Material<Isotropic>::tauS)? 0.0 : dH0(0)*std::pow(1.0 - std::pow(trss/Material<Isotropic>::tauS,1.0/q),1.0/p);
            return (trss>tauS)? 0.0 : DH*std::pow(1.0 - std::pow(trss/tauS,1.0/q),1.0/p);
        }
        
//        /**********************************************************************/
//        double transitionStress(const double& T) const
//        {/*!@param[in] T input temperature
//          *\returns the thermal transition stress at temperature T, in units of
//          * tauP.
//          *
//          * The transition stress is computed according to:
//          *	\f[
//          *		\sigma(T)=\left[1-\left(\frac{T}{T_a}\right)^p\right]^q
//          *	\f]
//          */
//            return (T>Ta)? 0.0 : tauP*std::pow(1.0 - std::pow(T/Ta,p),q);
//        }
        
        //        /**********************************************************************/
        //        double getPreFactor(const MatrixDim& sigma, const double& T)
        //        {// returns the pre-factor as function of current stress tensor "sigma" and temp T
        //
        //            double preFac = 1.0;
        //            
        //            return preFac;
        //        }
        
        
        
    };
    
}

#endif


