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
#include <Eigen/Geometry>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/DislocationDynamics/Materials/Material.h>

namespace model
{
    
    
    template <short unsigned int dim>  // non-type template
    class DislocationMobility
    {
        
//        const double k_b ;  // boltzman's constant [ev/K]
        
//        enum{Z=74};
        
        
    public:
        
        typedef Eigen::Matrix<double, dim, 1> VectorDim;
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        
//        const double& dH0;
        const VectorDim& glidePlaneNormal;
        const VectorDim& Burgers;
        const Eigen::Matrix<double,1,2> dH0;
        const double& tauP; // Peierls stress [Pa]
//        const double& Ta;
        const double& p;
        const double& q;
        const double Be;
        const double Bs;
        const double& Ta; // Peierls stress [Pa]
        
        
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
        /* init list */ Be(Material<Isotropic>::Be),
        /* init list */ Bs(Material<Isotropic>::Bs),
        /* init list */ Ta(Material<Isotropic>::Ta)
        {/*! Contructor initializes data members.
          */
            
//            std::cout<<"planeNormal="<<glidePlaneNormal.transpose()<<"\n"
//            <<", dH0="<<dH0<<"\n"
//            <<", k_b*T="<<Material<Isotropic>::kT*Material<Isotropic>::T<<"\n"
//            <<std::endl;
            
        }
        
        /**********************************************************************/
        VectorDim getVelocity(
                              //const VectorDim& pK,
                           const MatrixDim& sigma,
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
            
//            const double pKn(pK.norm()); // norm of pK
            
            const VectorDim pK=(sigma*Burgers).cross(t);
            const double pKnorm(pK.norm());
            
            if(pKnorm>0.0 && Material<Isotropic>::T>0.0)
            {
                v = pK/pKnorm; // unit vector in the direction of pK
                
                // 1- multiply by pre-exponential factor
                //v *= 1.0-exp(-sqrt(pKn/(A*Material<Isotropic>::T)));
                
                const double cos2Theta=std::pow(t.dot(Burgers.normalized()),2);
                
                const double vE=pKnorm/Be;
                      double vS=pKnorm/Bs;


                
                // v *= pKn/(A*Material<Isotropic>::T);
                // 2- multiply by exponential factor
                // 2.1- compute activation energy
                const double dG(deltaG(sigma,t));
                if (dG>0.0)
                {
                        vS *= exp(-dG/(Material<Isotropic>::kb*Material<Isotropic>::T));
                }
                
                v *= vS*cos2Theta+vE*(1.0-cos2Theta);

                
            }
            

            return v;
        }
        
        
        /**********************************************************************/
        double deltaG(const MatrixDim& sigma, const VectorDim& t) const
        {/*!
          */

//            const double trss = (sigma*Burgers.normalized()).dot(glidePlaneNormal);
            const double trss = std::fabs((sigma*Burgers.normalized()).dot(glidePlaneNormal));

            const VectorDim m = glidePlaneNormal.cross(Burgers.normalized());
            const double s=(sigma*Burgers.normalized()).dot(m);
            
            
            double thetaTwin=M_PI/3.0;
            
            const double BdotT = Burgers.dot(t);
            
            if(BdotT<0.0)
            {// A "negative" dislocation
                thetaTwin*=-1.0;
            }
            
            //const MatrixDim R = ;
            const VectorDim n1 = Eigen::AngleAxisd(thetaTwin, Burgers.normalized())*glidePlaneNormal;
            const double trss1 = (sigma*Burgers.normalized()).dot(n1);

            
            const VectorDim m1 = n1.cross(Burgers.normalized());
            const double s1=(sigma*Burgers.normalized()).dot(m1);

            
            const double a1=0.0;
            const double a2=0.56;
            const double a3=0.75;
            
            const double tauS = std::fabs(trss+a1*trss1+a2*s+a3*s1)/tauP;
            
            return dH0(0)*(std::pow(1.0-std::pow(tauS,p),q)-Material<Isotropic>::T/Ta);
            
        }
        
        
        
        
    };
    
}

#endif


