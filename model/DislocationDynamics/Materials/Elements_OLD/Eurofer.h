/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Element_Eurofer_H_
#define model_Element_Eurofer_H_

#include <deque>
#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/FCClattice.h>
#include <model/DislocationDynamics/Materials/BCClattice.h>
#include <model/DislocationDynamics/Materials/MaterialBase.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>
#include <model/IO/EigenDataReader.h>
namespace model
{
    using Eigen::Vector3d;
    
    template <>
    struct PeriodicElement<100,Isotropic>
    {
       // static const EigenDataReader EDR;
       // static const double Tempreal;
        typedef BCClattice<3> CrystalStructure;
        static constexpr int Z=100;
        static constexpr const char* name="Eurofer";
        static constexpr double nu=0.291;                // Poisson ratio [-]
      // static constexpr double mu=84.11e9;               // Shear modulus [Pa] 300K
        static constexpr double mu=78.536e9;               // Shear modulus [Pa] 573K
       // static constexpr double mu=-19.907e6*Tempreal+89942.76e6;               // Shear modulus [Pa]
        static constexpr double b=0.248e-9;            // Burgers vector[m]
        static constexpr double rho=7865.8;            // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        static constexpr double Tm=1181.0;              // melting temperature [K]
        
        
        static const DislocationMobility<BCClattice<3>> dm;
        
    };
    
               constexpr int    PeriodicElement<100,Isotropic>::Z;
        constexpr const char*   PeriodicElement<100,Isotropic>::name;
        constexpr double PeriodicElement<100,Isotropic>::nu;               // Poisson ratio [-]
        constexpr double PeriodicElement<100,Isotropic>::mu;              // Shear modulus [Pa]
        constexpr double PeriodicElement<100,Isotropic>::b;            // Burgers vector[m]
        constexpr double PeriodicElement<100,Isotropic>::rho;             // Mass density [kg/m^3]
        constexpr double PeriodicElement<100,Isotropic>::cs;        // Shear wave speed [m/s]
        constexpr double PeriodicElement<100,Isotropic>::Tm;              // melting temperature [K]

    
    //EDR.readScalarInFile("./DDInput.txt","temperature",Tempreal); // temperature   
    const DislocationMobility<BCClattice<3>> PeriodicElement<100,Isotropic>::dm=DislocationMobility<BCClattice<3>>(PeriodicElement<100,Isotropic>::b,
                                                                                              PeriodicElement<100,Isotropic>::mu,
                                                                                              PeriodicElement<100,Isotropic>::cs,
						                                              1.05e-04,0.0, // B0e [Pa*s], B1e [Pa*s/K]
						                                              9.8e-4,0.0,        // B0s [Pa*s], B1s [Pa*s/K]
						                                              8.3e-05,           // Bk [Pa*s]
						                                              0.86,              // dH0 [eV]
						                                              0.71,1.85,         // p,q
						                                              350.0,            // T0
						                                              542e6,            // tauC [Pa]
						                                              1.0,0.61,0.23,0.55,1.0 // non-Schmid coefficients
                                                                                              );
    
} // namespace model
#endif
