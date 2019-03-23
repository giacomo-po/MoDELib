/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Element_W_H_
#define model_Element_W_H_

#include <deque>
#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/FCClattice.h>
#include <model/DislocationDynamics/Materials/BCClattice.h>
#include <model/DislocationDynamics/Materials/MaterialBase.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>

namespace model
{
    using Eigen::Vector3d;
    
    template <>
    struct PeriodicElement<74,Isotropic>
    {
        typedef BCClattice<3> CrystalStructure;
        static constexpr int Z=74;
        static constexpr const char* name="Tungsten";
        static constexpr double nu=0.28;                // Poisson ratio [-]
        static constexpr double mu=161e9;               // Shear modulus [Pa]
        static constexpr double b=0.2722e-9;            // Burgers vector[m]
        static constexpr double rho=19250.0;            // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        static constexpr double Tm=3695.0;              // melting temperature [K]
        
        static const DislocationMobility<BCClattice<3>> dm;
        
    };
    
               constexpr int    PeriodicElement<74,Isotropic>::Z;
        constexpr const char*   PeriodicElement<74,Isotropic>::name;
        constexpr double PeriodicElement<74,Isotropic>::nu;               // Poisson ratio [-]
        constexpr double PeriodicElement<74,Isotropic>::mu;              // Shear modulus [Pa]
        constexpr double PeriodicElement<74,Isotropic>::b;            // Burgers vector[m]
        constexpr double PeriodicElement<74,Isotropic>::rho;             // Mass density [kg/m^3]
        constexpr double PeriodicElement<74,Isotropic>::cs;        // Shear wave speed [m/s]
        constexpr double PeriodicElement<74,Isotropic>::Tm;              // melting temperature [K]

    
    const DislocationMobility<BCClattice<3>> PeriodicElement<74,Isotropic>::dm=DislocationMobility<BCClattice<3>>(PeriodicElement<74,Isotropic>::b,
                                                                                              PeriodicElement<74,Isotropic>::mu,
                                                                                              PeriodicElement<74,Isotropic>::cs,
                                                                                              4.26e-04,0.87e-06, // B0e [Pa*s], B1e [Pa*s/K]
                                                                                              9.8e-4,0.0,        // B0s [Pa*s], B1s [Pa*s/K]
                                                                                              8.3e-05,           // Bk [Pa*s]
                                                                                              1.63,              // dH0 [eV]
                                                                                              0.86,1.69,         // p,q
                                                                                              0.8*Tm,            // T0
                                                                                              2.03e9,            // tauC [Pa]
                                                                                              1.2943,1.1702,4.9087,9.3352,0.3107 // non-Schmid coefficients
                                                                                              );
    
} // namespace model
#endif
