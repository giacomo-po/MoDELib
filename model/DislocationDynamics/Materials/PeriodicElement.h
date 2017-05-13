/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicElement_H_
#define model_PeriodicElement_H_

#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/Materials/HexLattice.h>

#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int Z, typename SymmetryType>
    struct PeriodicElement
    {
        // static_assert(0,"PERIODIC ELEMENT NOT IMPLEMENTED");
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<12,Isotropic>
    {
        typedef Hexagonal CrystalStructure;
        static constexpr int Z=12;
        static constexpr auto name="Magnesium";
        static constexpr double nu=0.29;                // Poisson ratio [-]
        static constexpr double mu=17.0e9;              // Shear modulus [Pa]
        static constexpr double b=0.321e-9;             // Burgers vector (basal lattice vector "a") [m]
        static const double c2a;                        // c/a ratio [-]
        static constexpr double rho=1738.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        static const DislocationMobility<Hexagonal> dm;
        
    };
    
    const DislocationMobility<Hexagonal> PeriodicElement<12,Isotropic>::dm=DislocationMobility<Hexagonal>(b,mu,cs,3.9e-08,7.5e-08);
    const double PeriodicElement<12,Isotropic>::c2a=1.624;              // c/a ratio [-]

    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<13,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=13;
        static constexpr auto name="Aluminum";
        static constexpr double nu=0.347;               // Poisson ratio [-]
        static constexpr double mu=26.0e9;              // Shear modulus [Pa]
        static constexpr double b=0.2851e-9;            // Burgers vector[m]
        static constexpr double rho=2700.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        //! FCC-mobility law with data from Olmsted MSMSE 13(3), 2005.
//        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,3.9e-08,7.5e-08);
        static const DislocationMobility<FCC> dm;

    };
    
    const DislocationMobility<FCC> PeriodicElement<13,Isotropic>::dm=DislocationMobility<FCC>(b,mu,cs,3.9e-08,7.5e-08);

    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<28,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=28;
        static constexpr auto name="Nickel";
        static constexpr double nu=0.31;                // Poisson ratio [-]
        static constexpr double mu=76e9;                // Shear modulus [Pa]
        static constexpr double b=0.2489e-9;            // Burgers vector[m]
        static constexpr double rho=8908.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        //! FCC mobility law with data from Olmsted MSMSE 13(3), 2005.
//        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,5.0e-08,6.4e-08);
        static const DislocationMobility<FCC> dm;

    };
    
    const DislocationMobility<FCC> PeriodicElement<28,Isotropic>::dm=DislocationMobility<FCC>(b,mu,cs,5.0e-08,6.4e-08);
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<29,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=29;
        static constexpr auto name="Copper";
        static constexpr double nu=0.34;                // Poisson ratio [-]
        static constexpr double mu=48e9;                // Shear modulus [Pa]
        static constexpr double b=0.2556e-9;            // Burgers vector[m]
        static constexpr double rho=8940.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
//        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,3.3333e-07,3.3333e-07);
        static const DislocationMobility<FCC> dm;
        
    };
    
    const DislocationMobility<FCC> PeriodicElement<29,Isotropic>::dm=DislocationMobility<FCC>(b,mu,cs,3.3333e-07,3.3333e-07);
    
    
    template <>
    struct PeriodicElement<26,Isotropic>
    {
        typedef BCC CrystalStructure;
        static constexpr int Z=26;
        static constexpr auto name="Iron";
        static constexpr double nu=0.291;                // Poisson ratio [-]
        static constexpr double mu=86e9;               // Shear modulus [Pa]
        static constexpr double b=0.2482e-9;            // Burgers vector[m]
        static constexpr double rho=7865.8;            // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        static constexpr double Tm=1181.0;              // melting temperature [K]
        
        static const DislocationMobility<BCC> dm;
        
    };
    
    const DislocationMobility<BCC> PeriodicElement<26,Isotropic>::dm=DislocationMobility<BCC>(b,mu,cs,
                                                                                              1.05e-04,0.0, // B0e [Pa*s], B1e [Pa*s/K]
                                                                                              9.8e-4,0.0,        // B0s [Pa*s], B1s [Pa*s/K]
                                                                                              8.3e-05,           // Bk [Pa*s]
                                                                                              0.86,              // dH0 [eV]
                                                                                              0.71,1.85,         // p,q
                                                                                              350.0,            // T0
                                                                                              542e6,            // tauC [Pa]
                                                                                              1.0,0.61,0.23,0.55,1.0 // non-Schmid coefficients
                                                                                              );
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<74,Isotropic>
    {
        typedef BCC CrystalStructure;
        static constexpr int Z=74;
        static constexpr auto name="Tungsten";
        static constexpr double nu=0.28;                // Poisson ratio [-]
        static constexpr double mu=161e9;               // Shear modulus [Pa]
        static constexpr double b=0.2722e-9;            // Burgers vector[m]
        static constexpr double rho=19250.0;            // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        static constexpr double Tm=3695.0;              // melting temperature [K]
        
        static const DislocationMobility<BCC> dm;
    
    };
    
    const DislocationMobility<BCC> PeriodicElement<74,Isotropic>::dm=DislocationMobility<BCC>(b,mu,cs,
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
