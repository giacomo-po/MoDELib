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
#include <model/DislocationDynamics/Materials/CrystalStructures.h>


namespace model {

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
    struct PeriodicElement<13,Isotropic>
    {
        enum{Z=13};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static constexpr double nu=0.347;               // Poisson ratio [-]

        static constexpr double mu=26.0e9;              // Shear modulus [Pa]
        static constexpr double b=0.2851e-9;            // Burgers vector[m]
        static constexpr double rho=2700.0;             // Mass density [kg/m^3]
        
        static constexpr double tauP=0.0;               // Peierls stress [Pa]
        static constexpr double p=1.0;                  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double q=1.69;                 // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double Ae=3.3333e-07;           // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double As=3.3333e-07;           // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ta=0.0;  // Athermal transition temperature [K]

        static const Eigen::Matrix<double,4,2> dH0;     // activation energy prefactor for kink nucleation [J]
    };
    
    const std::string PeriodicElement<13,Isotropic>::name="Aluminum";
    const Eigen::Matrix<double,4,2> PeriodicElement<13,Isotropic>::dH0=(Eigen::Matrix<double,4,2>()<<
                                                                        0.0,0.0, 
                                                                        0.0,0.0,
                                                                        0.0,0.0,
                                                                        0.0,0.0// screw, edge for <110> planes
                                                                        ).finished()*1.602e-19;
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<28,Isotropic>
    {
        enum{Z=28};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static constexpr double nu=0.31;                // Poisson ratio [-]
        static constexpr double mu=76e9;                // Shear modulus [Pa]
        static constexpr double b=0.2489e-9;            // Burgers vector[m]
        static constexpr double rho=8908.0;             // Mass density [kg/m^3]
        
        static constexpr double tauP=0.0;               // Peierls stress [Pa]
        static constexpr double p=1.0;                  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double q=1.69;                 // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double Ae=3.3333e-07;           // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double As=3.3333e-07;           // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ta=0.0;  // Athermal transition temperature [K]

        static const Eigen::Matrix<double,4,2> dH0;     // activation energy prefactor for kink nucleation [J]
    };
    
    const std::string PeriodicElement<28,Isotropic>::name="Nickel";
    const Eigen::Matrix<double,4,2> PeriodicElement<28,Isotropic>::dH0=(Eigen::Matrix<double,4,2>()<<
                                                                        0.0,0.0,
                                                                        0.0,0.0,
                                                                        0.0,0.0,
                                                                        0.0,0.0// screw, edge for <110> planes
                                                                        ).finished()*1.602e-19; 

    /**************************************************************************/
    /**************************************************************************/    
    template <>
    struct PeriodicElement<29,Isotropic>
    {
        enum{Z=29};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static constexpr double nu=0.34;                // Poisson ratio [-]
        static constexpr double mu=48e9;                // Shear modulus [Pa]
        static constexpr double b=0.2556e-9;            // Burgers vector[m]
        static constexpr double rho=8940.0;             // Mass density [kg/m^3]
        
        
        static constexpr double tauP=0.0; // Peierls stress [Pa]
        static constexpr double p=1.0;  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double q=1.69;  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double Ae=3.3333e-07;  // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double As=3.3333e-07;  // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ta=0.0;  // Athermal transition temperature [K]

        static const Eigen::Matrix<double,4,2> dH0; // activation energy prefactor for kink nucleation [J]
    };
    
    const std::string PeriodicElement<29,Isotropic>::name="Copper";
    const Eigen::Matrix<double,4,2> PeriodicElement<29,Isotropic>::dH0=(Eigen::Matrix<double,4,2>()<<
                                                                        0.0,0.0,
                                                                        0.0,0.0,
                                                                        0.0,0.0,
                                                                        0.0,0.0 // screw, edge for <111> plane
                                                                        ).finished()*1.602e-19;
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<26,Isotropic>
    {
        enum{Z=26};
        typedef BCC CrystalStructure;
        
        static const std::string name;
        static constexpr double nu=0.29;                // Poisson ratio [-]
        static constexpr double mu=82e9;                // Shear modulus [Pa]
        static constexpr double b=0.2482e-9;            // Burgers vector[m]
        static constexpr double rho=7874.0;             // Mass density [kg/m^3]
        
        static constexpr double tauP=420.0e6; // Peierls stress [Pa]
        static constexpr double p=1.0;  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double q=1.69;  // exponent in (1-(T/Ta)^p)^q [-]
        static constexpr double Ae=1.0e-6;  // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double As=1.0e-6;  // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ta=400.0;  // Athermal transition temperature [K]

        static const Eigen::Matrix<double,18,2> dH0; // activation energy prefactor for kink nucleation [J]
    };
    
    const std::string PeriodicElement<26,Isotropic>::name="Iron";
    const Eigen::Matrix<double,18,2> PeriodicElement<26,Isotropic>::dH0=(Eigen::Matrix<double,18,2>()<<
                                                                        0.8,0.08, // (0,1,1)
                                                                        0.8,0.08, // (1,0,1)
                                                                        0.8,0.08, // (1,-1,0)
                                                                        0.8,0.08, // (0,1,-1)
                                                                        0.8,0.08, // (1,1,0)
                                                                        0.8,0.08, // (1,0,-1)
                                                                        2.8,2.8, // (2,-1,1)
                                                                        2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8 // (2,-1,1)
                                                                        ).finished()*1.602e-19;
    
    /**************************************************************************/
    /**************************************************************************/    
    template <>
    struct PeriodicElement<74,Isotropic>
    {
        enum{Z=74};
        typedef BCC CrystalStructure;
        
        static const std::string name;
        static constexpr double nu=0.28;                // Poisson ratio [-]
        static constexpr double mu=161e9;                // Shear modulus [Pa]
        static constexpr double b=0.2722e-9;            // Burgers vector[m]
        static constexpr double rho=19250.0;    // Mass density [kg/m^3]
        
        static constexpr double tauP=910.0e6;   // Peierls stress [Pa]
        static constexpr double p=0.86;          //
        static constexpr double q=1.69;         //
//        static constexpr double A=1.0e-6;     // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ae=3.3333e-07;   // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double As=3.3333e-07;   // coefficient of v0=tau*b/(A*T) [Pa*s/K]
        static constexpr double Ta=800.0;       // Athermal transition temperature [K]

        static const Eigen::Matrix<double,18,2> dH0; // activation energy prefactor for kink nucleation [J]
    };
    
    const std::string PeriodicElement<74,Isotropic>::name="Tungsten";
    const Eigen::Matrix<double,18,2> PeriodicElement<74,Isotropic>::dH0=(Eigen::Matrix<double,18,2>()<<
                                                                         1.63,0.0, // (0,1,1)
                                                                         1.63,0.0, // (0,1,1)
                                                                         1.63,0.0, // (0,1,1)
                                                                         1.63,0.0, // (0,1,1)
                                                                         1.63,0.0, // (0,1,1)
                                                                         1.63,0.0, // (0,1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8, // (2,-1,1)
                                                                         2.8,2.8 // (2,-1,1)
                                                                         ).finished()*1.602e-19;
    
    

    
    /**************************************************************************/
} // namespace model
#endif
