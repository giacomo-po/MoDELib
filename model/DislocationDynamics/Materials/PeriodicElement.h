/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PERIODICELEMENT_H_
#define model_PERIODICELEMENT_H_

#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/CrystalStructures.h>


namespace model {

    /**************************************************************************/
    /**************************************************************************/
    template <int Z, typename SymmetryType>
    struct PeriodicElement {
       // static_assert(0,"PERIODIC ELEMENT NOT IMPLEMENTED");
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<13,Isotropic> {
    
        enum{Z=13};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
 
    };
    
    const std::string PeriodicElement<13,Isotropic>::name="Aluminum";
    const double PeriodicElement<13,Isotropic>::nu =0.347;        // Poisson ratio [-]
    const double PeriodicElement<13,Isotropic>::mu =26e9;         // Shear modulus [Pa]
    const double PeriodicElement<13,Isotropic>::b  =0.2851e-9;    // Burgers vector[m]
    const double PeriodicElement<13,Isotropic>::B  =1.0e-4;       // Dislocation drag coefficient [Pa*sec]
    const double PeriodicElement<13,Isotropic>::rho=2700.0;       // Mass density [kg/m^3]
	
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<28,Isotropic> {
        
        enum{Z=28};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
    
    const std::string PeriodicElement<28,Isotropic>::name="Nickel";
    const double PeriodicElement<28,Isotropic>::nu =0.31;       // Poisson ratio [-]
    const double PeriodicElement<28,Isotropic>::mu =76e9;       // Shear modulus [Pa]
    const double PeriodicElement<28,Isotropic>::b  =0.2489e-9;  // Burgers vector[m]
    const double PeriodicElement<28,Isotropic>::B  =1.0e-4;     // Dislocation drag coefficient [Pa*sec]
    const double PeriodicElement<28,Isotropic>::rho=8908.0;     // Mass density [kg/m^3]
    
    /**************************************************************************/
    /**************************************************************************/    
    template <>
    struct PeriodicElement<29,Isotropic> {
        
        enum{Z=29};
        typedef FCC CrystalStructure;
        
        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
    
    const std::string PeriodicElement<29,Isotropic>::name="Copper";
    const double PeriodicElement<29,Isotropic>::nu =0.34;    // Poisson ratio [-]
    const double PeriodicElement<29,Isotropic>::mu =48e9;    // Shear modulus [Pa]
    const double PeriodicElement<29,Isotropic>::b  =0.2556e-9; // Burgers vector[m]
    const double PeriodicElement<29,Isotropic>::B  =1.0e-4;  // Dislocation drag coefficient [Pa*sec]
    const double PeriodicElement<29,Isotropic>::rho=8940.0;  // Mass density [kg/m^3]
   
    /**************************************************************************/
    /**************************************************************************/    
    template <>
    struct PeriodicElement<74,Isotropic> {
        
        enum{Z=74};
        typedef BCC CrystalStructure;
        
        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
    
    const std::string PeriodicElement<74,Isotropic>::name="W";
    const double PeriodicElement<74,Isotropic>::nu =0.28;    // Poisson ratio [-]
    const double PeriodicElement<74,Isotropic>::mu =161e9;    // Shear modulus [Pa]
    const double PeriodicElement<74,Isotropic>::b  =0.2722e-9; // Burgers vector[m]
    const double PeriodicElement<74,Isotropic>::B  =1.0e-4;  // Dislocation drag coefficient [Pa*sec]
    const double PeriodicElement<74,Isotropic>::rho=19250.0;  // Mass density [kg/m^3]
    
    /**************************************************************************/
    /**************************************************************************/    
    template <>
    struct PeriodicElement<26,Isotropic> {
        
        enum{Z=26};
        typedef BCC CrystalStructure;
        
        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
    
    const std::string PeriodicElement<26,Isotropic>::name="Fe";
    const double PeriodicElement<26,Isotropic>::nu =0.29;    // Poisson ratio [-]
    const double PeriodicElement<26,Isotropic>::mu =82e9;    // Shear modulus [Pa]
    const double PeriodicElement<26,Isotropic>::b  =0.2482e-9; // Burgers vector[m]
    const double PeriodicElement<26,Isotropic>::B  =1.0e-4;  // Dislocation drag coefficient [Pa*sec]
    const double PeriodicElement<26,Isotropic>::rho=7874.0;  // Mass density [kg/m^3]

    /**************************************************************************/
} // namespace model
#endif
