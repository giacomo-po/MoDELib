/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_Nickel_H_
#define model_Nickel_H_

#include <string>
#include <model/Dislocations/Materials/MaterialSymmetry.h>
#include <model/Dislocations/Materials/CrystalStructures.h>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename SymmetryType>
    struct Nickel { };
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct Nickel<Isotropic> {
        
        typedef FCC CrystalStructure;
        
        static const std::string name;
        
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
        
    const std::string Nickel<Isotropic>::name="Nickel";
    const double Nickel<Isotropic>::nu =0.31;       // Poisson ratio [-]
    const double Nickel<Isotropic>::mu =76e9;       // Shear modulus [Pa]
    const double Nickel<Isotropic>::b  =0.2489e-9;  // Burgers vector[m]
    const double Nickel<Isotropic>::B  =1.0e-4;     // Dislocation drag coefficient [Pa*sec]
    const double Nickel<Isotropic>::rho=8908.0;     // Mass density [kg/m^3]
	
    /**************************************************************************/
} // namespace model 
#endif
