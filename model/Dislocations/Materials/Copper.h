/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_COPPER_H_
#define model_COPPER_H_

#include <model/Dislocations/Materials/MaterialSymmetry.h>
#include <model/Dislocations/Materials/CrystalStructures.h>


namespace model {
    
    template <typename SymmetryType>
    struct Copper { };
    
    template <>
    struct Copper<Isotropic> {
        
        typedef FCC CrystalStructure;

        static const std::string name;
        static const double nu;
        static const double mu;
        static const double b;
        static const double B;
        static const double rho;
        
    };
    
    const std::string Copper<Isotropic>::name="Copper";
    const double Copper<Isotropic>::nu =0.34;    // Poisson ratio [-]
    const double Copper<Isotropic>::mu =48e9;    // Shear modulus [Pa]
    const double Copper<Isotropic>::b  =0.2556e-9; // Burgers vector[m]
    const double Copper<Isotropic>::B  =1.0e-4;  // Dislocation drag coefficient [Pa*sec]
    const double Copper<Isotropic>::rho=8940.0;  // Mass density [kg/m^3]
	
	//////////////////////////////////////////////////////////////
} // namespace model 
#endif

