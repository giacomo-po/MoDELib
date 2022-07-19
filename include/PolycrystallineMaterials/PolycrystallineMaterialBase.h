/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_PolycrystallineMaterialBase_H_
#define model_PolycrystallineMaterialBase_H_

#include <string>
#include <numbers>
#include <TerminalColors.h> // defines mode::cout
#include <TextFileParser.h>

namespace model
{
    
    struct PolycrystallineMaterialBase
    {
        
        static constexpr double kB_SI=1.38064852e-23; // Boltzmann constant [J/K]
        static constexpr double eV2J=1.6021766208e-19;  // Convert [eV] to [J]
        const std::string materialFile;
        const std::string materialName;
        
        const std::string crystalStructure;
        
        // Material constants in SI units
        const double T;         // simulation temparature [K]
        const double Tm;        // Melting temparature [K]
        const double mu0_SI;    // [Pa]
        const double mu1_SI;    // [Pa/K]
        const double mu_SI;     // temperature-dependent shear modulus mu=mu0+mu1*T [Pa]
        const double nu;        // Poisson's ratio
        const double rho_SI;    // mass density [Kg/m^3]
        const double cs_SI;     // shear wave speed [m/s]
        const double b_SI;      // Burgers vector [m]
        
        // Material constants in code units
        const double kB;        // Boltzmann constant [-]
        const double mu;        // shear modulus [-]
        const double b;         // Burgers vector [-]
        const double cs;        // shear wave speed [-]
        const double C1;        // 1-nu
        const double C2;        // 1.0/(4.0*M_PI*C1)
        const double C3;        // 1.0-2.0*nu;
        const double C4;        // 0.5*C2;

        
        const double dOmegav;
        const double Ufv_SI;
        const double Ufv;
        const double Umv_SI;     // vacancy migration energy [eV]
        const double Umv;        // vacancy migration energy [-]
        const double D0v_SI;        // shear wave speed [-]
        const double Dv;        // shear wave speed [-]
        
        static const std::string& getMaterialFile(const std::string& fileName);
        PolycrystallineMaterialBase(const std::string& fileName,const double& absoluteTemperature);
        
    };
}
#endif
