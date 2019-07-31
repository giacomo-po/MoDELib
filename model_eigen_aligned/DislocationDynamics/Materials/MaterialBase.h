/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MaterialBase_H_
#define model_MaterialBase_H_

#include <string>
#include <TerminalColors.h> // defines mode::cout
#include <TextFileParser.h>

//#include <BCClattice.h>
//#include <FCClattice.h>
//#include <HEXlattice.h>

namespace model
{
    
    struct MaterialBase
    {
        
        static constexpr double kB_SI=1.38064852e-23; // Boltzmann constant [J/K]
        static constexpr double eV2J=1.6021766208e-19;  // Convert [eV] to [J]
        const std::string materialFile;
        const std::string materialName;
        
        const std::string crystalStructure;
        
        // Material constants in SI units
        const double T;         // Absolute temparature [K]
        const double Tm;        // Melting temparature [K]
        const double mu0_SI;    // [Pa]
        const double mu1_SI;    // [Pa/K]
        const double mu_SI;     // temperature-dependent shear modulus mu=m0+m1*T [Pa]
        const double nu;        // Poisson's ratio
        const double rho_SI;    // mass density [Kg/m^3]
        const double cs_SI;     // shear wave speed [m/s]
        const double b_SI;      // Burgers vector [m]
        
        // Material constants in code units
        const double kB;        // Boltzmann constant [-]
        const double mu;        // shear modulus [-]
        const double b;         // Burgers vector [-]
        const double cs;        // shear wave speed [-]
        
        
        // Vacancy concentration
        const double Omega;     // shear wave speed [-]
        const double Ufv_SI;
        const double Ufv;
        const double DVv;
        const double Udv_SI;        // shear wave speed [-]
        const double Udv;        // shear wave speed [-]
        const double D0v_SI;        // shear wave speed [-]
        const double Dv;        // shear wave speed [-]
//        const double c0;        // shear wave speed [-]

        
        /**********************************************************************/
        static const std::string& getMaterialFile(const std::string& fileName)
        {
            model::cout<<greenBoldColor<<"Creating material from file: "<<fileName<<defaultColor<<std::endl;
            return fileName;
        }
        
        /**********************************************************************/
        static double atomicVolume(const std::string& structure)
        {
            
            if(structure=="BCC")
            {
//                return BCClattice<3>::getLatticeBasis().determinant();
            }
            else if(structure=="FCC")
            {
 //               return FCClattice<3>::getLatticeBasis().determinant();
            }
            else if(structure=="HEX")
            {
//                return HEXlattice<3>::getLatticeBasis().determinant();
            }
            else
            {
                std::cout<<"Unknown crystal structure '"<<structure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
                return 0.0;
            }
            
            std::cout<<" !!!!!!!!! FINISH CALCULATION OF Omega !!!!!"<<std::endl;

            return 1.0;
        }
        
        /**************************************************************************/
        MaterialBase(const std::string& fileName) :
        /* init */ materialFile(getMaterialFile(fileName))
        /* init */,materialName(TextFileParser(materialFile).readString("materialName",true))
        /* init */,crystalStructure(TextFileParser(materialFile).readString("crystalStructure",true))
        /* init */,T(TextFileParser(materialFile).readScalar<double>("T",true))
        /* init */,Tm(TextFileParser(materialFile).readScalar<double>("Tm",true))
        /* init */,mu0_SI(TextFileParser(materialFile).readScalar<double>("mu0_SI",true))
        /* init */,mu1_SI(TextFileParser(materialFile).readScalar<double>("mu1_SI",true))
        /* init */,mu_SI(mu0_SI+mu1_SI*T)
        /* init */,nu(TextFileParser(materialFile).readScalar<double>("nu",true))
        /* init */,rho_SI(TextFileParser(materialFile).readScalar<double>("rho_SI",true))
        /* init */,cs_SI(sqrt(mu_SI/rho_SI))
        /* init */,b_SI(TextFileParser(materialFile).readScalar<double>("b_SI",true))
//        /* init */,Dv_SI(TextFileParser(materialFile).readScalar<double>("Dv_SI",true))
        /* init */,kB(kB_SI/mu_SI/std::pow(b_SI,3))
        /* init */,mu(1.0)
        /* init */,b(1.0)
        /* init */,cs(1.0)
        ///* init */,Omega(atomicVolume(crystalStructure))
        /* init */,Omega(TextFileParser(materialFile).readScalar<double>("Omega_SI",true) * 1e-30 /std::pow(b_SI,3))
        /* init */,Ufv_SI(TextFileParser(materialFile).readScalar<double>("Ufv_eV",true) * eV2J)
        /* init */,Ufv(Ufv_SI/mu_SI/std::pow(b_SI,3))
        /* init */,DVv(TextFileParser(materialFile).readScalar<double>("DVv",true))
        /* init */,Udv_SI(TextFileParser(materialFile).readScalar<double>("Udv_eV",true) * eV2J)
        /* init */,Udv(Udv_SI/mu_SI/std::pow(b_SI,3))
        /* init */,D0v_SI(TextFileParser(materialFile).readScalar<double>("D0v_SI",true))
        /* init */,Dv(D0v_SI/b_SI/cs_SI*exp(-Udv/kB/T))
        {
//            model::cout<<greenBoldColor<<"Reading material file: "<<materialFile<<defaultColor<<std::endl;

            std::cout<<"FINISH CALCULATION OF ATOMIC VOLUME AND REMOVE IT FROM INPUT FILE"<<std::endl;
            
        }
        
    };
}
#endif
