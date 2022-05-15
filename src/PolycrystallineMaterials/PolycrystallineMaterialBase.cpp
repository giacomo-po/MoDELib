/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_PolycrystallineMaterialBase_cpp_
#define model_PolycrystallineMaterialBase_cpp_

#include <string>
#include <TerminalColors.h> // defines mode::cout
#include <TextFileParser.h>
#include <PolycrystallineMaterialBase.h>

namespace model
{
    
        
        /**********************************************************************/
        const std::string& PolycrystallineMaterialBase::getMaterialFile(const std::string& fileName)
        {
            std::cout<<greenBoldColor<<"Reading material file "<<fileName<<defaultColor<<std::endl;
            return fileName;
        }
        
        PolycrystallineMaterialBase::PolycrystallineMaterialBase(const std::string& fileName,const double& absoluteTemperature) :
        /* init */ materialFile(getMaterialFile(fileName))
        /* init */,materialName(TextFileParser(materialFile).readString("materialName",true))
        /* init */,crystalStructure(TextFileParser(materialFile).readString("crystalStructure",true))
        /* init */,T(absoluteTemperature)
        /* init */,Tm(TextFileParser(materialFile).readScalar<double>("Tm",true))
        /* init */,mu0_SI(TextFileParser(materialFile).readScalar<double>("mu0_SI",true))
        /* init */,mu1_SI(TextFileParser(materialFile).readScalar<double>("mu1_SI",true))
        /* init */,mu_SI(mu0_SI+mu1_SI*T)
        /* init */,nu(TextFileParser(materialFile).readScalar<double>("nu",true))
        /* init */,rho_SI(TextFileParser(materialFile).readScalar<double>("rho_SI",true))
        /* init */,cs_SI(sqrt(mu_SI/rho_SI))
        /* init */,b_SI(TextFileParser(materialFile).readScalar<double>("b_SI",true))
        /* init */,kB(kB_SI/mu_SI/std::pow(b_SI,3))
        /* init */,mu(1.0)
        /* init */,b(1.0)
        /* init */,cs(1.0)
        /* init */,dOmegav(TextFileParser(materialFile).readScalar<double>("dOmegav",true))
        /* init */,Ufv_SI(TextFileParser(materialFile).readScalar<double>("Ufv_eV",true) * eV2J)
        /* init */,Ufv(Ufv_SI/mu_SI/std::pow(b_SI,3))
        /* init */,Umv_SI(TextFileParser(materialFile).readScalar<double>("Umv_eV",true) * eV2J)
        /* init */,Umv(Umv_SI/mu_SI/std::pow(b_SI,3))
        /* init */,D0v_SI(TextFileParser(materialFile).readScalar<double>("D0v_SI",true))
        /* init */,Dv(D0v_SI/b_SI/cs_SI*exp(-Umv/kB/T))
        {
            
        }
        

}
#endif
