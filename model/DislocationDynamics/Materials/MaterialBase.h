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

namespace model
{
    
    struct MaterialBase
    {
        
        static constexpr double kB_SI=1.38064852e-23; // Boltzmann constant [J/K]
        const std::string materialFile;
        const std::string materialName;
        
        static const std::string& getMaterialFile(const std::string& fileName)
        {
            model::cout<<greenBoldColor<<"Reading material file: "<<fileName<<defaultColor<<std::endl;
            return fileName;
        }
        
        /**************************************************************************/
        MaterialBase(const std::string& fileName) :
        /* init */ materialFile(getMaterialFile(fileName)),
        /* init */ materialName(TextFileParser(materialFile).readString("materialName",true))
        {
//            model::cout<<greenBoldColor<<"Reading material file: "<<materialFile<<defaultColor<<std::endl;

        }
        
    };
}
#endif
