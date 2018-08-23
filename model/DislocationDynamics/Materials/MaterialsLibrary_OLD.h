/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_MaterialsLibrary_H_
#define model_MaterialsLibrary_H_

#include <cmath>
#include <string>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/TerminalColors.h> // defines mode::cout



namespace model
{
    

    
    /**************************************************************************/
    /**************************************************************************/
//    template<>
    struct MaterialsLibrary
    {
        
        
//        static int selectedMaterial;
        
//        static const PeriodicElement<13,Isotropic> Al;
//        static const PeriodicElement<28,Isotropic> Ni;
//        static const PeriodicElement<29,Isotropic> Cu;
//        static const PeriodicElement<74,Isotropic>  W;
        
        
        static const Copper cu;

        static std::map<std::string,const MaterialBase* const> getMaterials()
        {
            std::map<std::string,const MaterialBase* const> temp;
            temp.emplace( cu.name(),&MaterialsLibrary::cu);
            
            std::cout<<"Materials Library:"<<std::endl;
            for(const auto& material : temp)
            {
//                FINISH HERE
                std::cout<<"    "<<material.second->name()<<std::endl;
            }
            
            
            return temp;
        }
        
    public:
        
        static const std::map<std::string,const MaterialBase* const> materialsMap;

        
    };
    
    const Copper MaterialsLibrary::cu;
    
//    const PeriodicElement<13,Isotropic> MaterialsLibrary::Al=PeriodicElement<13,Isotropic>();
//    const PeriodicElement<28,Isotropic> MaterialsLibrary::Ni=PeriodicElement<28,Isotropic>();
//    const PeriodicElement<29,Isotropic> MaterialsLibrary::Cu=PeriodicElement<29,Isotropic>();
//    const PeriodicElement<74,Isotropic> MaterialsLibrary::W =PeriodicElement<74,Isotropic>();
    const std::map<std::string,const MaterialBase* const> MaterialsLibrary::materialsMap=MaterialsLibrary::getMaterials();

    
}
#endif
