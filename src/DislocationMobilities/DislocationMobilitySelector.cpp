/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilitySelector_cpp_
#define _model_DislocationMobilitySelector_cpp_

#include <DislocationMobilitySelector.h>

namespace model
{

        DislocationMobilitySelector::DislocationMobilitySelector(const std::string& defaultStr_in):
        /* init */ defaultStr(defaultStr_in)
        {
            
        }


        std::shared_ptr<DislocationMobilityBase> DislocationMobilitySelector::getMobility(const std::string& dislocationMobilityType,
                                                                                          const PolycrystallineMaterialBase& material) const
        {
            
            if(dislocationMobilityType=="FCC")
            {
                return std::shared_ptr<DislocationMobilityBase>(new DislocationMobilityFCC(material));
            }
            else if(dislocationMobilityType=="BCC")
            {
                return std::shared_ptr<DislocationMobilityBase>(new DislocationMobilityBCC(material));
            }
            else if(dislocationMobilityType=="HEXbasal")
            {
                return std::shared_ptr<DislocationMobilityBase>(new DislocationMobilityHEXbasal(material));
            }
            else if(dislocationMobilityType=="HEXprismatic")
            {
                return std::shared_ptr<DislocationMobilityBase>(new DislocationMobilityHEXprismatic(material));
            }
            else if(dislocationMobilityType=="HEXpyramidal")
            {
                return std::shared_ptr<DislocationMobilityBase>(new DislocationMobilityHEXpyramidal(material));
            }
            else
            {
                return getMobility(defaultStr,material);
            }
        }
}
#endif
