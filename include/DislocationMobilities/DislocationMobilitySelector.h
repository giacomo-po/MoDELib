/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilitySelector_h_
#define _model_DislocationMobilitySelector_h_

#include <DislocationMobilityBase.h>
#include <DislocationMobilityFCC.h>
#include <DislocationMobilityBCC.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXpyramidal.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct DislocationMobilitySelector
    {
        
        const std::string defaultStr;
        
        DislocationMobilitySelector(const std::string& defaultStr_in);
                
        std::shared_ptr<DislocationMobilityBase> getMobility(const std::string& dislocationMobilityType,
                                                             const PolycrystallineMaterialBase& material) const;
        
    };
    
    
}
#endif
