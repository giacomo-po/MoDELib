/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_InteractionBase_h
#define _model_InteractionBase_h

#include <vector>

namespace model {
    
    template<typename _DataType, int _DataPerParticle>
    struct InteractionBase
    {
        enum{DataPerParticle=_DataPerParticle};
        
        static std::vector<_DataType> resultVector;


        
    };
    
    
    // declare statica data
    template<typename _DataType, int _DataPerParticle>
    std::vector<_DataType> InteractionBase<_DataType,_DataPerParticle>::resultVector;


    
} // end namespace
#endif
