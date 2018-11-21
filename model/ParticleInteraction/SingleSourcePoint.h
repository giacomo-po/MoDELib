/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_SingleSourcePoint_h
#define _model_SingleSourcePoint_h

#include <SpatialCellParticle.h>


namespace model
{

    
    
    // template specialization: one or more fields
    template<typename Derived,
    /*    */ typename FieldType>
    struct SingleSourcePoint
    {
        
        bool enabled;

        SingleSourcePoint(const bool& enb) :
        /* base init */ enabled(enb)
        {
            
        }
    };
    
    
    
} // end namespace
#endif
