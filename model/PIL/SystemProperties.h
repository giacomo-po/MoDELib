/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _SystemProperties_h
#define _SystemProperties_h

namespace pil {
//    template <typename ...UserProperties>
//    struct SystemProperties{
//        
//    
//    };
    
    
    // Varaidic Predeclaration
    template <typename ...AdditionalUserProperties>
    struct SystemProperties {};
    
    // Template Specialization: one ore more base types
    template <typename UserProperty, typename ...AdditionalUserProperties>
    struct SystemProperties<UserProperty,AdditionalUserProperties...> :
    /* inheritance */ public UserProperty,
    /* inheritance */ public SystemProperties<AdditionalUserProperties...>
    {
        
    };
    
    // Template Specialization: one base types
    template <typename UserProperty>
    struct SystemProperties<UserProperty> :
    /* inheritance */ public UserProperty
    {
        
    };
    
    // Template Specialization: zero base types
    template <>
    struct SystemProperties<>
    {
        
    };
}

#endif

