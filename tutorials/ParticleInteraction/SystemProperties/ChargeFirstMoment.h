/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _ChargeFirstMoment_h
#define _ChargeFirstMoment_h

#include <../ChargedParticles/ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem


struct ChargeFirstMoment{
    
    
    
    
    typedef ChargedParticle::PositionType PropertyType;
    
    
    PropertyType chargeMoment;
        
    /*****************************************/
    ChargeFirstMoment() : chargeMoment(PropertyType::Zero())
    {/*! Constructor initializes ChargeFirstMoment to 0.0
      */
    }
    
    /*****************************************/
    operator const PropertyType& () const
    {/*! Cast operator returns ChargeFirstMoment
      */
        return chargeMoment;
    }

};


#endif

