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

#ifndef _CoulombEnergy_h
#define _CoulombEnergy_h

#include <iostream>
#include <tutorials/ChargedParticles/ChargedParticle.h>

struct CoulombEnergy{
    
    typedef ChargedParticle::PositionType PositionType;
    
    
    /*****************************************/
    CoulombEnergy(ChargedParticle& cp1,ChargedParticle& cp2)
    {/*! Constructor with two ChargedParticle(s). Compu
      */
        
        //std::cout<<"Computing CoulombEnergy"<<std::endl;
        
        const PositionType R(cp1.P()-cp2.P());
        const double r(R.norm()); // distance between particles
        if(r!=0.0)
        {
            const double e(cp1.q()*cp2.q()/r); // the force on particle 2
            cp2.energy+=e;
            cp1.energy+=e;
        }
        
    }
    
    /*****************************************/
    static void reset(ChargedParticle& cp)
    {/*! Set energy of cp to zero
      */
        cp.energy=0.0;
    }
    
};


#endif

