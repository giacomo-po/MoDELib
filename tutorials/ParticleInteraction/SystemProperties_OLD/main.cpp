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

#include <pil/ParticleSystem.h> // the main object from pil library
#include <pil/SystemProperties.h> // the main object from pil library
#include <pil/Utilities/TerminalColors.h> // the main object from pil library

#include <../ChargedParticles/ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem
#include <TotalCharge.h>
#include <ChargeFirstMoment.h>

int main (int argc, char * const argv[]) {
    
    // 0- define additional system properties
    typedef pil::SystemProperties<TotalCharge,ChargeFirstMoment> ChargedSystemProperties;
    
    // 1- define the type of ParticleSystem specifying the type of particles
    typedef pil::ParticleSystem<ChargedParticle,ChargedSystemProperties> ChargedParticleSystem;
    
    // 2- create a particleSystem of ChargedParticle(s)
    std::cout<<pil::blueBoldColor<<"CREATING INITIAL PARTICLES"<<pil::defaultColor<<std::endl;
    ChargedParticleSystem particleSystem;

    // 3- add some ChargedParticle(s) to the particleSystem with random position and different charges
    typedef typename ChargedParticleSystem::PositionType PositionType; // helper
    particleSystem.addParticle(PositionType::Random(), 1.0);
    particleSystem.addParticle(PositionType::Random(),-1.0);
    particleSystem.addParticle(PositionType::Random(), 2.0);
    
    // ?- Get TotalCharge
    std::cout<<pil::blueBoldColor<<"TOTAL SYSTEM CHARGE"<<pil::defaultColor<<std::endl;
    std::cout<< particleSystem.getProperty<TotalCharge>() <<std::endl;

    // ?- Get ChargeFirstMoment
    std::cout<<pil::blueBoldColor<<"SYSTEM CHARGE FIRST MOMENT"<<pil::defaultColor<<std::endl;
    std::cout<< particleSystem.getProperty<ChargeFirstMoment>().transpose() <<std::endl;
	
    return 0;
}


