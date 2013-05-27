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

#include <iostream>
#include <iomanip>

#include <pil/ParticleSystem.h> // the main object from pil library
#include <pil/Utilities/TerminalColors.h> // the main object from pil library

#include <ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem
#include <CoulombForce.h> // a user-defined type of force interaction between ChargedParticle objects
#include <CoulombEnergy.h> // a user-defined type of energy interaction between ChargedParticle objects


int main (int argc, char * const argv[]) {
    
    // 0- define the type of ParticleSystem specifying the type of particles
    typedef pil::ParticleSystem<ChargedParticle> ChargedParticleSystem;
    
    // 1- create a particleSystem of ChargedParticle(s)
    std::cout<<pil::blueBoldColor<<"CREATING INITIAL PARTICLES"<<pil::defaultColor<<std::endl;
    ChargedParticleSystem particleSystem;
    
    // 2- add some ChargedParticle(s) to the particleSystem with random position
    //    and different charges
    typedef ChargedParticleSystem::PositionType PositionType; // helper
    particleSystem.addParticle(PositionType::Random(), 1.0);
    particleSystem.addParticle(PositionType::Random(),-1.0);
    particleSystem.addParticle(PositionType::Random(), 2.0);
    
    // 3- compute all binary CoulombForce interactions
    particleSystem.computeInteraction<CoulombForce>();
    
    // 4- compute all binary CoulombEnergy interactions
    particleSystem.computeInteraction<CoulombEnergy>();
    
    // 5- List all particles again
    std::cout<<pil::blueBoldColor<<"PARTICLE SYSTEM AFTER COMPUTATION:"<<pil::defaultColor<<std::endl;
    std::cout<<particleSystem<<std::endl;

    
    // 6- Reset all computation for next iteration
    particleSystem.resetInteraction<CoulombForce>();
    particleSystem.resetInteraction<CoulombEnergy>();
    
    // 7- List all particles again
    std::cout<<pil::blueBoldColor<<"PARTICLE SYSTEM AFTER RESET:"<<pil::defaultColor<<std::endl;
    std::cout<<particleSystem<<std::endl;
	
    return 0;
}


