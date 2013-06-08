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


#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/Utilities/TerminalColors.h> 

#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem

//#include <tutorials/ParticleInteraction/ChargedParticles/CoulombEnergy.h> // a user-defined type of energy interaction between ChargedParticle objects


using namespace model;

int main (int argc, char * argv[]) {
    
    // 0- define the type of ParticleSystem specifying the type of particles
    typedef model::ParticleSystem<ChargedParticle> ChargedParticleSystem;
    
    
//    std::cout<<" bytes in ChargedParticle="<<sizeof(ChargedParticle)<<std::endl;
    
    // 1- create a particleSystem of ChargedParticle(s)
    double cellSize=1.0;
    ChargedParticleSystem particleSystem(argc,argv,cellSize);
    
    // 2- add some ChargedParticle(s) to the particleSystem with random position
    //    in [-10,10] and charge=1.0
    std::cout<<model::blueBoldColor<<"CREATING RANDOM INITIAL PARTICLES"<<model::defaultColor<<std::endl;
    typedef  ChargedParticleSystem::PositionType PositionType; // helper
    for (size_t k=0;k<500000;++k)
    {
        // note that PositionType::Random() returns values in [-1, 1]
        particleSystem.addParticle(PositionType::Random()*10.0, 1.0);
    }
    // Add a more particles 
//    for (int k=0;k<500;++k)
//    {
////        particleSystem.addParticle(PositionType::Random()*3.0+PositionType::Ones()*2.0, 1.0); // cluster
//                particleSystem.addParticle(PositionType::Random()*10.0, 1.0); // uniform
//    }
    
    
    


//
//
    
//    std::cout<<model::blueBoldColor<<"PARTITIONING"<<model::defaultColor<<std::endl;
	//particleSystem.partionSystem();
    
    
//    SequentialOutputFile<'P',true> pFile0;
//    pFile0<<particleSystem.particles()<<std::endl;

//    particleSystem.MPIoutput();

    std::cout<<model::blueBoldColor<<"COMPUTING INTERACTION"<<model::defaultColor<<std::endl;

    //    // 3- compute all binary CoulombForce interactions
    typedef ChargedParticle::CoulombForceInteraction CoulombForceInteraction;
    particleSystem.computeInteraction<CoulombForceInteraction>();
    
    
    SequentialOutputFile<'P',true> pFile1;
    pFile1<<particleSystem.particles()<<std::endl;


    ////
////
////    particleSystem.getInteractionResult<CoulombForce>(0);
////
////    
//    particleSystem.MPIoutput();
    
    return 0;
}


