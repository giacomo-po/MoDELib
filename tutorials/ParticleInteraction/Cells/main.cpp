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


// run
// mpiexec -np 4 particles

#include <model/ParticleInteraction/ParticleSystem.h> // the main object from pil library
#include <model/Utilities/TerminalColors.h> // the main object from pil library

#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem
#include <tutorials/ParticleInteraction/ChargedParticles/CoulombForce.h> // a user-defined type of force interaction between ChargedParticle objects
#include <tutorials/ParticleInteraction/ChargedParticles/CoulombEnergy.h> // a user-defined type of energy interaction between ChargedParticle objects


using namespace model;

int main (int argc, char * argv[]) {
    
    // 0- define the type of ParticleSystem specifying the type of particles
    typedef model::ParticleSystem<ChargedParticle> ChargedParticleSystem;
    
    // 1- create a particleSystem of ChargedParticle(s)
    double cellSize=6.0;
    ChargedParticleSystem particleSystem(argc,argv,cellSize);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // 2- add some ChargedParticle(s) to the particleSystem with random position
    //    and different charges
    std::cout<<model::blueBoldColor<<"CREATING RANDOM INITIAL PARTICLES"<<model::defaultColor<<std::endl;
    typedef  ChargedParticleSystem::PositionType PositionType; // helper
    for (int k=0;k<500;++k)
    {
        particleSystem.addParticle(PositionType::Random()*10.0, 1.0);
    }
    // Add a more particles 
    for (int k=0;k<500;++k)
    {
//        particleSystem.addParticle(PositionType::Random()*3.0+PositionType::Ones()*2.0, 1.0); // cluster
                particleSystem.addParticle(PositionType::Random()*10.0, 1.0); // uniform
    }
    
    
    
    // 4- compute all binary CoulombEnergy interactions
    //particleSystem.computeInteraction<CoulombEnergy>();
    
    // 5- Output all particles to a file
    //std::cout<<pil::blueBoldColor<<"Outputting particles to P-file and cells to C-file."<<pil::defaultColor<<std::endl;
//    pil::SequentialOutputFile<'P',1> pFile;
//    pFile<<particleSystem<<std::endl;

    particleSystem.MPIoutput();

    
	particleSystem.partionCells();
    
    // 3- compute all binary CoulombForce interactions
    particleSystem.computeInteraction<CoulombForce>();


    particleSystem.getInteractionResult<CoulombForce>(0);

    
    particleSystem.MPIoutput();
    
    return 0;
}


