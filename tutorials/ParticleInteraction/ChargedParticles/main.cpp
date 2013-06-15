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

//#include <model/MPI/MPIcout.h>

//#include <tutorials/ParticleInteraction/ChargedParticles/CoulombEnergy.h> // a user-defined type of energy interaction between ChargedParticle objects



using namespace model;

int main (int argc, char * argv[]) {
    
    // 0- define the type of ParticleSystem specifying the type of particles
    typedef model::ParticleSystem<ChargedParticle> ChargedParticleSystem;
    
    // 1- create a particleSystem of ChargedParticle(s)
    ChargedParticleSystem particleSystem(argc,argv); // initialized constructor: for both serial and parallel
    //ChargedParticleSystem particleSystem; // uninitialized constructor: for both serial and parallel
    //ChargedParticleSystem::initMPI(argc,argv); // if uninitialized constructor is used, call initMPI (parallel only)
    
    double cellSize=1.0;
    particleSystem.setCellSize(cellSize);
    
    // 2- add some ChargedParticle(s) to the particleSystem with random position
    //    in [-10,10] and charge=1.0
    std::cout<<"Creating particles..."<<std::endl;
    typedef typename ChargedParticleSystem::PositionType PositionType; // helper
    for (size_t k=0;k<500000;++k)
    {
        // note that PositionType::Random() returns values in [-1, 1]
        particleSystem.addParticle(PositionType::Random()*10+PositionType::Ones()*0.0*cellSize, 1.0);
    }
//    std::cout<<model::greenColor<<" done."<<model::defaultColor<<std::endl;
//    std::cout<<model::greenColor<<" done."<<model::defaultColor<<std::endl;

    // Add a more particles
//    for (int k=0;k<500;++k)
//    {
////        particleSystem.addParticle(PositionType::Random()*3.0+PositionType::Ones()*2.0, 1.0); // cluster
//                particleSystem.addParticle(PositionType::Random()*10.0, 1.0); // uniform
//    }
    
    
    std::cout<<"There are "<<particleSystem.cells().size()<<" cells"<<std::endl;


//
//
    
//    std::cout<<model::blueBoldColor<<"PARTITIONING"<<model::defaultColor<<std::endl;
	//particleSystem.partionSystem();
    
    
//    SequentialOutputFile<'P',true> pFile0;
//    pFile0<<particleSystem.particles()<<std::endl;

//    particleSystem.MPIoutput();
    
//    std::cout<<model::blueBoldColor<<"COMPUTING INTERACTION"<<model::defaultColor<<std::endl;

    
    // -3 computation of the CoulombForceInteraction
    typedef typename ChargedParticle::CoulombForceInteraction CoulombForceInteraction;
    // -3.1  reset CoulombForceInteraction
 //   particleSystem.resetInteraction<CoulombForceInteraction>();

    
//    particleSystem.computeMoment0<CoulombForceInteraction>();

    
    // -3.2a compute all binary CoulombForce interactions
    std::cout<<"Computing nearest-neighbor interaction..."<<std::endl;
    particleSystem.computeNeighborInteraction<CoulombForceInteraction>();
//    std::cout<<model::greenColor<<" done."<<model::defaultColor<<std::endl;

    // -3.2b
//    typedef typename ChargedParticle::TotalCellCharge CellCharge;
//    particleSystem.computeCellProperty<CellCharge>();
//    particleSystem.computeFarInteraction<CoulombForceInteraction>();

    // -4 output
    std::cout<<"Writing output file P/P_0.txt..."<<std::endl;
    SequentialOutputFile<'P',true> pFile1;
    pFile1<<particleSystem.particles()<<std::endl;
//    std::cout<<model::greenColor<<" done."<<model::defaultColor<<std::endl;


    ////
////
////    particleSystem.getInteractionResult<CoulombForce>(0);
////
////    
//    particleSystem.MPIoutput();
    
    return 0;
}


