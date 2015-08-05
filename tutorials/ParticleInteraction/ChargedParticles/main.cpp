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


#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/SequentialOutputFile.h>

#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h> // a user-defined type of particle to be inserted in ParticleSystem
#include <model/ParticleInteraction/ParticleSystem.h>

#include <model/MPI/MPIcout.h>

//#include <tutorials/ParticleInteraction/ChargedParticles/CoulombEnergy.h> // a user-defined type of energy interaction between ChargedParticle objects



using namespace model;

int main (int argc, char * argv[]) {

    // 0- define the type of ParticleSystem specifying the type of particles
    enum{dim=3}; // We work in three dimensions
    typedef ChargedParticle<dim> ChargedParticleType; // the type of particle
    typedef ParticleSystem<ChargedParticleType> ChargedParticleSystem;
    
    // 1- create a particleSystem of ChargedParticle(s)
    ChargedParticleSystem particleSystem(argc,argv); // initialized constructor: for both serial and parallel
    //ChargedParticleSystem particleSystem; // uninitialized constructor: for both serial and parallel
    //ChargedParticleSystem::initMPI(argc,argv); // if uninitialized constructor is used, call initMPI (parallel only)
    
    double cellSize=1.0;
    particleSystem.setCellSize(cellSize);
    
    // 2- add some ChargedParticle(s) to the particleSystem with random position
    //    in [-10,10] and charge=1.0
    model::cout<<"Creating particles..."<<std::endl;
    typedef typename ChargedParticle<dim>::VectorDimD VectorDimD; // helper
    for (size_t k=0;k<500000;++k)
    {
        const VectorDimD P(VectorDimD::Random()*10.0*cellSize); // PositionType::Random() returns values in [-1, 1]
        const VectorDimD V(VectorDimD::Random()*0.0);
        const double q(1.0);
        particleSystem.addParticle(P,V,q);
    }

    // Add a more particles
//    for (int k=0;k<500;++k)
//    {
////        particleSystem.addParticle(PositionType::Random()*3.0+PositionType::Ones()*2.0, 1.0); // cluster
//                particleSystem.addParticle(PositionType::Random()*10.0, 1.0); // uniform
//    }
    
    
    model::cout<<"There are "<<particleSystem.particles().size()
    <<" particles occupying "<<particleSystem.cells().size()<<" cells."<<std::endl;
    
    // -3.2a compute all binary CoulombForce interactions
    model::cout<<"Computing electric field (nearest-neighbor)..."<<std::endl;
    typedef typename ChargedParticleType::Efield Efield;
    particleSystem.computeNeighborField<Efield>();

    model::cout<<"Computing magnetic field (nearest-neighbor)..."<<std::endl;
    typedef typename ChargedParticleType::Bfield Bfield;
    particleSystem.computeNeighborField<Bfield>();

    
    // -4 output
    model::cout<<"Writing output file P/P_0.txt..."<<std::endl;
    SequentialOutputFile<'P',true> pFile1;
    pFile1<<particleSystem<<std::endl;

    
    return 0;
}


