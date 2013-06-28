/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create a DislocationNetwork object
    DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DN(argc,argv);
    // Run the simulation
    
//    const char* profile_output = std::getenv("CPUPROFILE");

    
//    ProfilerStart(profile_output);
    DN.runSteps();
    
//    ProfilerStop();
    
//    for(int k=0;k<1000;++k)
//    {
//        DN.assembleAndSolve();
//    }
    
    return 0;
}

