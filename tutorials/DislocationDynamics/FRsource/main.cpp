#include <model/Quadrature/UniformOpen.h> 
#include <model/DislocationDynamics/DislocationNetwork.h>

//#include <cstdlib>
//#include <gperftools/profiler.h>

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

