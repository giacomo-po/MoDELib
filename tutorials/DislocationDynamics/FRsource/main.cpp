#include <model/Quadrature/UniformOpen.h> 
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create a DislocationNetwork object
    DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN(argc,argv);
    // Run the simulation
    DN.runSteps();
    return 0;
}

