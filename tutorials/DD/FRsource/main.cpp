#define VERBOSELEVEL 4 // suppress most of the verbose output
#define customUserOutputs "myOutputs.h" // declare the custom output file
#include <model/quadrature/UniformOpen.h>

#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char * const argv[]) {
    // Create a DislocationNetwork object with:
    // dim=3, tangentContinuity=1, CatmullRom splines, centripetal parametrization
    // 16 quadrature points per segment Uniformly distributed
    //DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN(argc,argv);
//     DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre> DN(argc,argv); // alternatively use GaussLegendre quadrature
    DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN(argc,argv); // alternatively use GaussLegendre quadrature

    // Main simulation loop
//    DN.Nsteps=1;
//	DN.shared.externalStress.setZero();
	
        DN.runByStep(true);

    return 0;
}

