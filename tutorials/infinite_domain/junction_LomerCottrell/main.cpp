// TO DO:

// - Using Variadics implement integrator


// - move all the geometric stuff that is dislocation-only related to the Dislocation layer (especially quadrature)

#define VERBOSELEVEL 4

//#define UserStressFile	"mmdl/Dislocations/ExternalStresses/shearYZ.h"
//#define UpdateBoundaryConditionsFile "examples/Indentation_5m/UpdateBoundaryConditions.h"

//#define Contact_Loading       

//#define customUserOutputs "examples/Indentation_5m/myOutputs.h"

#include <mmdl/Dislocations/DislocationNetwork.h>
#include <mmdl/Dislocations/DislocationSharedObjects.h>
#include <stdio.h>
#include <math.h>

//#include "mmdl/BVP/UpdateBoundaryConditions.h"


using namespace mmdl;


int main (int argc, char * const argv[]) {  
  
 DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> DN;
	
	double t0=clock();
	
	if (argc==1){
		// argv[0] is by default the name of the executable so use default name for inputfile 
		DN.read("./","DDinput.txt");
	}
	else{
		// argv[1] is assumed to be the filename with working directory the current directory 
		DN.read("./",argv[1]);
	}

	
	
	// 1: uses Nsteps from DDinput
	DN.runByStep(); 

	
	
	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;

	

	
    return 0;
}

