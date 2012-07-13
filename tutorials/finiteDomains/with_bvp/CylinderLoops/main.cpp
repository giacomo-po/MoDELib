// TO DO:

// - Using Variadics implement integrator


// - move all the geometric stuff that is dislocation-only related to the Dislocation layer (especially quadrature)

#define VERBOSELEVEL 4

//#define UserStressFile	"mmdl/Dislocations/ExternalStresses/shearYZ.h"
//#define UpdateBoundaryConditionsFile "/Users/giacomo/Documents/UCLA/Research/C++/Mmdl/examples/TamerLoop/UpdateBoundaryConditions.h"
#define customUserOutputs "myOutputs.h"
//#define customUserOutputs "/Users/giacomo/Documents/UCLA/Research/C++/Mmdl/examples/GiacomoDeformed/myOutputs.h"
#include <mmdl/Dislocations/DislocationNetwork.h>





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


//	for(int i=0; i<10;i++){
        // change BC
        // compute Nstps or timeWindow
    //    DN.Nsteps=200;

        DN.runByStep(true);
//    }
    
    
//    // ....
//    
//    DN.runByStep();
//    
//    
//    DN.timeWidow=1234.234527;
//    
    
	
	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;

	
    return 0;
}
