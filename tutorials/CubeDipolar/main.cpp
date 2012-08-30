/*! \page tutorial1 Formation of a glisile junction in FCC copper
 
 
 
 
 */

#define VERBOSELEVEL 4

#include <iostream>
#include <model/Dislocations/DislocationNetwork.h>



using namespace model;


int main (int argc, char * const argv[]) {  
  
// DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> DN;
    DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN;
 //   DislocationNetwork<3,1,CatmullRom,chordal,32,UniformOpen,Copper> DN;

	
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

