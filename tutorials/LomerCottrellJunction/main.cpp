/*! \page tutorial1 Formation of a glisile junction in FCC copper
 
 
 
 
 */

#define VERBOSELEVEL 4

#include <iostream>
#include <model/DislocationDynamics/DislocationNetwork.h>



using namespace model;


int main (int argc, char * const argv[]) {  
  
// DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> DN;
    DislocationNetwork<3,1,CatmullRom,centripetal,32,UniformOpen> DN(argc,argv);
//    DislocationNetwork<3,1,CatmullRom,chordal,64,UniformOpen,Copper> DN;

	
	double t0=clock();

	


	
	
	// 1: uses Nsteps from DDinput
	DN.runByStep(); 


	
	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;

	

	
    return 0;
}

