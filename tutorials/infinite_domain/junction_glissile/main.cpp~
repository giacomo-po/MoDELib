/*! \page tutorial1 Formation of a glisile junction in FCC copper

*/

#define VERBOSELEVEL 4


#include <model/Dislocations/DislocationNetwork.h>
#include <model/Dislocations/DislocationSharedObjects.h>
#include <stdio.h>
#include <math.h>

#include <time.h>


using namespace model;


int main (int argc, char * const argv[]) {  
  
 DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> DN;
	
//	double t0=clock();
time_t t0;	
std::time (&t0);

	
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

	time_t t1;	
std::time (&t1);

	std::cout<<std::setprecision(8)<<std::scientific<<"simulation time="<<std::difftime (t1,t0)<<std::endl;

	
//	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;

	

	
    return 0;
}

