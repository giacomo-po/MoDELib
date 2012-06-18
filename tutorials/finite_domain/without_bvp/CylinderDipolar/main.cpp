


#define VERBOSELEVEL 4

#define customUserOutputs "myOutputs.h"

#include <model/Dislocations/DislocationNetwork.h>

using namespace model;

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


	double nu=DN.shared.material.nu;
    DN.Nsteps=1;
	DN.shared.externalStress.setZero();
	double eDot33=0.2*1.0e-10;
	const double R=2000.0;
	const double H=12000.0;
	const double V=M_PI*R*R*H;
	
	
	for(int i=0; i<10000;i++){
        DN.runByStep(true);

		Eigen::Matrix<double,3,3> eDotP=DN.plasticStrainRate()/V;
		double dt=DN.get_dt();
		UniqueOutputFile<'S'> standard_output;
		standard_output<<i<<"  "<<DN.get_totalTime()<<" "<<dt<<"  "<<DN.network_length()<<" "<<DN.shared.externalStress(3.3)<<"  "<<eDotP.row(0)<<"  "<<eDotP.row(1)<<"  "<<eDotP.row(2)<<std::endl;

		DN.shared.externalStress(3.3)+=2.0*(1+nu)*(eDot33-eDotP(2,2))*dt;
    }
    
    
    
	
	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;

	
    return 0;
}
