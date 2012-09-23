#define VERBOSELEVEL 4 // suppress most of the verbose output
#define customUserOutputs "myOutputs.h" // declare the custom output file
#include <model/DislocationDynamics/DislocationNetwork.h> 

using namespace model;

int main (int argc, char * const argv[]) {
    // Create a DislocationNetwork object with:
    // dim=3, tangentContinuity=1, CatmullRom splines, centripetal parametrization
    // 16 quadrature points per segment Uniformly distributed
    DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN(argc,argv);
    // DislocationNetwork<3,1,CatmullRom,centripetal,16,GaussLegendre> DN; // alternatively use GaussLegendre quadrature

    // Main simulation loop
    DN.Nsteps=1;
	DN.shared.externalStress.setZero();
	double eDot33=0.2*1.0e-9;
	const double L=4000.0;
	const double V=L*L*L;
	
	
	for(int i=0; i<1000;i++){
        DN.runByStep(true);
        
		Eigen::Matrix<double,3,3> eDotP=DN.plasticStrainRate()/V;
		double dt=DN.get_dt();
		UniqueOutputFile<'S'> standard_output;
		standard_output<<i<<"  "<<DN.get_totalTime()<<" "<<dt<<"  "<<DN.network_length()<<" "<<DN.shared.externalStress(2,2)<<"  "<<eDotP.row(0)<<"  "<<eDotP.row(1)<<"  "<<eDotP.row(2)<<std::endl;
        
		DN.shared.externalStress(2,2)+=2.0*(1+Material<Isotropic>::nu)*(eDot33-eDotP(2,2))*dt;
    }


	
//	std::cout<<"simulation time= "<<(clock()-t0)/CLOCKS_PER_SEC<<" [sec]"<<std::endl;
    return 0;
}




   // double nu=Material<Isotropic>::nu;

//	if (argc==1){
//		// argv[0] is by default the name of the executable so use default name for inputfile
//		DN.read("./","DDinput.txt");
//	}
//	else{
//		// argv[1] is assumed to be the filename with working directory the current directory
//		DN.read("./",argv[1]);
//	}
