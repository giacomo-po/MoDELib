//#define customUserOutputs "./myOutputs.h" // declare the custom output file
#include <model/DislocationDynamics/DislocationNetwork.h> 

using namespace model;

int main (int argc, char* argv[]) {
    // Create a DislocationNetwork object with:
    // dim=3, tangentContinuity=1, CatmullRom splines, centripetal parametrization
    // 16 quadrature points per segment Uniformly distributed
    //DislocationNetwork<3,1,CatmullRom,centripetal,16,UniformOpen> DN(argc,argv);
     DislocationNetwork<3,1,CatmullRom,16,GaussLegendre> DN(argc,argv); // alternatively use GaussLegendre quadrature

	const double nu(Material<Isotropic>::nu);
    DN.Nsteps=1;
	DN.shared.externalStress.setZero();

	double eDot33=0.2*1.0e-9;
	const double R= 2127.0*0.5;;
	const double H=6.0*R;
	const double V=M_PI*R*R*H;
		
	for(int i=0; i<10000;i++){
        DN.runSteps();
        
		Eigen::Matrix<double,3,3> eDotP=DN.plasticStrainRate()/V;
		double dt=DN.get_dt();
		UniqueOutputFile<'S'> standard_output;
		standard_output<<i<<"  "<<DN.get_totalTime()<<" "<<dt<<"  "<<DN.network_length()<<" "<<DN.shared.externalStress(2,2)<<"  "<<eDotP.row(0)<<"  "<<eDotP.row(1)<<"  "<<eDotP.row(2)<<std::endl;
        
		DN.shared.externalStress(2,2)+=2.0*(1.0+nu)*(eDot33-eDotP(2,2))*dt;
    }
	
    return 0;
}
