/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// Define the non-singluar method used for calculations
#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

// Select the external load controller (if nothing is defined DummyExternalLoadController.h is used)
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/UniformExternalLoadController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/SampleUserStressController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/SequentialTorsionTensionController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/ClockIndentationController.h>

//#include <model/DislocationDynamics/DislocationNetwork.h>
#include <DefectiveCrystal.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
//    DislocationNetwork<3,0,Hermite> DN(argc,argv);
    DefectiveCrystal<3,0,Hermite> DC(argc,argv);
    
    // Run time steps
    DC.runSteps();

    DC.DN->bvpSolver->vacancyConcentration()=1.0;
    DC.DN->vacancyConcentrationAssembly();
    
    std::ofstream file("concentration.txt");
    int N=500;
    double dx(500.0/N);
    for(int k=-N;k<N+1;++k)
    {
        Eigen::Matrix<double,3,1> x(k*dx,0,0);
        file<<x.transpose()<<" "<<DC.DN->vacancyConcentration(x)<<"\n";
    }
    
    return 0;
}
