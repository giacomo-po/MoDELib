/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

/* Define the non-singluar method used for calculations
 * _MODEL_NON_SINGULAR_DD_ = 0 classical theory
 * _MODEL_NON_SINGULAR_DD_ = 1 Cai's regularization method
 * _MODEL_NON_SINGULAR_DD_ = 2 Lazar's regularization method
 */
#define _MODEL_NON_SINGULAR_DD_ 0
#define userLoadController "./LoadController.h"
//#define DislocationNucleationFile
//#define userBVPfile "./bvpFile.h"
//#define userOutputFile "./myOutputs.h" // declare the custom output file
//#include <model/FEM/Boundaries/AtXmin.h>
//#include <model/FEM/Boundaries/AtXmax.h>
//#include <./FlatPunch.h>
//#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
//    // Read punch data
//    model::EigenDataReader EDR;
//    EDR.readScalarInFile("./punchInput.txt","punchSize",FlatPunch::punchSize);
//
//    EDR.readScalarInFile("./punchInput.txt","strainRate",FlatPunch::strainRate);
//    EDR.readScalarInFile("./punchInput.txt","initialDisplacement",FlatPunch::initialDisplacement);
//    EDR.readScalarInFile("./punchInput.txt","relaxSteps",FlatPunch::relaxSteps);
    
    // Run time steps
    //DN.runSteps();
    DN.runDDD();
    
    return 0;
}
