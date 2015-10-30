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
#define userBVPfile "./bvpFile.h"
#define userOutputFile "./myOutputs.h" // declare the custom output file

#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <./TensionTorsioner.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
    // Read TensionTorsioner data
    model::EigenDataReader EDR;

    EDR.readScalarInFile("./loadInput.txt","epsilonDot",TensionTorsioner::epsilonDot);
    EDR.readScalarInFile("./loadInput.txt","gammaDot",TensionTorsioner::gammaDot);
    EDR.readScalarInFile("./loadInput.txt","initialDisplacement",TensionTorsioner::initialDisplacement);
    EDR.readScalarInFile("./loadInput.txt","initialTwist_Rad",TensionTorsioner::initialTwist_Rad);
    EDR.readScalarInFile("./loadInput.txt","relaxSteps",TensionTorsioner::relaxSteps);
    EDR.readScalarInFile("./loadInput.txt","apply_tension",TensionTorsioner::apply_tension);
    EDR.readScalarInFile("./loadInput.txt","apply_torsion",TensionTorsioner::apply_torsion);

    // Initialize TensionTorsioner
    TensionTorsioner::init(DN.runningID(),DN.userOutputColumn());

    // Run time steps
    DN.runSteps();
    
    return 0;
}
