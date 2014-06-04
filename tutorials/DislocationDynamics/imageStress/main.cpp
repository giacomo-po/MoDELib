/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#define _MODEL_NON_SINGULAR_DD_ 1 /* 1 = Cai's non-singular theory, 2 = Lazar's non-singular gradient theory */

//#define customUserOutputs "./myOutputs.h" // declare the custom output file
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/FEM/BC/LowerCorner.h>
#include <model/FEM/BC/OnMaxAxis.h>
#include <model/FEM/BC/Fix.h>


using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
    // Set up boundary conditions
    const size_t id0=DN.shared.bvpSolver.finiteElement().createNodeList<OnMaxAxis<0>>();
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id0,1);
    
    const size_t id1=DN.shared.bvpSolver.finiteElement().createNodeList<OnMaxAxis<1>>();
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id1,2);
    
    const size_t id2=DN.shared.bvpSolver.finiteElement().createNodeList<OnMaxAxis<2>>();
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id2,0);
    
    const size_t id3=DN.shared.bvpSolver.finiteElement().createNodeList<LowerCorner>();
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id3,0);
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id3,1);
    DN.shared.bvpSolver.displacement().addDirechletCondition<Fix>(id3,2);

    // Run time steps
    DN.runSteps();
    
    return 0;
}
