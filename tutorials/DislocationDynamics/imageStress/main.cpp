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
#include <model/FEM/Boundaries/LowerCorner.h>
#include <model/FEM/Boundaries/OnMaxAxis.h>
#include <model/FEM/Boundaries/AtXmax.h>

#include <model/FEM/BoundaryConditions/Fix.h>


using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
    if(DN.shared.use_bvp)
    {
        Fix fix;
        // Set up boundary conditions
        auto nodeList_0=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<0>>();
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_0,1);
        
        auto nodeList_1=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<1>>();
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_1,2);
        
        auto nodeList_2=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<2>>();
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_2,0);
        
        auto nodeList_3=DN.shared.bvpSolver.finiteElement().getNodeList<LowerCorner>();
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,0);
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,1);
        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,2);
    }
    // Run time steps
    DN.runSteps();
    
    
    auto topBnd=DN.shared.bvpSolver.finiteElement().boundary<AtXmax<2>,3,GaussLegendre>();
    
    Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Zero());
    
    typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
    
    topBnd.integrate(&DN.shared.bvpSolver,temp,&BvpSolverType::traction,DN);
    
    
    return 0;
}
