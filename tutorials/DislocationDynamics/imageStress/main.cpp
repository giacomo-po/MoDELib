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

//#define customUserOutputs "./myOutputs.h" // declare the custom output file
//#include <model/FEM/Boundaries/LowerCorner.h>
//#include <model/FEM/Boundaries/OnMaxAxis.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
//    if(DN.shared.use_bvp)
//    {
//        Fix fix;
//        // Set up boundary conditions
//        //        auto nodeList_0=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<0>>();
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_0,1);
//        //
//        //        auto nodeList_1=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<1>>();
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_1,2);
//        //
//        //        auto nodeList_2=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<2>>();
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_2,0);
//        //
//        //        auto nodeList_3=DN.shared.bvpSolver.finiteElement().getNodeList<LowerCorner>();
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,0);
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,1);
//        //        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,2);
//        
//        auto nodeList_0=DN.shared.bvpSolver.finiteElement().getNodeList<AtXmin<2>>();
//        DN.shared.bvpSolver.addDirichletCondition(fix,nodeList_0,0); // fix x-component
//        DN.shared.bvpSolver.addDirichletCondition(fix,nodeList_0,1); // fix y-component
//        DN.shared.bvpSolver.addDirichletCondition(fix,nodeList_0,2); // fix z-component
//        
//        
//    }
    
    // Run time steps
    DN.runSteps();
    
    
    if(DN.shared.use_bvp)
    {
        auto topBnd=DN.shared.bvpSolver.finiteElement().boundary<AtXmax<2>,3,GaussLegendre>();
        Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Zero());
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        topBnd.integrate(&DN.shared.bvpSolver,temp,&BvpSolverType::traction,DN);
    }
    
    
    return 0;
}
