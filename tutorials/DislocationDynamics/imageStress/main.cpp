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

#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/DislocationDynamics/DislocationNetwork.h>

using namespace model;

int main (int argc, char* argv[])
{
    // Create the DislocationNetwork object
    typedef DislocationNetwork<3,1,CatmullRom,UniformOpen> DislocationNetworkType;
    DislocationNetworkType DN(argc,argv);
    
    // Run time steps
    DN.runSteps();
    
    // Compute the total correction force on the top face
    if(DN.shared.use_bvp)
    {
        auto topBnd=DN.shared.bvpSolver.finiteElement().boundary<AtXmax<2>,3,GaussLegendre>();
        Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Zero());
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
//        topBnd.integrate(&DN.shared.bvpSolver,temp,&BvpSolverType::bvpTraction,DN);
        topBnd.integrate(&DN.shared.bvpSolver,temp,&BvpSolverType::bvpTraction);
    }
    
    return 0;
}
