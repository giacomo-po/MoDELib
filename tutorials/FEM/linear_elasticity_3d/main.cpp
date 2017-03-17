/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more
//#define EIGEN_STACK_ALLOCATION_LIMIT 1000000


#include <iostream>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/FEM/FiniteElement.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/FEM/Domains/SubDomain.h>

using namespace model;



int main(int argc, char** argv)
{
    
    // Take meshID as a user input
    int meshID(0);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    // Create a 3d-SimplicialMesh object
    SimplicialMesh<3> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    
    // Create a FiniteElement object on the mesh
    typedef LagrangeElement<3,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    FiniteElementType fe(mesh);
    
    /**************************************************************************/
    // Define some constants for the 3d elastic problem
    const double mu =75.6;  // GPa (for Cu)
    const double lam=119.9; // GPa (for Cu)
    const double C11(lam+2.0*mu);
    const double C12(lam);
    const double C44(mu);
    
    Eigen::Matrix<double,6,6> C;
    C  << C11, C12, C12, 0.0, 0.0, 0.0,
    /***/ C12, C11, C12, 0.0, 0.0, 0.0,
    /***/ C12, C12, C11, 0.0, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
    C/=mu; // make dimensionless
    
    /**************************************************************************/
    // Define trial function (displacement field) u and related expressions
    auto u=fe.trial<'u',3>();       // displacement field u=[u1; u2; u3]
    auto b=grad(u);             // displacement gradient b=[u1,1; u1,2; u2,1; u2,2]
    auto e=def(u);              // engineering strain e=[u1,1; u2,2; u3,3; u1,2+u2,1; u2,3+u3,2; u1,3+u3,1]
    auto s=C*e;                 // stress field s=[s11; s22; s33; s12; s23; s13]
    
    /**************************************************************************/
    // Create the BilinearWeakForm bWF_u=int(test(e)^T*s)dV
    auto dV=fe.domain<EntireDomain,4,GaussLegendre>();
    auto bWF_u=(test(e),s)*dV;
    
    /**************************************************************************/
    // Create the LinearWeakForm lWF_1=int(test(u)^T*f)ndA
    // f is a constant traction vector.
//    Eigen::Matrix<double,3,1> f((Eigen::Matrix<double,3,1>()<<0.0,0.000,0.00).finished());
    auto f=make_constant((Eigen::Matrix<double,3,1>()<<0.0,0.000,0.00).finished());
    auto dA_1=fe.boundary<AtXmax<2>,3,GaussLegendre>();
    auto lWF_1=(test(u),f)*dA_1;
    
    /**************************************************************************/
    // Create the LinearWeakForm lWF_2=int(test(u)^T*p)ndA
    // p is a constant hydrostatic tensor. The traction vector will be t=p*n;
    Eigen::Matrix<double,3,3> p(-0.01*Eigen::Matrix<double,3,3>::Identity());
//    auto p=make_constant(-0.01*Eigen::Matrix<double,3,3>::Identity());
    auto ndA_2=fe.boundary<ExternalBoundary,3,GaussLegendre>();
    auto lWF_2=(test(u),make_constant(p))*ndA_2;
    
    /**************************************************************************/
    // Create the LinearWeakForm lWF_1=int(test(u)^T*f)ndA
    // a is a constant boby force vector.
//    Eigen::Matrix<double,3,1> a((Eigen::Matrix<double,3,1>()<<0.0,0.000,-0.00001).finished());
    auto a=make_constant((Eigen::Matrix<double,3,1>()<<0.0,0.000,-0.00001).finished());
    auto lWF_3=(test(u),a)*dV;
    
    /**************************************************************************/
    // Create the WeakProblem
    auto weakProblem(bWF_u==lWF_1+lWF_2+lWF_3); //  weak problem
    
    /**************************************************************************/
    // Set up Dirichlet boundary conditions
    // Create a list of nodes having x(0)=x0_min, where x0_min is the minimum value among the fe nodes
    const size_t nodeList_0=fe.createNodeList<AtXmin<2>>();
    // Fix those those nodes
    Fix fix;
    u.addDirichletCondition(nodeList_0,fix,{1, 1, 1}); // fix u1, u2, u3
    
    // Create a list of nodes having x(2)=x2_max, where x2_max is the maximum value among the fe nodes
    //const size_t nodeList_1=fe.createNodeList<AtXmax<2>>();
    //u.addDirichletCondition(nodeList_1,fix,{0, 0, 1}); // fix only u3
    
    
    /**************************************************************************/
    // Assemble
    //    weakProblem.assembleWithLagrangeConstraints();
    //    weakProblem.assembleWithPenaltyConstraints(1000000.0);
    weakProblem.assembleWithMasterSlaveConstraints();
    
    /**************************************************************************/
    // Solve for u using current value as guess solution
    const double solverTolerance=1.0e-6;
    u=weakProblem.solveWithGuess(solverTolerance,u);
    
    /**************************************************************************/
    // Output displacement and stress on external mesh faces
    SequentialOutputFile<'U',1> uFile;
    uFile<<u.onBoundary();
    
    SequentialOutputFile<'S',1> sFile;
    sFile<<s.onBoundary();
    
    return 0;
}




