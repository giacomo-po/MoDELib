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
#include <Eigen/Dense>
#include <model/FEM/FiniteElement.h>
#include <model/FEM/Elements/LagrangeElement.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>

using namespace model;

struct MyFunction : public EvalFunction<MyFunction>
{
    
    constexpr static int rows=1;
    constexpr static int cols=1;
    
    
//    const double dt;
//    
//    MyFunction(const double& dt_in) : dt(dt_in){}
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    double operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the constant c.
      */
        
        const Eigen::Matrix<double,ElementType::dim,1> P=ele.position(bary);
        
//        return dt*(abs(P(0))<0.1);
        return P(0);

    }
    
};


int main(int argc, char** argv)
{
    
    
    // Take meshID as a user input
    int meshID(0);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    const int dim=1;
    const int SForder=1;
    
    // Create mesh and read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    SimplicialMesh<dim> mesh;
    mesh.readMesh(meshID);
    
    // Create a FiniteElement object on the mesh
    typedef LagrangeElement<dim,SForder> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    
    FiniteElementType fe(mesh);
    
    
    auto f=fe.trial<'f',1>();
//    auto gradf=grad(f);
    
    auto dV=fe.domain<EntireDomain,3,GaussLegendre>();
    auto lhs=(test(grad(f)),grad(f))*dV;
    MyFunction g;
    auto rhs=(test(f),g)*dV;

    //    auto g=make_constant(1.0);
    //MyFunction g(1.0);
    
    


    auto weakProblem(lhs==rhs); //  weak problem

    Fix fix;
    
    const size_t nodeList_0=fe.createNodeList<AtXmin<0>>();
    f.addDirichletCondition(nodeList_0,fix,{1}); // set f(x_min)=0

    const size_t nodeList_1=fe.createNodeList<AtXmax<0>>();
    f.addDirichletCondition(nodeList_1,fix,{1}); // set f(x_min)=0

    
    weakProblem.assembleWithMasterSlaveConstraints();
    
    const double solverTolerance=1.0e-6;
    f=weakProblem.solveWithGuess(solverTolerance,f);

    SequentialOutputFile<'D',1> uFile;
    uFile<<f.onDomain();
    
    return 0;
}


