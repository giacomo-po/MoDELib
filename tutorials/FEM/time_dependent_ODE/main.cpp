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

using namespace model;



/**************************************************************************/
/**************************************************************************/
struct MyFunction //: public EvalExpression<Constant<T,_rows,_cols> >
{
    
    constexpr static int rows=1;
    constexpr static int cols=1;
    
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    double operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the constant c.
      */
        
        const Eigen::Matrix<double,ElementType::dim,1> P=ele.position(bary);
        
        return 1.0;
    }
    
    /**********************************************************************/
    EvalExpression<MyFunction> eval() const
    {
        return EvalExpression<MyFunction>(*this);
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
    
    // Create a 3d-SimplicialMesh object
    SimplicialMesh<3> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    
    // Create a FiniteElement object on the mesh
    typedef LagrangeElement<3,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    FiniteElementType fe(mesh);
    
    
    auto rho=fe.trial<1>();
    rho=1.0;
    
    
    auto rhoDot=fe.trial<1>();
    
    //auto a=2.0*rhoDot;
    
    MyFunction f;
    
    /**************************************************************************/
    // Create lhs (BilinearWeakForm) and rhs (LinearWeakForm)
    auto dV=fe.domain<EntireDomain,5,GaussLegendre>();
    auto bwf=(rhoDot.test(),rhoDot)*dV;
    auto lwf=(rhoDot.test(),f.eval())*dV;
    
    // Create the WeakProblem
    auto weakProblem(bwf=lwf); //  weak problem
    
    // Assemble WeakProblem
    weakProblem.assembleWithMasterSlaveConstraints();
    
    /**************************************************************************/
    // Solve WeakProblem
    const double solverTolerance=1.0e-6;
    rhoDot=weakProblem.solve(solverTolerance);
    
//    rho += rhoDot*dt;
    
    /**************************************************************************/
    // Output displacement and stress on external mesh faces
    SequentialOutputFile<'D',1> uFile;
    uFile<<rhoDot;
    
    
    auto rho1 = 2.0*(2.0*rho);
    
    Eigen::Matrix<double,3,1> P;
    P<<0.5,0.5,0.5;
    
    std::cout<<rho1(P)<<std::endl;
    
    return 0;
}


