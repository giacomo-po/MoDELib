#include <iostream>

#include <model/Utilities/SequentialOutputFile.h>

#include <model/FEM/FiniteElement.h>
#include <model/FEM/Domains/TopBoundary.h>
#include <model/FEM/BC/Fix.h>
#include <model/FEM/BC/AtXmin.h>

using namespace model;

struct PushTop //: public DirichletCondition<TrialFunctionType>
{
    
    template <typename NodeType>
    std::pair<bool,double> operator()(const NodeType& node) const
    {
        return std::pair<bool,double>(node.p0(1)==0.0 && node.p0(0)<=50.0,-5);
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
    
    
    /**********************/
    const double mu =75.6;  // GPa (for Cu)
    const double lam=119.9; // GPa (for Cu)
    const double C11(lam+2.0*mu);
    const double C12(lam);
    const double C44(2.0*mu); // this multiplies a true strain, so 2 is necessary
    
    Eigen::Matrix<double,6,6> C;
    C  << C11, C12, C12, 0.0, 0.0, 0.0,
    /***/ C12, C11, C12, 0.0, 0.0, 0.0,
    /***/ C12, C12, C11, 0.0, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
    C/=mu; // make dimensionless
    
    
    auto u=fe.trial<3>();       // displacement field u=[u1; u2; u3]
    auto b=grad(u);             // displacement gradient b=[u1,1; u1,2; u2,1; u2,2]
    auto e=def(u);              // strain e=[u1,1; u2,2; u3,3; u1,2+u2,1; u2,3+u3,12; u1,3+u3,1]
    auto s=C*e;                 // stress field s=[s11; s12; s21; s22]

    
    // Create the BilinearWeakForm bWF_u=int(test(e)^T*s)dV
    ElementaryDomain<3,4,GaussLegendre>  dV;
    auto bWF_u=(e.test(),s)*dV;
    
    Eigen::Matrix<double,3,1> f;
    f<<0.0,0.0,0.1;
    Eigen::Matrix<double,3,3> p(0.1*Eigen::Matrix<double,3,3>::Identity());
    
//    auto ndA=fe.externalBoundary<TopBoundary,3,GaussLegendre>();
    //auto  dV=fe.domain<Top,3,GaussLegendre>();
//    auto ndA=fe.boundary<ExternalBoundary,3,GaussLegendre>();
    auto ndA=fe.boundary<TopBoundary,3,GaussLegendre>();
  
    
    // Create the LinearWeakForm lWF_u=int(test(u)^T*f)ndS
    auto lWF_u=(u.test(),f)*ndA;
    
    // Create the WeakProblem
    auto wp_u(bWF_u=lWF_u); //  weak problem
    
    
    /*******************************************************/
    // Set up Dirichlet boundary conditions
    
    // Create a list of nodes having x(0)=x0_min, where x0_min is the minimum value among the fe nodes
    const size_t nodeListID_0=fe.createNodeList<AtXmin<0> >();
    // Fix the 0-th component of displacement for those nodes
    u.addDirechletCondition<Fix>(nodeListID_0,0);

    // Create a list of nodes having x(1)=x1_min, where x1_min is the minimum value among the fe nodes
    const size_t nodeListID_1=fe.createNodeList<AtXmin<1> >();
    // Fix the 1-st component of displacement for those nodes
    u.addDirechletCondition<Fix>(nodeListID_1,1);
    
    // Create a list of nodes having x(2)=x2_min, where x2_min is the minimum value among the fe nodes
    const size_t nodeListID_2=fe.createNodeList<AtXmin<2> >();
    // Fix the 2-nd component of displacement for those nodes
    u.addDirechletCondition<Fix>(nodeListID_2,2);

    

    /*******************************************************/
    // Solve
    
    wp_u.assembleWithLagrangeConstraints();
//    wp_u.assembleWithPenaltyConstraints(1000.0);
    
    wp_u.solve(0.0001);
    //wp_u.output();
    u.dofContainer=wp_u.x.segment(0,u.dofContainer.size());


    SequentialOutputFile<'U',1> uFile;
    uFile<<u;
    
    SequentialOutputFile<'S',1> sFile;
    sFile<<s;
    
    
	return 0;
}




