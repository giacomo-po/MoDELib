#include <iostream>
#include <float.h>     /* DBL_MAX */
#include <stdlib.h>    /* atoi */

#include <model/Utilities/SequentialOutputFile.h>

#include <model/FEM/FiniteElement.h>
#include <model/FEM/Domains/TopBoundary.h>
#include <model/FEM/Domains/FixMin.h>

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
    
    int meshID(7);
    
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    
    SimplicialMesh<3> mesh;
    mesh.readMesh(meshID);
    
    
    typedef LagrangeElement<3,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    
    FiniteElementType fe(mesh);
    
    
    /**********************/
    const double mu =75.6;  // GPa (for Cu)
    const double lam=119.9; // GPa (for Cu)
    const double C11(lam+2.0*mu);
    const double C12(lam);
    const double C44(mu); // this multiplies a true strain, so 2 is necessary
    
    Eigen::Matrix<double,6,6> C;
    C  << C11, C12, C12, 0.0,     0.0,     0.0,
    /***/ C12, C11, C12, 0.0,     0.0,     0.0,
    /***/ C12, C12, C11, 0.0,     0.0,     0.0,
    /***/ 0.0, 0.0, 0.0, 2.0*C44, 0.0,     0.0,
    /***/ 0.0, 0.0, 0.0, 0.0,     2.0*C44, 0.0,
    /***/ 0.0, 0.0, 0.0, 0.0,     0.0,     2.0*C44;
    C/=mu; // make dimensionless
    
    
    auto u=fe.trial<3>();       // displacement field u=[u1; u2; u3]
    auto b=grad(u);             // displacement gradient b=[u1,1; u1,2; u2,1; u2,2]
    auto e=def(u);              // strain e=[u1,1; u2,2; u3,3; u1,2+u2,1; u2,3+u3,12; u1,3+u3,1]
    auto s=C*e;                 // stress field s=[s11; s12; s21; s22]
    
    ElementaryDomain<3,4,GaussLegendre>  dV;
    
    
    auto bWF_u=(e.test(),s)*dV;    // bilinear weak form int(test(b)^T*s)dV
    
    Eigen::Matrix<double,3,1> f;
    f<<0.0,0.0,0.1;
    Eigen::Matrix<double,3,3> p(0.1*Eigen::Matrix<double,3,3>::Identity());
    
//    auto ndA=fe.externalBoundary<TopBoundary,3,GaussLegendre>();
    //auto  dV=fe.domain<Top,3,GaussLegendre>();
//    auto ndA=fe.boundary<ExternalBoundary,3,GaussLegendre>();
    auto ndA=fe.boundary<TopBoundary,3,GaussLegendre>();
  
    
    auto lWF_u=(u.test(),f)*ndA;

    auto wp_u(bWF_u=lWF_u); //  weak problem
    
    
//    DirichletCondition<FixBottom> d1(fe);
    
    FixMin<0> fm0(fe);
    FixMin<1> fm1(fe);
    FixMin<2> fm2(fe);
    
    u.addDirechletCondition(fm0,0);
    u.addDirechletCondition(fm1,1);
    u.addDirechletCondition(fm2,2);

    
//    wp_u.assembleWithLagrangeConstraints();
    wp_u.assembleWithPenaltyConstraints(1000.0);
    
    wp_u.solve(0.0001);
    //wp_u.output();
    u.dofContainer=wp_u.x.segment(0,u.dofContainer.size());


    SequentialOutputFile<'U',1> uFile;
    uFile<<u;
    
    SequentialOutputFile<'S',1> sFile;
    sFile<<s;
    
    
	return 0;
}




