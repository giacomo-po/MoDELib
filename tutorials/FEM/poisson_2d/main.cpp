#include <iostream>
#include <stdlib.h>     /* atoi */

#include <model/Utilities/SequentialOutputFile.h>

#include<model/FEM/FiniteElement.h>




using namespace model;

struct FixXmin //: public DirichletCondition<TrialFunctionType>
{
    
    template <typename NodeType>
    std::pair<bool,double> operator()(const NodeType& node) const
    {
        return std::pair<bool,double>(node.p0(0)==0,0.0);
    }
    
};

struct FixXmax //: public DirichletCondition<TrialFunctionType>
{
    
    template <typename NodeType>
    std::pair<bool,double> operator()(const NodeType& node) const
    {
        return std::pair<bool,double>(node.p0(0)==300.0,node.p0(1));
    }
    
};


struct FixBottom //: public DirichletCondition<TrialFunctionType>
{
    
    template <typename NodeType>
    std::pair<bool,double> operator()(const NodeType& node) const
    {
        return std::pair<bool,double>(node.p0(1)==-550.0,0.0);
    }
    
};

struct FixLeft //: public DirichletCondition<TrialFunctionType>
{
    
    template <typename NodeType>
    std::pair<bool,double> operator()(const NodeType& node) const
    {
        return std::pair<bool,double>(node.p0(0)==0.0,0.0);
//        return std::pair<bool,double>(node.p0(0)==0.0 && node.p0(1)>-300.0,0.0);

    }
    
};

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
    
    int meshID(2);
    
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    //TO DO:
    //    -TrialFunction and TestFunction "could" use different shapefunctions
    
    auto kappa=make_constant(3.0);
//    Eigen::Matrix<double,2+1,1> a;
//    a<<0.5,0.5,0;
    
    //-1 Create mesh and read from file
    SimplicialMesh<2> mesh;
    mesh.readMesh(meshID);
    
    //-2 Create a FiniteElement on the mesh
    typedef LagrangeElement<2,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    
    FiniteElementType fe(mesh);

//    // Test outNormal
//    const double L0(0.0);
//    const int np=10;
//    const double dL(1.0/np);
//
//    SequentialOutputFile<'A',1> aFile;
//
//    Eigen::Matrix<double,2,1> dx(-0.1,0.1);
//    fe.element(0).Xe.col(2)+=dx;
//    
//    for (int k=0;k<np+1;++k)
//    {
//        const double L1(k*dL);
//        const double L2(1.0-L1);
//        const Eigen::Matrix<double,3,1> bary(L1,L0,L2);
//
//        aFile<<fe.elementBegin()->position(bary).transpose()<<" "
//             <<fe.elementBegin()->outNormal(bary,1).transpose()<<"\n";
//    }
  
    /**********************/
    //-3 Create a TrialFunction from the FiniteElement
    auto T=fe.trial<1>();
    
    auto q=grad(T);
    
    //-4 Create a weak form
    auto wf_T=(q.test(),q);
    
    //FixXmin fixMin;
    
    
    // assemble the weak form
    wf_T.assembleOnDomain<3,GaussLegendre>();
    //wf_T.assembleMatrix();
    //wf_T.assembleConstrainedMatrix();
    //    wf_T.prune();
    
//    auto wp_T(wf_T=1);
//    
//    T.addDirechletCondition(FixXmin(),0);
//    T.addDirechletCondition(FixXmax(),0);
//    wp_T.assembleWithLagrangeConstraints();
//    wp_T.solve();
//    wp_T.output();
//    
//    T.test(),1.0;
    
    
    /**********************/
    const double mu =75.6;  // GPa (for Cu)
    const double lam=119.9; // GPa (for Cu)
    Eigen::Matrix<double,4,4> C;
    C<< lam+2.0*mu, 0.0, 0.0,        lam,
    /**/       0.0,  mu,  mu,        0.0,
    /**/       0.0,  mu,  mu,        0.0,
    /**/       lam, 0.0, 0.0, lam+2.0*mu;
    C/=mu; // make dimensionless
    
    Eigen::Matrix<double,2,1> f;
    f<<0.0,-0.001;

    
    auto u=fe.trial<2>();       // displacement field u=[u1; u2]
    auto b=grad(u);             // displacement gradient b=[u1,1; u1,2; u2,1; u2,2]
    auto s=C*b;                 // stress field s=[s11; s12; s21; s22]
    auto bWF_u=(b.test(),s);    // bilinear weak form int(test(b)^T*s)dV
    bWF_u.assembleOnDomain<3,GaussLegendre>();

    
    
    auto lWF_u=(u.test(),f);
    lWF_u.assembleOnDomain<3,GaussLegendre>(); // linear weak form int(test(u)^T*f)dV
//    lWF_u.assembleOnExternalBoundary<4,GaussLegendre>();

    
    auto wp_u(bWF_u=lWF_u); //  weak problem

    u.addDirechletCondition(FixBottom(),1);
    u.addDirechletCondition(FixBottom(),0);

//    u.addDirechletCondition(FixLeft(),0);
//    u.addDirechletCondition(PushTop(),1);

    
    wp_u.assembleWithLagrangeConstraints();
    
    wp_u.solve();
    wp_u.output();
//
//    
    
//    b(x);
//    s(x);
    
    
    //    auto kappa=make_constant(3.0);
    //    Eigen::Matrix<double,2+1,1> b;
    //    b<<0.5,0.5,0;
    //
    //    std::cout<<T.sfm(*fe.elementBegin(),b)<<std::endl;
    //
    //    Eigen::Matrix<double,2,2> m(Eigen::Matrix<double,2,2>::Ones());
    //
    //    Eigen::Matrix<double,2,1> a(Eigen::Matrix<double,2,1>::Ones());
    //
    //
    //    auto M=make_constant(m);
    //
    //    auto r=3*(2*grad(T)+m*grad(2.0*T+T)+a*T);
    //
    //
    ////    std::cout<<"Computing sfm of first element:"<<std::endl;
    ////    std::cout<<fe.elementBegin()->sfm(b)<<std::endl;
    ////
    ////
    ////
    ////    std::cout<<"Computing diff of of first element:"<<std::endl;
    ////    std::cout<<fe.elementBegin()->diff(b)<<std::endl;
    //

    
//    SequentialOutputFile<'P',1> pFile;
//    for (typename FiniteElementType::NodeContainerType::const_iterator nIter=fe.nodeBegin();nIter!=fe.nodeEnd();++nIter)
//    {
//        pFile<<nIter->p0.transpose()<<std::endl;
//    }
    
    SimplicialMesh<3> mesh3;
    mesh.readMesh(7);
    
    
	return 0;
}




