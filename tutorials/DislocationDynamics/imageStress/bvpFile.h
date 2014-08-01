Fix fix;
// Set up boundary conditions
//        auto nodeList_0=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<0>>();
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_0,1);
//
//        auto nodeList_1=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<1>>();
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_1,2);
//
//        auto nodeList_2=DN.shared.bvpSolver.finiteElement().getNodeList<OnMaxAxis<2>>();
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_2,0);
//
//        auto nodeList_3=DN.shared.bvpSolver.finiteElement().getNodeList<LowerCorner>();
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,0);
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,1);
//        DN.shared.bvpSolver.displacement().addDirichletCondition(fix,nodeList_3,2);

//const bool useDisplacementMultipole=true;

auto nodeList_0=finiteElement().template getNodeList<AtXmin<2>>();
addDirichletCondition(fix,nodeList_0,0,DN); // fix x-component
addDirichletCondition(fix,nodeList_0,1,DN); // fix y-component
addDirichletCondition(fix,nodeList_0,2,DN); // fix z-component


// Call solver, passing the rhs vector and the guess
displacement()=solve(-lwf.globalVector(),displacement());