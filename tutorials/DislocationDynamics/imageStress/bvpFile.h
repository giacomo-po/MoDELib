Fix fix;

const size_t nodeListID=finiteElement().template createNodeList<AtXmin<2>>();
addDirichletCondition(nodeListID,fix,{1,1,1},DN); // fix x-component


// Call solver, passing the rhs vector and the guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());