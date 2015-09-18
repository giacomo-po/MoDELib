Fix fix;

const size_t nodeListID=finiteElement().template createNodeList<AtXmin<1>>();
addDirichletCondition(nodeListID,fix,{1,1,1},DN); // fix u1,u2,u3


// Call solver, passing the rhs vector and the guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());