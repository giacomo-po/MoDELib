Fix fix;

const size_t nodeListID=finiteElement().template createNodeList<AtXmin<2>>();
displacement().addDirichletCondition(nodeListID,fix,{1,1,1}); // fix u1,u2,u3

// Correct for Dislocation displacement
modifyDirichletConditions(DN);

// Call solver, passing the rhs vector and the guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());