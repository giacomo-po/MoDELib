// BOUNDARY CONDITIONS
// Fix nodes at z=0
const size_t nodeList_0=finiteElement().template createNodeList<AtXmin<2>>();
Fix fix;
addDirichletCondition(nodeList_0,fix,{1,1,1},DN); // fix u1 u2 u3

// Displace top nodes using flat punch
const size_t nodeList_1=finiteElement().template createNodeList<FlatPunch>();
FlatPunch punch(finiteElement());
punch.updateDisplacement(DN);
addDirichletCondition(nodeList_1,punch,{0,0,1},DN); // prescribe only u3

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());