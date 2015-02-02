// BOUNDARY CONDITIONS
// Fix nodes at z=0
auto nodeList_0=finiteElement().template getNodeList<AtXmin<2>>();
Fix fix;
addDirichletCondition(fix,nodeList_0,0,DN); // fix x-component of displacement
addDirichletCondition(fix,nodeList_0,1,DN); // fix y-component of displacement
addDirichletCondition(fix,nodeList_0,2,DN); // fix z-component of displacement

// Displace top nodes using flat punch
auto nodeList_1=finiteElement().template getNodeList<FlatPunch>();
FlatPunch punch(finiteElement());
punch.updateDisplacement(DN);
addDirichletCondition(punch,nodeList_1,2,DN); // prescribe z-component of displacement

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());