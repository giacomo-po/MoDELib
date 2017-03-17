// BOUNDARY CONDITIONS
// Fix nodes at z=0
const size_t nodeList_0=finiteElement().template createNodeList<AtXmin<2>>();
Fix fix;
displacement().addDirichletCondition(nodeList_0,fix,{1,1,1}); // fix u1 u2 u3

// Displace top nodes using flat punch
if(DN.runningID()>=FlatPunch::relaxSteps)
{
const size_t nodeList_1=finiteElement().template createNodeList<FlatPunch>();
FlatPunch punch(finiteElement());
punch.updateDisplacement(DN);
displacement().addDirichletCondition(nodeList_1,punch,{0,0,1}); // prescribe only u3
}

// Correct for Dislocation displacement
modifyDirichletConditions(DN);

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());