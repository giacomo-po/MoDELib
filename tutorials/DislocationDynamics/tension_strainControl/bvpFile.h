// BOUNDARY CONDITIONS
// Fix nodes at z=0
const size_t nodeList_0=finiteElement().template createNodeList<AtXmin<2>>();
Fix fix;
displacement().addDirichletCondition(nodeList_0,fix,{1,1,1}); // fix u1 u2 u3

// Displace top nodes using Tensioner
if(DN.runningID()>=Tensioner::relaxSteps)
{
    const size_t nodeList_1=finiteElement().template createNodeList<Tensioner>();
    Tensioner tt(finiteElement());
    tt.updateDisplacement(DN);
    
    if (Tensioner::apply_tension)
    {
        displacement().addDirichletCondition(nodeList_1,tt,{0,0,1}); // prescribe only u3
    }
    else
    {
        std::cout<<"Tensioner IS NOT APPLYING ANY LOADS"<<std::endl;
    }
}
        
// Correct for Dislocation displacement
modifyDirichletConditions(DN);

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());