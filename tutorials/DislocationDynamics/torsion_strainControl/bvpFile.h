// Fix nodes at z=0 in u1 and u2
const size_t nodeList_0=finiteElement().template createNodeList<AtXmin<2>>();
Fix fix;
displacement().addDirichletCondition(nodeList_0,fix,{1,1,0}); // fix u1 u2

// Fix origin in u1, u2, u3
const size_t nodeList_1=finiteElement().template createNodeList<AtPoint>(Eigen::Vector3d::Zero());
displacement().addDirichletCondition(nodeList_1,fix,{1,1,1}); // fix u1 u2 u3

// Displace top nodes using Torsioner
if(DN.runningID()>=Torsioner::relaxSteps)
{
    const size_t nodeList_2=finiteElement().template createNodeList<Torsioner>();
    Torsioner tt(finiteElement());
    tt.updateDisplacement(DN);
    
    if(Torsioner::apply_torsion)
    {
        displacement().addDirichletCondition(nodeList_2,tt,{1,1,0}); // prescribe only u1 & u2
    }
    else
    {
        std::cout<<"Torsioner IS NOT APPLYING ANY LOADS"<<std::endl;
    }
}
        
// Correct for Dislocation displacement
modifyDirichletConditions(DN);

// CALL SOLVER, passing the negative dislocation traction vector and the solution guess
displacement()=solve(-dislocationTraction.globalVector(),displacement());