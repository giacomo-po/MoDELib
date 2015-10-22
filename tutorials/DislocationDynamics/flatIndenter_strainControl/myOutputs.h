/******************************************************************************/
// COMPUTATION OF DISPLACEMENT
double disp = 0.0;

// Obtain a list of nodes under the punch
const size_t nodeListID=DN.shared.bvpSolver.finiteElement().template createNodeList<FlatPunch>();

// Sum FEM displacement of those nodes
for(auto node : DN.shared.bvpSolver.finiteElement().nodeList(nodeListID))
{
    const Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().template segment<dim>(dim*node->gID);
    disp += nodeDisp(2);
}

// Compute dislocation displacement for those nodes
typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
typedef typename FieldPointType::DisplacementField DisplacementField;
std::deque<FieldPointType> fieldPoints; // the container of field points
for (auto node : DN.shared.bvpSolver.finiteElement().nodeList(nodeListID)) // range-based for loop (C++11)
{
    fieldPoints.emplace_back(*node);
}
DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);

// Sum dislocation displacement to disp
for(auto node : fieldPoints)
{
    Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>(); // line integral part of displacement
    
    if (DN.shared.use_virtualSegments) // solid-angle jump of virtual segments
    {
        for(const auto& segment : DN.links())
        {
            segment.second.addToSolidAngleJump(node.P,node.S,nodeDisp);
        }
    }
        
    disp += nodeDisp(2);
}

double avgdisp = disp/DN.shared.bvpSolver.finiteElement().nodeList(nodeListID).size();


/******************************************************************************/
// COMPUTATION OF LOAD
auto boundaryUnderPunch = FlatPunch::boundary(DN.shared.bvpSolver.finiteElement());
Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::bvpTraction); // integrate the bvp correction
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::ddTraction,DN); // integrate the dd traction

/******************************************************************************/
// OUTPUT TO F_0.txt
f_file<<avgdisp<<"  "<<force(2);



