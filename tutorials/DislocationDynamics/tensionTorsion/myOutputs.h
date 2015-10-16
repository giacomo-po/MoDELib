double avgDispZ = 0.0;
double avgRotZ =  0.0;
Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
Eigen::Matrix<double,3,1> torque(Eigen::Matrix<double,3,1>::Zero());

if(DN.shared.use_bvp)
{
    /******************************************************************************/
    // COMPUTATION OF DISPLACEMENT
    double dispZ = 0.0;
    double  rotZ = 0.0;
    
    TensionTorsioner tt(DN.shared.bvpSolver.finiteElement());
    
    // Obtain a list of nodes under the TensionTorsioner
    const size_t nodeListID=DN.shared.bvpSolver.finiteElement().template createNodeList<TensionTorsioner>();
    
    // Sum FEM displacement of those nodes
    for(auto node : DN.shared.bvpSolver.finiteElement().nodeList(nodeListID))
    {
        const Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().template segment<dim>(dim*node->gID);
        dispZ += nodeDisp(2);
        
        const Eigen::Matrix<double,3,1> v0=(node->P0-tt.pivot).normalized();
        const Eigen::Matrix<double,3,1> v1=((Eigen::Matrix<double,3,1>()<<(node->P0+nodeDisp-tt.pivot).template segment<2>(0),0.0).finished()).normalized();
        const double sinTheta=(v1-v1.dot(v0)*v0).dot(Eigen::Matrix<double,3,1>::UnitZ().cross(v0));
        if(sinTheta>=-1.0 && sinTheta<=1.0)
        {
            rotZ += std::asin(sinTheta);
        }
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
    
    // Sum dislocation displacement to dispZ and rotZ
    for(auto node : fieldPoints)
    {
        Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>();
        if (DN.shared.use_virtualSegments) // solid-angle jump of virtual segments
        {
            for(const auto& segment : DN.links())
            {
                segment.second.addToSolidAngleJump(node.P,node.S,nodeDisp);
            }
        }
        
        dispZ += nodeDisp(2);
        
        
        const Eigen::Matrix<double,3,1> v0=(node.P-tt.pivot).normalized();
        const Eigen::Matrix<double,3,1> v1=((Eigen::Matrix<double,3,1>()<<(node.P+nodeDisp-tt.pivot).template segment<2>(0),0.0).finished()).normalized();
        const double sinTheta=(v1-v1.dot(v0)*v0).dot(Eigen::Matrix<double,3,1>::UnitZ().cross(v0));
        if(sinTheta>=-1.0 && sinTheta<=1.0)
        {
            rotZ += std::asin(sinTheta);
        }
    }
    
    avgDispZ = dispZ/DN.shared.bvpSolver.finiteElement().nodeList(nodeListID).size();
    avgRotZ =  rotZ/DN.shared.bvpSolver.finiteElement().nodeList(nodeListID).size();
    
    /******************************************************************************/
    // COMPUTATION OF LOAD
    auto loadedBnd = TensionTorsioner::boundary(DN.shared.bvpSolver.finiteElement());

    typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
    loadedBnd.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::bvpTraction);    // integrate the bvp correction
    loadedBnd.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
    
    loadedBnd.integrate(&DN.shared.bvpSolver,torque,&BvpSolverType::bvpMoment,tt.pivot);   // integrate bvp moment about tt.pivot
    loadedBnd.integrate(&DN.shared.bvpSolver,torque,&BvpSolverType::ddMoment ,tt.pivot,DN); // integrate  dd moment about tt.pivot
}

/******************************************************************************/
// PLASTIC STRAIN RATE
Eigen::Matrix<double,dim,dim> pDR(DN.plasticDistortionRate());
const auto length=DN.networkLength();
/******************************************************************************/
// OUTPUT TO F_0.txt
UniqueOutputFile<'F'> f_file;
model::cout<<", F/F_0"<<std::flush;
f_file<< DN.runningID()<<" "<<DN.get_totalTime()<<" "<<DN.get_dt()<<"  "<<length.first<<" "<<length.second<<" "<<avgDispZ<<" "<<avgRotZ<<"  "<<force(2)<<" "<<torque(2)<<"  "<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<std::endl;




