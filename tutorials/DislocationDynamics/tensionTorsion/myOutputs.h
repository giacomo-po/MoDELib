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
        Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().template segment<dim>(dim*node->gID);
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
        // Compute S vector
        Eigen::Matrix<double,dim,1> s(Eigen::Matrix<double,dim,1>::Zero());
        for(auto ele : *node)
        {
            const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(node->P0));
            for(int k=0;k<dim+1;++k)
            {
                if (std::fabs(bary(k))<FLT_EPSILON && ele->simplex.child(k).isBoundarySimplex())
                {
                    s += ele->simplex.nda.col(k).normalized();
                }
            }
        }
        const double sNorm(s.norm());
        assert(sNorm>0.0 && "s-vector has zero norm.");
        fieldPoints.emplace_back(*node,s/sNorm);
    }
    DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);
    
    // Sum dislocation displacement to dispZ and rotZ
    for(auto node : fieldPoints)
    {
        Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>();
        dispZ += nodeDisp(2);
        
        
        const Eigen::Matrix<double,3,1> v0=(node.P-tt.pivot).normalized();
        const Eigen::Matrix<double,3,1> v1=((Eigen::Matrix<double,3,1>()<<(node.P+nodeDisp-tt.pivot).template segment<2>(0),0.0).finished()).normalized();
        const double sinTheta=(v1-v1.dot(v0)*v0).dot(Eigen::Matrix<double,3,1>::UnitZ().cross(v0));
        if(sinTheta>=-1.0 && sinTheta<=1.0)
        {
            rotZ += std::asin(sinTheta);
        }
        
        //    const double dotp((node.P-tt.pivot).template segment<2>(0).normalized().dot((node.P+nodeDisp-tt.pivot).template segment<2>(0).normalized()));
        //    if (dotp>0.1 && dotp<=1.0)
        //    {
        //        rotZ += acos(dotp);
        //    }
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




