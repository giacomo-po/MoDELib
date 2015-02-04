
//const double zSurf = 1.0e3;
//const double cubelength = 1000.0;
//const double V = cubelength*cubelength*cubelength;

/******************************************************************************/
// COMPUTATION OF DISPLACEMENT
double disp = 0.0;

// Obtain a list of nodes under the punch
auto nodeList=DN.shared.bvpSolver.finiteElement().template getNodeList<FlatPunch>();

// Sum FEM displacement of those nodes
for(auto node : nodeList)
{
    Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().template segment<dim>(dim*node->gID);
    disp += nodeDisp(2);
}

// Compute dislocation displacement for those nodes
typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
typedef typename FieldPointType::DisplacementField DisplacementField;
std::deque<FieldPointType> fieldPoints; // the container of field points
for (auto node : nodeList) // range-based for loop (C++11)
{
    // Compute S vector
    Eigen::Matrix<double,dim,1> s(Eigen::Matrix<double,dim,1>::Zero());
    for(auto ele : *node)
    {
        const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(node->P0));
        //std::map<double,int> baryIDmap;
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

// Sum dislocation displacement to disp
for(auto node : fieldPoints)
{
    Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>();
    disp += nodeDisp(2);
}

double avgdisp = disp/nodeList.size();


/******************************************************************************/
// COMPUTATION OF LOAD
auto boundaryUnderPunch = FlatPunch::boundary(DN.shared.bvpSolver.finiteElement());
Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::bvpTraction); // integrate the bvp correction
boundaryUnderPunch.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::ddTraction,DN); // integrate the dd traction


/******************************************************************************/
// PLASTIC STRAIN RATE
Eigen::Matrix<double,dim,dim> pSR(DN.plasticStrainRate());

/******************************************************************************/
// OUTPUT TO F_0.txt
UniqueOutputFile<'F'> f_file;
f_file<< DN.runningID()<<" "<<DN.get_dt()<<"  "<<DN.networkLength()<<" "<<avgdisp<<"  "<<force(2)<<"  "<<pSR(0,0)<<"  "<<pSR(0,1)<<"  "<<pSR(0,2)<<"  "<<pSR(1,1)<<"  "<<pSR(1,2) <<"  "<<pSR(2,2)<<"  "<<std::endl;



