
#ifndef _LoadController_h_
#define _LoadController_h_

#include <iostream>
#include <sstream>      // std::stringstream
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <model/Utilities/TerminalColors.h>
#include <model/MPI/MPIcout.h>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/Utilities/EigenDataReader.h>

using namespace model;

/******************************************************************************/
/******************************************************************************/
struct TopNodeSelector
{
    const Eigen::VectorXd& U;
    const double zFace;
    const int radius;
    const int radiusSquared;
    const double depth;
    const Eigen::Vector3d C;
    const double tol;
    
public:
    
    /**************************************/
    template <typename FiniteElementType>
    TopNodeSelector(const FiniteElementType& fe,
                    const Eigen::VectorXd& u_in,
                    const int& _radius,
                    const double& _depth,
                    const double& tol_in=FLT_EPSILON) :
    U(u_in),
    /* init list */ zFace(fe.xMax()(2)),
    /* init list */ radius(_radius),
    /* init list */ radiusSquared(_radius*_radius),
    /* init list */ depth(_depth),
    /* init list */ C((Eigen::Vector3d()<<(fe.xMax()(0)+fe.xMin()(0))*0.5,(fe.xMax()(1)+fe.xMin()(1))*0.5,zFace+radius-depth).finished()),
    /* init list */ tol(tol_in)
    {
        
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {
        const size_t& gID(node.gID);
        return this->operator()(node.P0,U.template segment<3>(3*gID));
    }
    
    /**************************************/
    bool operator()(const Eigen::Matrix<double,3,1>& P0,const Eigen::Matrix<double,3,1>& u) const
    {
        return (P0+u-C).squaredNorm()<radiusSquared;
    }
    
};

/******************************************************************************/
/******************************************************************************/
template <typename TrialFunctionType>
struct LoadController
{
    typedef typename TrialFunctionType::FiniteElementType FiniteElementType;
    typedef typename FiniteElementType::ElementType ElementType;
    static constexpr int dim=TrialFunctionType::dim;
    typedef IntegrationDomain<FiniteElementType,1,3,GaussLegendre> IntegrationDomainType;
    
    TrialFunctionType& u;
    
    //! Height of sample in z direction
    const double Lz;
    
//    int punchShape;
    
    double indenterRadius;
    
    //! Applied axial strain rate
    double epsilonDot;
    
    //! Initial value of the punch displacement
    double initialDepth;
    
    // Number of DD steps before applying strain.
    int relaxSteps;
    
    double last_update_time;
    
    double deltaDepth;
    
    bool apply_load;
    
    const size_t nodeList_bottom;
    
    Eigen::Vector3d C;

    
//    size_t nodeList_top;
//    IntegrationDomainType indentedBnd;
//    double projectedArea;
    
    /**************************************************************************/
    LoadController(TrialFunctionType& u_in) :
    /* init list */ u(u_in),
    /* init list */ Lz(u.fe().xMax()(2)-u.fe().xMin()(2)),
//    /* init list */ punchShape(0),
    /* init list */ indenterRadius(0.0),
    /* init list */ epsilonDot(1.0e-9),
    /* init list */ initialDepth(0.0),
    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaDepth(0.0),
    /* init list */ apply_load(true),
    /* init list */ nodeList_bottom(u.fe().template createNodeList<AtXmin<2>>()),
    /* init list */ C((Eigen::Vector3d()<<(u.fe().xMax()(0)+u.fe().xMin()(0))*0.5,(u.fe().xMax()(1)+u.fe().xMin()(1))*0.5,u.fe().xMax()(2)+indenterRadius-(initialDepth+deltaDepth)).finished())
//    /* init list */ nodeList_top(0),
//    /* init list */ projectedArea(0.0)
    {
        
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void init(const DislocationNetworkType& DN)
    {
        
        const long int runID=DN.runningID();
        const unsigned int userOutputColumn=DN.userOutputColumn();
        
        model::cout<<greenColor<<"Initializing SphericalIndenter at runID="<<runID<<defaultColor<<std::endl;
        
        
        model::EigenDataReader EDR;
        
        EDR.readScalarInFile("./loadInput.txt","indenterRadius",LoadController::indenterRadius);
//        EDR.readScalarInFile("./loadInput.txt","intenderRadius",LoadController::punchSize);
        EDR.readScalarInFile("./loadInput.txt","epsilonDot",epsilonDot);
        EDR.readScalarInFile("./loadInput.txt","initialDepth",initialDepth);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",LoadController::relaxSteps);
        EDR.readScalarInFile("./loadInput.txt","apply_load",LoadController::apply_load);
        
//        nodeList_top=u.fe().template createNodeList<TopNodeSelector>(punchShape,punchSize,0.01),
//        indentedBnd=boundaryUnderIndenter(u.fe());
//        projectedArea=indentedBnd.volume();
//        model::cout<<greenColor<<"indenter area="<<projectedArea<<defaultColor<<std::endl;
        
        // Read restart information
        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                
                initialDepth=iter->second(userOutputColumn-1);
                model::cout<<"initialDepth="<<initialDepth<<std::endl;
                last_update_time=iter->second(0);
                model::cout<<"last_update_time="<<last_update_time<<std::endl;
                
                C<<(u.fe().xMax()(0)+u.fe().xMin()(0))*0.5,(u.fe().xMax()(1)+u.fe().xMin()(1))*0.5,u.fe().xMax()(2)+indenterRadius-(initialDepth+deltaDepth);

                
            }
            else
            {
                //                assert(0 && "LoadController::init runID not found inf F file");
            }
        }
        else
        {
            model::cout<<"LoadController: F/F_0.txt cannot be opened."<<std::endl;
        }
        
    }
    
    /**************************************************************************/
    Eigen::VectorXd globalVector() const
    {/*!\returns the zero-vector, since no Neumann boundary-conditions are applied
      */
        return Eigen::VectorXd::Zero(u.gSize());
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void addDirichletConditions(const DislocationNetworkType& DN)
    {
        
        // Fix bottom nodes
        Fix fix;
        u.addDirichletCondition(nodeList_bottom,fix,{1,1,1}); // fix u1 u2 u3
        
        
        size_t nodeList_top=u.fe().template createNodeList<TopNodeSelector>(u.dofVector(),indenterRadius,initialDepth+deltaDepth,0.01);
        IntegrationDomainType indentedBnd=boundaryUnderIndenter(u.fe());
        double projectedArea=indentedBnd.volume();
        model::cout<<greenColor<<"C="<<C.transpose()<<std::endl;
        model::cout<<greenColor<<"indenter area="<<projectedArea<<defaultColor<<std::endl;

        
        // Displace top nodes using operator()
        if(DN.runningID()>=relaxSteps)
        {
            if (apply_load)
            {
                u.addDirichletCondition(nodeList_top,*this,{1,1,1}); // prescribe all components of displacement
            }
            else
            {
                std::cout<<"Indenter IS NOT APPLYING ANY LOADS"<<std::endl;
            }
        }
    }
    
    /**************************************/
    template <typename NodeType,int dofPerNode>
    Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType& node,
                                                   Eigen::Matrix<double,dofPerNode,1>& val) const
    {
        
        //value of new displacement C-P0+radius*d;
        const size_t& gID(node.gID);
        const Eigen::Vector3d d=(node.P0+u.dofs(node)-C).normalized();
        //        const Eigen::Vector3d d=(P0+u-C).normalized();

        val=C-node.P0+indenterRadius*d;
        return val;
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            deltaDepth += epsilonDot*deltaT*Lz;
            last_update_time += deltaT;
            C<<(u.fe().xMax()(0)+u.fe().xMin()(0))*0.5,(u.fe().xMax()(1)+u.fe().xMin()(1))*0.5,u.fe().xMax()(2)+indenterRadius-(initialDepth+deltaDepth);
        }
    }
    
    /**************************************/
    template <typename FiniteElementType>
    IntegrationDomainType boundaryUnderIndenter(const FiniteElementType& fe) const
    {
        
        TopNodeSelector tns(u.fe(),u.dofVector(),indenterRadius,initialDepth+deltaDepth,0.01);
        
        IntegrationDomainType temp;
        
        // loop ever Elements
        //        for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
        //             /*                                                            */ eIter!=fe.elementEnd();
        //             /*                                                            */ eIter++)
        for(const auto& eIter : fe.elements())
        {
            if(eIter.second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop over boundary faces of the current element
                {
                    bool isFaceUnderPunch(true);
                    //                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter.second.simplex.child(boundaryFaces[f]).vertices();
                    const auto vertices=eIter.second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isFaceUnderPunch *= tns(*u.fe().mesh2femIDmap()[vertices[v]->xID(0)]); // check if the current vertices satisfies operator()
                    }
                    if(isFaceUnderPunch)
                    {
                        temp.emplace_back(&eIter.second,boundaryFaces[f]);
                    }
                }
            }
        }
        
        return temp;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    std::string output(const DislocationNetworkType& DN) const
    {
        size_t nodeList_top=u.fe().template createNodeList<TopNodeSelector>(u.dofVector(),indenterRadius,initialDepth+deltaDepth,0.01);
        IntegrationDomainType indentedBnd=boundaryUnderIndenter(u.fe());
        double projectedArea=indentedBnd.volume();

        
        // OUTPUT contact triangles
        SequentialOutputFile<'I',1>::set_count(DN.runningID());
        //        SequentialOutputFile<'I',1>::set_increment(DN.outputFrequency);
        SequentialOutputFile<'I',1> indenterFile;
        
        int kk=0;
        for (const auto& pair : indentedBnd)
        {
            const ElementType& ele(*pair.first);
            const int& boundaryFaceID(pair.second);
            indenterFile<<kk<<" "<<ele.simplex.child(boundaryFaceID).xID<<"\n";
            ++kk;
        }
        
        
        // COMPUTATION OF DISPLACEMENT
        double dispZ = 0.0;
        // Sum FEM displacement of nodes in nodeList_top
        for(auto node : u.fe().nodeList(nodeList_top))
        {
            const Eigen::Matrix<double,dim,1> nodeDisp = DN.shared.bvpSolver.displacement().dofs(*node);
            dispZ += nodeDisp(2);
        }
        
        // Compute dislocation displacement for nodes in nodeList_top
        typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
        typedef typename FieldPointType::DisplacementField DisplacementField;
        std::deque<FieldPointType> fieldPoints; // the container of field points
        for (auto node : u.fe().nodeList(nodeList_top)) // range-based for loop (C++11)
        {
            fieldPoints.emplace_back(*node);
        }
        DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);
        
        // Sum dislocation displacement to dispZ
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
        }
        
        // Average displacement
        const double avgDispZ = dispZ/u.fe().nodeList(nodeList_top).size();
        
        // COMPUTATION OF LOAD
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        
        Eigen::Matrix<double,3,1> FEMforce(Eigen::Matrix<double,3,1>::Zero());
        indentedBnd.integrate(&DN.shared.bvpSolver,FEMforce,&BvpSolverType::bvpTraction);    // integrate the bvp correction
        
        Eigen::Matrix<double,3,1> DDforce(Eigen::Matrix<double,3,1>::Zero());
        indentedBnd.integrate(&DN.shared.bvpSolver,DDforce,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
        
        std::stringstream os;
        os<<std::setprecision(15)<<std::scientific<<avgDispZ<<" "<<FEMforce(2)/projectedArea<<" "<<DDforce(2)/projectedArea<<" ";
        return os.str();
    }
    
    
    
};

#endif


//    /**************************************/
//    template <typename NodeType>
//    bool operator()(const NodeType& node) const
//    {/*!@param[in] node a finite element node
//      *\returns true if the node coordinate is under the LoadController
//      */
//
//        return (std::fabs(node.P0(2)-u.fe().xMax()(2))<geometricTol);
//
//    }
