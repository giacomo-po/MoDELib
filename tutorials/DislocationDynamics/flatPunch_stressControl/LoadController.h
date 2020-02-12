
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
    const double zFace;
    const double centerX;
    const double centerY;
    
    const int punchShape;
    const double punchSize;
    const double tol;
    
public:
    
    /**************************************/
    template <typename FiniteElementType>
    TopNodeSelector(const FiniteElementType& fe,
                    const int& _punchShape,
                    const double& _punchSize,
                    const double& tol_in=FLT_EPSILON) :
    /* init list */ zFace(fe.xMax()(2)),
    /* init list */ centerX((fe.xMax()(0)+fe.xMin()(0))*0.5),
    /* init list */ centerY((fe.xMax()(1)+fe.xMin()(1))*0.5),
    /* init list */ punchShape(_punchShape),
    /* init list */ punchSize(_punchSize),
    /* init list */ tol(tol_in)
    {
        
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {
        return this->operator()(node.P0);
    }
    
    /**************************************/
    bool operator()(const Eigen::Matrix<double,3,1>& P0) const
    {
        switch (punchShape)
        {
            case 0: // square punch
                return     std::fabs(P0(2)-zFace)<tol
                /*   */ && std::fabs(P0(0)-centerX)<(0.5*punchSize+tol)
                /*   */ && std::fabs(P0(1)-centerY)<(0.5*punchSize+tol);
                break;
                
            case 1: // circular punch
                return     std::fabs(P0(2)-zFace)<tol
                /*   */ && sqrt(std::pow(P0(0)-centerX,2)+std::pow(P0(1)-centerY,2))<(0.5*punchSize+tol);
                break;
                
            default:
                assert(0 && "Punch Shape not implemented.");
                return false;
                break;
        }
        
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
    
    int punchShape;
    
    double punchSize;
    
    //! Applied axial strain rate
    double sigmaDot;
    
    //! Initial value of the punch displacement
    double initialStress;
    
    // Number of DD steps before applying strain.
    int relaxSteps;
    
    double last_update_time;
    
    double deltaStress;
    
    bool apply_load;
    
    const size_t nodeList_bottom;
    
    size_t nodeList_top;
    IntegrationDomainType punchBnd;
    double punchArea;
    
    /**************************************************************************/
    LoadController(TrialFunctionType& u_in) :
    /* init list */ u(u_in),
    /* init list */ Lz(u.fe().xMax()(2)-u.fe().xMin()(2)),
    /* init list */ punchShape(0),
    /* init list */ punchSize(0.0),
    /* init list */ sigmaDot(1.0e-9),
    /* init list */ initialStress(0.0),
    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaStress(0.0),
    /* init list */ apply_load(true),
    /* init list */ nodeList_bottom(u.fe().template createNodeList<AtXmin<2>>()),
    /* init list */ nodeList_top(0),
    /* init list */ punchArea(0.0)
    {
        
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void init(const DislocationNetworkType& DN)
    {
        
        const long int runID=DN.runningID();
        const unsigned int userOutputColumn=DN.userOutputColumn();
        
        model::cout<<greenColor<<"Initializing Indenter at runID="<<runID<<defaultColor<<std::endl;
        
        
        model::EigenDataReader EDR;
        
        EDR.readScalarInFile("./loadInput.txt","punchShape",punchShape);
        EDR.readScalarInFile("./loadInput.txt","punchSize",punchSize);
        EDR.readScalarInFile("./loadInput.txt","sigmaDot",sigmaDot);
        EDR.readScalarInFile("./loadInput.txt","initialStress",initialStress);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",relaxSteps);
        EDR.readScalarInFile("./loadInput.txt","apply_load",apply_load);
        
        nodeList_top=u.fe().template createNodeList<TopNodeSelector>(punchShape,punchSize,0.01),
        punchBnd=boundaryUnderPunch(u.fe());
        punchArea=punchBnd.volume();
        model::cout<<greenColor<<"indenter area="<<punchArea<<defaultColor<<std::endl;
        
        // Read restart information
        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                
                initialStress=iter->second(userOutputColumn);
                model::cout<<"initialStress="<<initialStress<<std::endl;
                last_update_time=iter->second(0);
                model::cout<<"last_update_time="<<last_update_time<<std::endl;
                
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
        
        
        Eigen::VectorXd gv(Eigen::VectorXd::Zero(u.gSize()));
        
        if (apply_load)
        {
            auto f=make_constant((Eigen::Matrix<double,dim,1>()<<0.0,0.0,initialStress+deltaStress).finished());
//            auto dA=u.fe().template boundary<AtXmax<2>,3,GaussLegendre>();
            auto lWF=(test(u),f)*punchBnd;
            gv=lWF.globalVector();
        }
        else
        {
            std::cout<<"FlatPunch IS NOT APPLYING ANY LOADS"<<std::endl;
        }
        
        return gv;
        
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void addDirichletConditions(const DislocationNetworkType& DN)
    {
        
        // Fix bottom nodes
        Fix fix;
        u.addDirichletCondition(nodeList_bottom,fix,{1,1,1}); // fix u1 u2 u3
        
//        // Displace top nodes using operator()
//        if(DN.runningID()>=relaxSteps)
//        {
//            if (apply_load)
//            {
//                u.addDirichletCondition(nodeList_top,*this,{0,0,1}); // prescribe only u3
//            }
//            else
//            {
//                std::cout<<"Indenter IS NOT APPLYING ANY LOADS"<<std::endl;
//            }
//        }
    }
    
//    /**************************************/
//    template <typename NodeType,int dofPerNode>
//    Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType&,
//                                                   Eigen::Matrix<double,dofPerNode,1>& val) const
//    {
//        val=(Eigen::Matrix<double,3,1>()<<0.0,0.0,initialDisplacement+deltaDisplacement).finished();
//        return val;
//    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
//            deltaDisplacement += epsilonDot*deltaT*Lz;
            deltaStress += sigmaDot*deltaT;
            last_update_time += deltaT;
        }
    }
    
    /**************************************/
    template <typename FiniteElementType>
    IntegrationDomainType boundaryUnderPunch(const FiniteElementType& fe) const
    {
        
        TopNodeSelector tns(fe,punchShape,punchSize,0.01);
        
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
                        isFaceUnderPunch *= tns(*vertices[v]); // check if the current vertices satisfies operator()
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
        // OUTPUT contact triangles
        SequentialOutputFile<'I',1>::set_count(DN.runningID());
        //        SequentialOutputFile<'I',1>::set_increment(DN.outputFrequency);
        SequentialOutputFile<'I',1> indenterFile;
        
        int kk=0;
        for (const auto& pair : punchBnd)
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
        punchBnd.integrate(&DN.shared.bvpSolver,FEMforce,&BvpSolverType::bvpTraction);    // integrate the bvp correction
        
        Eigen::Matrix<double,3,1> DDforce(Eigen::Matrix<double,3,1>::Zero());
        punchBnd.integrate(&DN.shared.bvpSolver,DDforce,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
        
        std::stringstream os;
        os<<std::setprecision(15)<<std::scientific<<avgDispZ<<" "<<initialStress+deltaStress<<" "<<FEMforce(2)/punchArea<<" "<<DDforce(2)/punchArea<<" ";
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
