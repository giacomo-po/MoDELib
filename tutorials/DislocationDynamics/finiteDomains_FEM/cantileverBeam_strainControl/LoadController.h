
#ifndef _LoadController_h_
#define _LoadController_h_

#include <iostream>
#include <sstream>      // std::stringstream
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <TerminalColors.h>
#include <MPIcout.h>
#include <IntegrationDomain.h>
#include <AtXmin.h>
#include <AtXmax.h>
#include <Fix.h>
#include <TextFileParser.h>
#include <IDreader.h>

using namespace model;

template <typename TrialFunctionType>
struct LoadController
{
    typedef typename TrialFunctionType::FiniteElementType FiniteElementType;
   
    static constexpr int dim=TrialFunctionType::dim;
    
    TrialFunctionType& u;
        
    //! Height of sample in z direction
    const double Lz;
    
    //! Height of sample in z direction
    const double geometricTol;
    
    //! Applied axial strain rate
    double epsilonDot;
    
    //! Initial value of the punch displacement
    double initialDisplacement;
    
    // Number of DD steps before applying strain.
    int relaxSteps;
    
    double last_update_time;
    
    double deltaDisplacement;
    
    bool enable;
    
    const size_t nodeList_bottom;
    
    const size_t nodeList_top;
    
    const IntegrationDomain<FiniteElementType,1,3,GaussLegendre> loadedBnd;
    
    const double topArea;
    
    /**************************************************************************/
    LoadController(TrialFunctionType& u_in) :
    /* init list */ u(u_in),
    //    /* init list */ zFace(u.fe().xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(u.fe().xMax()(2)-u.fe().xMin()(2)),
    /* init list */ geometricTol(0.01),
    /* init list */ epsilonDot(1.0e-9),
    /* init list */ initialDisplacement(0.0),
    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaDisplacement(0.0),
    /* init list */ enable(true),
    /* init list */ nodeList_bottom(u.fe().template createNodeList<AtXmin<2>>()),
    /* init list */ nodeList_top(u.fe().template createNodeList<AtXmax<2>>()),
    /* init list */ loadedBnd(topBoundary(u.fe())),
    /* init list */ topArea(loadedBnd.volume())
    {
        std::cout<<"LoadController: topArea="<<topArea<<std::endl;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void init(const DislocationNetworkType& DN)
    {
        
        const long int runID=DN.simulationParameters.runID;
        model::cout<<greenColor<<"Initializing LoadController at runID="<<runID<<defaultColor<<std::endl;
        
        epsilonDot=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<double>("epsilonDot",true);
        initialDisplacement=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<double>("initialDisplacement",true);
        relaxSteps=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<int>("relaxSteps",true);
        enable=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<int>("enable",true);
        
        if(runID>0)
        {// a restart
            IDreader<'F',1,200,double> vReader;
            vReader.readLabelsFile("F/F_labels.txt");
            if (vReader.isGood(0,true))
            {
                vReader.read(0,true);
                const auto iter=vReader.find(runID);
                if (iter!=vReader.end())
                {
                    
                    initialDisplacement=vReader(runID,"displacement [b]");
                    model::cout<<"initialDisplacement="<<initialDisplacement<<std::endl;

                    last_update_time=vReader(runID,"time [b/cs]");
                    model::cout<<"last_update_time="<<last_update_time<<std::endl;
                }
                else
                {
                    model::cout<<"LoadController::init runID="<<runID<<" not found in F file. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                model::cout<<"LoadController: F/F_0.txt cannot be opened."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    Eigen::VectorXd globalVector(const DislocationNetworkType&) const
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
        
        // Displace top nodes using operator()

            //            const size_t nodeList_top=u.fe().template createNodeList<AtXmax<2>>();
            //            Tensioner tt(finiteElement());
            
            if (enable)
            {
                if(DN.simulationParameters.runID>=relaxSteps)
                {
                    u.addDirichletCondition(nodeList_top,*this,{1,0,0}); // prescribe only u3
                }
            }
            else
            {
                std::cout<<"LoadController is DISABLED"<<std::endl;
            }
        
    }
    
    /**************************************/
    template <typename NodeType,int dofPerNode>
    Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType&,
                                                   Eigen::Matrix<double,dofPerNode,1>& val) const
    {
        val=(Eigen::Matrix<double,3,1>()<<initialDisplacement+deltaDisplacement,0.0,0.0).finished();
        return val;
    }
        
    /**************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
        const double deltaT = DN.simulationParameters.totalTime - last_update_time;
        last_update_time += deltaT;

        if(DN.simulationParameters.runID>=relaxSteps)
        {
            deltaDisplacement += epsilonDot*deltaT*Lz;
        }
    }
    
    /**************************************/
    template <typename FiniteElementType>
    IntegrationDomain<FiniteElementType,1,3,GaussLegendre> topBoundary(const FiniteElementType& fe) const
    {
        IntegrationDomain<FiniteElementType,1,3,GaussLegendre> temp;
        
        for(const auto& eIter : fe.elements())
        {
            if(eIter.second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter.second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= (std::fabs(vertices[v]->P0(2)-u.fe().xMax()(2))<geometricTol); // check if the current vertices satisfies operator()
                    }
                    if(isExternalBoundaryFace)
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
    void output(const DislocationNetworkType& DN,
                       const long int& runID,
                       UniqueOutputFile<'F'>& f_file,
                       std::ofstream& F_labels) const
    {
        
        // COMPUTATION OF DISPLACEMENT
        double dispX = 0.0;
        // Sum FEM displacement of nodes in nodeList_top
        for(const auto& node : u.fe().nodeList(nodeList_top))
        {
            const Eigen::Matrix<double,dim,1> nodeDisp = DN.bvpSolver->displacement().dofs(*node);
            dispX += nodeDisp(0);
        }
        
        // Compute dislocation displacement for nodes in nodeList_top
        std::vector<FEMnodeEvaluation<typename DislocationNetworkType::ElementType,dim,1>> fieldPoints; // the container of field points
        fieldPoints.reserve(DN.bvpSolver->finiteElement().nodeList(nodeList_top).size());
        for (const auto& node : u.fe().nodeList(nodeList_top)) // range-based for loop (C++11)
        {
            fieldPoints.emplace_back(node->gID,node->P0);
        }
        DN.displacement(fieldPoints);
        
        // Sum dislocation displacement to dispX
        for(const auto& node : fieldPoints)
        {
            dispX += node(0);
        }
        
        // Average displacement
        const double avgdispX = dispX/u.fe().nodeList(nodeList_top).size();
        
        // COMPUTATION OF LOAD
        Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        loadedBnd.integrate(DN.bvpSolver.get(),force,&BvpSolverType::bvpTraction);    // integrate the bvp correction
        loadedBnd.integrate(DN.bvpSolver.get(),force,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
        
        f_file<<avgdispX<<" "<<force(0)<<" ";
        
        if(runID==0)
        {
            F_labels<<"displacement [b]\n";
            F_labels<<"force [mu b^2]\n";
        }
        

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
