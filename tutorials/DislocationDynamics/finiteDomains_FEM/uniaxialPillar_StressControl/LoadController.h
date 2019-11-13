
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
    double sigmaDot;
    
    //! Initial value of the punch displacement
    double initialStress;
    
    // Number of DD steps before applying strain.
    int relaxSteps;
    
    double last_update_time;
    
    double deltaStress;
    
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
    /* init list */ sigmaDot(1.0e-9),
    /* init list */ initialStress(0.0),
    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaStress(0.0),
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
        
        sigmaDot=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<double>("sigmaDot",true);
        initialStress=TextFileParser("./inputFiles/loadControllerInput.txt").readScalar<double>("initialStress",true);
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
                    
                    initialStress=vReader(runID,"stress [mu]");//iter->second(userOutputColumn-1);
                    model::cout<<"initialStress="<<initialStress<<std::endl;
                    last_update_time=vReader(runID,"time [b/cs]");
                    
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
    Eigen::VectorXd globalVector(const DislocationNetworkType& DN) const
    {/*!\returns the zero-vector, since no Neumann boundary-conditions are applied
      */
        
        Eigen::VectorXd gv(Eigen::VectorXd::Zero(u.gSize()));
        if(enable)
        {
            if(DN.simulationParameters.runID>=relaxSteps)
            {
                auto f=make_constant((Eigen::Matrix<double,dim,1>()<<0.0,0.0,initialStress+deltaStress).finished());
                auto dA=u.fe().template boundary<AtXmax<2>,3,GaussLegendre>();
                auto lWF=(test(u),f)*dA;
                gv=lWF.globalVector();
            }
        }
        else
        {
            std::cout<<"LoadController is DISABLED"<<std::endl;
        }
        
        return gv;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void addDirichletConditions(const DislocationNetworkType&)
    {
        // Fix bottom nodes
        Fix fix;
        u.addDirichletCondition(nodeList_bottom,fix,{1,1,1}); // fix u1 u2 u3
    }
    
    
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
        const double deltaT = DN.simulationParameters.totalTime - last_update_time;
        last_update_time += deltaT;
        
        if(DN.simulationParameters.runID>=relaxSteps)
        {
            deltaStress += sigmaDot*deltaT;
        }
    }
    
    /**************************************************************************/
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
        double dispZ = 0.0;
        // Sum FEM displacement of nodes in nodeList_top
        for(const auto& node : u.fe().nodeList(nodeList_top))
        {
            const Eigen::Matrix<double,dim,1> nodeDisp = DN.bvpSolver->displacement().dofs(*node);
            dispZ += nodeDisp(2);
        }
        
        // Compute dislocation displacement for nodes in nodeList_top
        std::vector<FEMnodeEvaluation<typename DislocationNetworkType::ElementType,dim,1>> fieldPoints; // the container of field points
        fieldPoints.reserve(DN.bvpSolver->finiteElement().nodeList(nodeList_top).size());
        for (const auto& node : u.fe().nodeList(nodeList_top)) // range-based for loop (C++11)
        {
            fieldPoints.emplace_back(node->gID,node->P0);
        }
        DN.displacement(fieldPoints);
        
        // Sum dislocation displacement to dispZ
        for(const auto& node : fieldPoints)
        {
            dispZ += node(2);
        }
        
        // Average displacement
        const double avgdispZ = dispZ/u.fe().nodeList(nodeList_top).size();
        
        // COMPUTATION OF LOAD
        Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        loadedBnd.integrate(DN.bvpSolver.get(),force,&BvpSolverType::bvpTraction);    // integrate the bvp correction
        loadedBnd.integrate(DN.bvpSolver.get(),force,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
        
        f_file<<avgdispZ<<" "<<force(2)/topArea<<" ";
        
        if(runID==0)
        {
            F_labels<<"displacement [b]\n";
            F_labels<<"stress [mu]\n";
        }
        
        
    }
    
};

#endif

