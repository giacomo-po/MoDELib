
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

//#include <model/Utilities/EigenDataReader.h>


using namespace model;

template <typename TrialFunctionType>
struct LoadController
{
    typedef typename TrialFunctionType::FiniteElementType FiniteElementType;
    
    static constexpr int dim=TrialFunctionType::dim;
    
    TrialFunctionType& u;
    
    
    //! Height of sample in z direction
    const double Lz;
    
    const Eigen::Matrix<double,3,1> pivot;
    
    //! Rate of rotation
    double thetaDot;
    
    //! Tolerance used to detect nodes on zFace
    const double geometricTol;
    
    //! Initial value of rotation
    double initialTwist_Rad;
    
    // Number of DD steps before applying strain.
    int relaxSteps;
    
    double last_update_time;
    
    double deltaTheta;
    
    bool apply_torsion;
    
    const size_t nodeList_bottom;
    
    const size_t nodeList_top;
    
    const IntegrationDomain<FiniteElementType,1,3,GaussLegendre> loadedBnd;
    
    //    const double topArea;
    
    /**************************************************************************/
    LoadController(TrialFunctionType& u_in) :
    /* init list */ u(u_in),
    /* init list */ Lz(u.fe().xMax()(2)-u.fe().xMin()(2)),
    /* init list */ pivot((u.fe().xMax()(0)+u.fe().xMin()(0))*0.5,(u.fe().xMax()(1)+u.fe().xMin()(1))*0.5,u.fe().xMax()(2)),
    /* init list */ thetaDot(1.0e-12),
    /* init list */ geometricTol(0.01),
    /* init list */ initialTwist_Rad(0.0),
    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaTheta(0.0),
    /* init list */ apply_torsion(true),
    /* init list */ nodeList_bottom(u.fe().template createNodeList<AtXmin<2>>()),
    /* init list */ nodeList_top(u.fe().template createNodeList<AtXmax<2>>()),
    /* init list */ loadedBnd(topBoundary(u.fe()))
    //    /* init list */ topArea(loadedBnd.volume())
    {
        //        std::cout<<"LoadController: topArea="<<topArea<<std::endl;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void init(const DislocationNetworkType& DN)
    {
        const long int runID=DN.simulationParameters.runID;
        model::cout<<greenColor<<"Initializing LoadController at runID="<<runID<<defaultColor<<std::endl;

//        if (runID%DN.bvpSolver->stepsBetweenBVPupdates)
//        {
//            model::cout<<"        LoadController can only be restarted with runID being a multiple of "<<DN.bvpSolver->stepsBetweenBVPupdates<<". EXITING "<<std::endl;
//            exit(EXIT_FAILURE);
//        }
        
        thetaDot=TextFileParser("./loadInput.txt").readScalar<double>("thetaDot",true);
        initialTwist_Rad=TextFileParser("./loadInput.txt").readScalar<double>("initialTwist_Rad",true);
        relaxSteps=TextFileParser("./loadInput.txt").readScalar<int>("relaxSteps",true);
        apply_torsion=TextFileParser("./loadInput.txt").readScalar<int>("apply_torsion",true);
        
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
                    model::cout<<"Initializing LoadController at runID="<<runID<<std::endl;
                    
                    initialTwist_Rad=vReader(runID,"twist angle [rad]");
                    //                initialTwist_Rad=iter->second(userOutputColumn-1);
                    model::cout<<"initialTwist_Rad="<<initialTwist_Rad<<std::endl;
                    
                    //                last_update_time=iter->second(0);
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
        
        // Displace top nodes using operator()
        if(DN.simulationParameters.runID>=relaxSteps)
        {
            if (apply_torsion)
            {
                u.addDirichletCondition(nodeList_top,*this,{1,1,0}); // prescribe only u3
            }
            else
            {
                std::cout<<"LoadController IS NOT APPLYING ANY LOADS"<<std::endl;
            }
        }
    }
    
    /**************************************/
    template <typename NodeType,int dofPerNode>
    Eigen::Matrix<double,dofPerNode,1>& operator()(const NodeType& node,
                                                   Eigen::Matrix<double,dofPerNode,1>& val) const
    {
        Eigen::Matrix<double,3,3> rot(Eigen::AngleAxisd(initialTwist_Rad+deltaTheta, Eigen::Vector3d::UnitZ()));
        val=(Eigen::Matrix<double,3,1>()<<(rot*(node.P0-pivot)-(node.P0-pivot)).template segment<2>(0),0.0).finished();
        return val;
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {/*!@param[in] DN the DislocationNetwork
      *
      * Computes the change in rotation angle based on the time allapsed
      * since last update.
      */
        if(DN.simulationParameters.runID>=relaxSteps)
        {
            const double deltaT = DN.simulationParameters.totalTime - last_update_time;
            deltaTheta += thetaDot*deltaT;
            last_update_time += deltaT;
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
        double  rotZ = 0.0;
        
        // Sum FEM displacement of those nodes
        for(const auto& node : DN.bvpSolver->finiteElement().nodeList(nodeList_top))
        {
            const Eigen::Matrix<double,dim,1> nodeDisp = DN.bvpSolver->displacement().dofs(*node);
            
            const Eigen::Matrix<double,dim,1> v0=(node->P0-pivot).normalized();
            const Eigen::Matrix<double,dim,1> v1=((Eigen::Matrix<double,3,1>()<<(node->P0+nodeDisp-pivot).template segment<2>(0),0.0).finished()).normalized();
            const double sinTheta=(v1-v1.dot(v0)*v0).dot(Eigen::Matrix<double,3,1>::UnitZ().cross(v0));
            if(sinTheta>=-1.0 && sinTheta<=1.0)
            {
                rotZ += std::asin(sinTheta);
            }
        }
        
        // Compute dislocation displacement for those nodes
        std::vector<FEMnodeEvaluation<typename DislocationNetworkType::ElementType,dim,1>> fieldPoints; // the container of field points
        fieldPoints.reserve(DN.bvpSolver->finiteElement().nodeList(nodeList_top).size());
        for (const auto& node : DN.bvpSolver->finiteElement().nodeList(nodeList_top)) // range-based for loop (C++11)
        {
            fieldPoints.emplace_back(node->gID,node->P0);
        }
        DN.displacement(fieldPoints);
        
        // Sum dislocation displacement to dispZ and rotZ
        for(const auto& node : fieldPoints)
        {
            
            const Eigen::Matrix<double,dim,1> v0=(node.P-pivot).normalized();
            const Eigen::Matrix<double,dim,1> v1=((Eigen::Matrix<double,3,1>()<<(node.P+node-pivot).template segment<2>(0),0.0).finished()).normalized();
            const double sinTheta=(v1-v1.dot(v0)*v0).dot(Eigen::Matrix<double,3,1>::UnitZ().cross(v0));
            if(sinTheta>=-1.0 && sinTheta<=1.0)
            {
                rotZ += std::asin(sinTheta);
            }
        }
        
        double avgRotZ =  rotZ/DN.bvpSolver->finiteElement().nodeList(nodeList_top).size();
        
        // COMPUTATION OF LOAD
        //        auto loadedBnd = topBoundary(DN.shared.bvpSolver->finiteElement());
        
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        Eigen::Matrix<double,dim,1> torque(Eigen::Matrix<double,dim,1>::Zero());
        loadedBnd.integrate(DN.bvpSolver.get(),torque,&BvpSolverType::bvpMoment,pivot);   // integrate bvp moment about tt.pivot
        loadedBnd.integrate(DN.bvpSolver.get(),torque,&BvpSolverType::ddMoment ,pivot,DN); // integrate  dd moment about tt.pivot
        
        
        
        //        std::stringstream os;
        f_file<<" "<<avgRotZ<<" "<<" "<<torque(2);
        //        return os.str();
        
        if(runID==0)
        {
            F_labels<<"twist angle [rad]\n";
            F_labels<<"torque [mu b^3]\n";
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
//        return (std::fabs(node.P0(2)-zFace)<geometricTol);
//
//    }
