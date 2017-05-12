
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
#include <model/FEM/Boundaries/MidPlane.h>
#include <model/FEM/Boundaries/BoundaryIntersection.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/FEM/BoundaryConditions/AxisRotation.h>
#include <model/Utilities/EigenDataReader.h>

using namespace model;

///**************************************************************************/
//template<int comp,int gradDir>
//struct LinearTraction
//{
//    
//    constexpr static int rows=3;
//    constexpr static int cols=1;
//    
//    
//    const double stressGradient;
//    const double refCoordinate;
//
//    /**********************************************************************/
//    LinearTraction(const double& stressGradient_in,
//                   const double& refCoordinate_in=0.0) :
//    stressGradient(stressGradient_in),
//    refCoordinate(refCoordinate_in)
//    {
//        
//    }
//    
//    /**********************************************************************/
//    template<typename ElementType, typename BaryType>
//    const Eigen::Matrix<double,3,1> operator() (const ElementType& ele, const BaryType& bary) const
//    {/*!@param[in] elem the element
//      * @param[in] bary the barycentric cooridinate
//      *\returns the constant c.
//      */
//        Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Zero());
//        temp(comp)=stressGradient*(ele.simplex.bary2pos(bary)(gradDir)-refCoordinate);
//        return temp;
//    }
//    
//    /**********************************************************************/
//    EvalExpression<LinearTraction<comp,gradDir>> eval() const
//    {
//        return EvalExpression<LinearTraction<comp,gradDir>>(*this);
//    }
//    
//};

/**************************************************************************/
template <typename TrialFunctionType>
struct LoadController
{
    
    typedef typename TrialFunctionType::FiniteElementType FiniteElementType;
    
    static constexpr int dim=TrialFunctionType::dim;
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    
    TrialFunctionType& u;
        
    //! Height of sample in z direction
    const double Lz;
    
    //! Height of sample in z direction
    const double geometricTol;
    
    //! Applied axial strain rate
    double thetaDot;
    
    //! Initial value of the punch displacement
    double initialRotation;
    
    // Number of DD steps before applying strain.
//    int relaxSteps;
    
    double last_update_time;
    
    double deltaTheta;
    
    bool apply_bending;
    
    const size_t nodeList_left;
    const size_t nodeList_right;
    const VectorDim Pleft;
    const VectorDim Pright;
    
    const size_t nodeList_top;
    
    const IntegrationDomain<FiniteElementType,1,3,GaussLegendre> loadedBnd;
    
    const double topArea;
    
    /**************************************************************************/
    LoadController(TrialFunctionType& u_in) :
    /* init list */ u(u_in),
    //    /* init list */ zFace(u.fe().xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(u.fe().xMax()(2)-u.fe().xMin()(2)),
    /* init list */ geometricTol(0.01),
    /* init list */ thetaDot(1.0e-9),
    /* init list */ initialRotation(0.0),
//    /* init list */ relaxSteps(0),
    /* init list */ last_update_time(0.0),
    /* init list */ deltaTheta(0.0),
    /* init list */ apply_bending(true),
//    /* init list */ nodeList_left(u.fe().template createNodeList<BoundaryIntersection<AtXmin<1>,MidPlane<2>>>()),
//    /* init list */ nodeList_right(u.fe().template createNodeList<BoundaryIntersection<AtXmax<1>,MidPlane<2>>>()),
    /* init list */ nodeList_left (u.fe().template createNodeList<AtXmin<1>>()),
    /* init list */ nodeList_right(u.fe().template createNodeList<AtXmax<1>>()),
    /* init list */ Pleft((VectorDim()<<u.fe().xMean(0),u.fe().xMin(1),u.fe().xMean(2)).finished()),
    /* init list */ Pright((VectorDim()<<u.fe().xMean(0),u.fe().xMax(1),u.fe().xMean(2)).finished()),
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
        
        const long int runID=DN.runningID();
        const unsigned int userOutputColumn=DN.userOutputColumn();
        
        model::cout<<greenColor<<"Initializing LoadController at runID="<<runID<<defaultColor<<std::endl;
        
        model::EigenDataReader EDR;
        
        EDR.readScalarInFile("./loadInput.txt","thetaDot",thetaDot);
        EDR.readScalarInFile("./loadInput.txt","initialRotation",initialRotation);
//        EDR.readScalarInFile("./loadInput.txt","relaxSteps",LoadController::relaxSteps);
        EDR.readScalarInFile("./loadInput.txt","apply_bending",LoadController::apply_bending);
        
        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                
                initialRotation=iter->second(userOutputColumn);
                model::cout<<"initialRotation="<<initialRotation<<std::endl;
                //                initialTwist_Rad=iter->second(userOutputColumn);
                //                model::cout<<"initialTwist_Rad="<<initialTwist_Rad<<std::endl;
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
        
        // Displace top nodes using operator()

            
//            if (apply_bending)
//            {
//                LinearTraction<1,2> lt0(-(initialRotation+deltaTheta)/(u.fe().xMax()(2)-u.fe().xMin()(2)));
//                auto dA0=u.fe().template boundary<AtXmin<1>,3,GaussLegendre>();
//                auto lWF0=(u.test(),lt0.eval())*dA0;
//
//                LinearTraction<1,2> lt1((initialRotation+deltaTheta)/(u.fe().xMax()(2)-u.fe().xMin()(2)));
//                auto dA1=u.fe().template boundary<AtXmax<1>,3,GaussLegendre>();
//                auto lWF1=(u.test(),lt1.eval())*dA1;
//                gv=lWF0.globalVector()+lWF1.globalVector();
//            }
//            else
//            {
//                std::cout<<"Tensioner IS NOT APPLYING ANY LOADS"<<std::endl;
//            }
        
        return gv;
    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void addDirichletConditions(const DislocationNetworkType&)
    {
        // Fix bottom nodes
//        Fix fix;
//        u.addDirichletCondition(nodeList_left,fix,{1,1,1}); // fix u1 u2 u3
//        u.addDirichletCondition(nodeList_right,fix,{0,0,1}); // fix u3
        
        const VectorDim axis(VectorDim::UnitX());
        AxisRotation arl(Pleft,axis,-(initialRotation+deltaTheta));
        u.addDirichletCondition(nodeList_left,arl,{1,1,1}); // rotate u1 u2 u3

        AxisRotation arr(Pright,axis,initialRotation+deltaTheta);
        u.addDirichletCondition(nodeList_right,arr,{1,1,1}); // rotate u1 u2 u3

    }
    

    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
//        if(DN.runningID()>=relaxSteps)
//        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            deltaTheta += thetaDot*deltaT;
            //            deltaTheta += thetaDot*deltaT;
            last_update_time += deltaT;
//        }
    }
    
    /**************************************************************************/
    template <typename FiniteElementType>
    IntegrationDomain<FiniteElementType,1,3,GaussLegendre> topBoundary(const FiniteElementType& fe) const
    {
        IntegrationDomain<FiniteElementType,1,3,GaussLegendre> temp;
        
        // loop ever Elements
        for (const auto& eIter : fe.elements())
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
    std::string output(const DislocationNetworkType& DN) const
    {
        
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
//        auto loadedBnd = topBoundary(u.fe);
        Eigen::Matrix<double,3,1> force(Eigen::Matrix<double,3,1>::Zero());
        typedef typename DislocationNetworkType::BvpSolverType BvpSolverType;
        loadedBnd.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::bvpTraction);    // integrate the bvp correction
        loadedBnd.integrate(&DN.shared.bvpSolver,force,&BvpSolverType::ddTraction,DN);  // integrate the dd traction
        
        std::stringstream os;
        os<<std::setprecision(15)<<std::scientific<<avgDispZ<<" "<<force(2)/topArea<<" ";
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

//    /**************************************/
//    template <typename NodeType>
//    Eigen::Matrix<double,3,1> operator()(const NodeType&, const TrialFunctionType&) const
//    {
//        return (Eigen::Matrix<double,3,1>()<<0.0,0.0,initialRotation+deltaTheta).finished();
//    }
