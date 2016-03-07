
#ifndef _Tensioner_h_
#define _Tensioner_h_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/Utilities/EigenDataReader.h>

using namespace model;

//template <typename DislocationNetworkType>
struct Tensioner
{
    
    //! z-coordinate of face at which Tensioner is applied
    const double zFace;
    
    //! Height of sample in z direction
    const double Lz;

//    const double R;
//    
//    const Eigen::Matrix<double,3,1> pivot;
//    
//    const double thetaDot;
    
    //! Height of sample in z direction
    const double geometricTol;
    
    //! Applied axial strain rate
    static double epsilonDot;

//    static double gammaDot;

    
    //! Initial value of the punch displacement
    static double initialDisplacement;
//    static double initialTwist_Rad;
    
    // Number of DD steps before applying strain.
    static int relaxSteps;
    
    static double last_update_time;
    
    static double deltaDisplacement;
//    static double deltaTheta;

    static bool apply_tension;
//    static bool apply_torsion;

    
    /**************************************/
    template <typename FiniteElementType>
    Tensioner(const FiniteElementType& fe) :
    //    /* init list */ DN(DN_in),
    /* init list */ zFace(fe.xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(fe.xMax()(2)-fe.xMin()(2)),
//    /* init list */ R((fe.xMax()(0)-fe.xMin()(0))*0.5),
//    /* init list */ pivot((fe.xMax()(0)+fe.xMin()(0))*0.5,(fe.xMax()(1)+fe.xMin()(1))*0.5,fe.xMax()(2)),
//    /* init list */ thetaDot(Lz/R*gammaDot),
    /* init list */ geometricTol(0.01)
    {
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {/*!@param[in] node a finite element node
      *\returns true if the node coordinate is under the Tensioner
      */

        return (std::fabs(node.P0(2)-zFace)<geometricTol);
        
    }
    
    /**************************************/
    template <typename NodeType,typename TrialFunctionType>
    Eigen::Matrix<double,3,1> operator()(const NodeType&, const TrialFunctionType&) const
    {
//        Eigen::Matrix<double,3,3> rot(Eigen::AngleAxisd(initialTwist_Rad+deltaTheta, Eigen::Vector3d::UnitZ()));
        return (Eigen::Matrix<double,3,1>()<<0.0,0.0,initialDisplacement+deltaDisplacement).finished();
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void updateDisplacement(const DislocationNetworkType& DN)
    {
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            deltaDisplacement += epsilonDot*deltaT*Lz;
//            deltaTheta += thetaDot*deltaT;
            last_update_time += deltaT;
        }
    }
    
    /**************************************/
    template <typename FiniteElementType>
    static IntegrationDomain<FiniteElementType,1,3,GaussLegendre> boundary(const FiniteElementType& fe)
    {
        IntegrationDomain<FiniteElementType,1,3,GaussLegendre> temp;
        
        // loop ever Elements
        for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
             /*                                                            */ eIter!=fe.elementEnd();
             /*                                                            */ eIter++)
        {
            if(eIter->second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= Tensioner(fe)(*vertices[v]); // check if the current vertices satisfies operator()
                    }
                    if(isExternalBoundaryFace)
                    {
                        temp.emplace_back(&eIter->second,boundaryFaces[f]);
                    }
                }
            }
        }
        
        return temp;
    }
    
    /**************************************/
    static void init(const long int& runID,
                     const unsigned int& userOutputColumn)
    {
        
        model::EigenDataReader EDR;
        
        EDR.readScalarInFile("./loadInput.txt","epsilonDot",epsilonDot);
//        EDR.readScalarInFile("./loadInput.txt","gammaDot",Tensioner::gammaDot);
        EDR.readScalarInFile("./loadInput.txt","initialDisplacement",initialDisplacement);
//        EDR.readScalarInFile("./loadInput.txt","initialTwist_Rad",Tensioner::initialTwist_Rad);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",Tensioner::relaxSteps);
        EDR.readScalarInFile("./loadInput.txt","apply_tension",Tensioner::apply_tension);
//        EDR.readScalarInFile("./loadInput.txt","apply_torsion",Tensioner::apply_torsion);

        
        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                model::cout<<"Initializing Tensioner at runID="<<runID<<std::endl;

                initialDisplacement=iter->second(userOutputColumn-1);
                model::cout<<"initialDisplacement="<<initialDisplacement<<std::endl;
//                initialTwist_Rad=iter->second(userOutputColumn);
//                model::cout<<"initialTwist_Rad="<<initialTwist_Rad<<std::endl;
                last_update_time=iter->second(0);
                model::cout<<"last_update_time="<<last_update_time<<std::endl;
                
            }
            else
            {
                //                assert(0 && "Tensioner::init runID not found inf F file");
            }
        }
        else
        {
            model::cout<<"Tensioner F/F_0.txt cannot be opened."<<std::endl;
        }
        
    }
    
};

// Static data
double Tensioner::epsilonDot=1.0e-9;
//double Tensioner::gammaDot=1.0e-9;
double Tensioner::initialDisplacement=0.0; // 0.0;
//double Tensioner::initialTwist_Rad=0.0;
int Tensioner::relaxSteps=0;
double Tensioner::last_update_time=0.0;
double Tensioner::deltaDisplacement=0.0;
//double Tensioner::deltaTheta=0.0;
bool Tensioner::apply_tension=true;
//bool Tensioner::apply_torsion=true;
//double Tensioner::target_VonMisesStrain=0.01;

#endif

